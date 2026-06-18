library(qs2)
library(openxlsx)
library(tidyverse)
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)


# \\\\\\\\\\\\\\\\\\\\\\\\\\ #
# \\\\\\\ PREPARATION \\\\\\ #
# \\\\\\\\\\\\\\\\\\\\\\\\\\ #


# Constants ----
prefix <- '20260616_default_protocol'
out.dir <- 'pathway_enrichment_analysis'
de.file <- 'default_results/20260209_default_protocol_all_contrasts_differential_expression_results_edgeR4.qs2'
counts.file <- 'default_results/dge_matrix.qs2'
design.file <- 'default_results/design_matrix.qs2'

# Open results ----
df <- qs_read(de.file)
counts <- qs_read(counts.file)
design_matrix <- qs_read(design.file)

# MSigDB datasets ----
H_genes <- read.gmt('MSigDB/h.all.v2026.1.Hs.entrez.gmt')
C2_genes <- read.gmt('MSigDB/c2.all.v2026.1.Hs.entrez.gmt') # Includes KEGG and Reactome
C5_genes <- read.gmt('MSigDB/c5.all.v2026.1.Hs.entrez.gmt')
C6_genes <- read.gmt('MSigDB/c6.all.v2026.1.Hs.entrez.gmt')

# Camera index ----
prepare_camera_idx <- function(gmt.genes){
  tmp <- gmt.genes %>% 
    dplyr::group_by(term) %>%
    dplyr::summarise(genes = list(gene), .groups = "drop") %>%
    tibble::deframe()
  
  ids2indices(gene.sets = tmp,
                    identifiers = counts$genes$entrezid)
}
idxH <- prepare_camera_idx(H_genes)
idxC2 <- prepare_camera_idx(C2_genes)
idxC5 <- prepare_camera_idx(C5_genes)
idxC6 <- prepare_camera_idx(C6_genes)


# \\\\\\\\\\\\\\\\\\\\\\\\\\ #
# \\\\\\\ ENRICHMENT \\\\\\\ #
# \\\\\\\\\\\\\\\\\\\\\\\\\\ #


# Enrichment analysis ----
safe_enrichGO <- possibly(\(gene_list, direction, ont) {
    tmp <- enrichGO(
      gene = gene_list, 
      OrgDb = org.Hs.eg.db, 
      keyType = "ENTREZID", 
      ont = ont,
      pAdjustMethod = "BH", 
      pvalueCutoff = 1, 
      readable = TRUE,
      minGSSize = 10,
      maxGSSize = 500,
      pool = TRUE
    )
    tmp <- enrichplot::pairwise_termsim(tmp)
    suppressWarnings(clusterProfiler::simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min)) %>% 
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }, otherwise = NULL)
safe_enricher <- possibly(\(gene_list, direction, term) {
  enricher(gene = gene_list,
           TERM2GENE = term,
           pAdjustMethod = "BH",
           pvalueCutoff = 1,
           minGSSize = 10,
           maxGSSize = 500) %>% 
    setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
}, otherwise = NULL)
enrich2tb <- function(enrichResult, direction) {
  ont <- ifelse(enrichResult@ontology == "unknown", NA_character_, enrichResult@ontology)
  enrichResult %>% 
    as_tibble() %>%
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::mutate(GO = ont,
                  GeneDirection = direction) %>%
    dplyr::rename(pvalue_BH = p.adjust) %>%
    dplyr::select(ID, GO, Description, GeneDirection,GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, pvalue_BH, Count, geneID)
}
enrichment_function <- function(de_results, contrast, mode=c("normal", "adjusted")){
  mode <- match.arg(mode)
  if(!(dir.exists(glue::glue("{out.dir}/enrichment_analysis")))){
    dir.create(glue::glue("{out.dir}/enrichment_analysis"),
               recursive = TRUE)
  }
  
  # Remove genes without EntrezID
  tb <- de_results %>% 
    tidyr::drop_na(entrez_id) %>% 
    dplyr::distinct(entrez_id, .keep_all = TRUE)
  
  # Annotate according to pvalue
  if(mode == "normal"){
    all_genes <- tb %>% 
      dplyr::filter(pvalue < 0.05) %>% 
      dplyr::pull(entrez_id) 
    up_genes <- tb %>% 
      dplyr::filter(pvalue < 0.05 & log2fc > 0) %>% 
      dplyr::pull(entrez_id)
    down_genes <- tb %>% 
      dplyr::filter(pvalue < 0.05 & log2fc < 0) %>% 
      dplyr::pull(entrez_id)
  } else {
    all_genes <- tb %>% 
      dplyr::filter(pvalue_BH < 0.05) %>% 
      dplyr::pull(entrez_id) 
    up_genes <- tb %>% 
      dplyr::filter(pvalue_BH < 0.05 & log2fc > 0) %>% 
      dplyr::pull(entrez_id)
    down_genes <- tb %>% 
      dplyr::filter(pvalue_BH < 0.05 & log2fc < 0) %>% 
      dplyr::pull(entrez_id)
  }
  
  # MSigDB Hallmarks collection
  eH <- imap(list(All = all_genes, Up = up_genes, Down = down_genes), \(gene_list, direction) {
    print(glue::glue("MSigDB - Hallmarks - {direction} enrichment analysis for comparison {contrast}"))
    safe_enricher(gene_list, direction, H_genes)
    })
  eH_tb <- eH %>% 
    imap(\(x, idx) enrich2tb(x, idx)) %>% 
    list_rbind() %>% 
    dplyr::arrange(pvalue)
  
  # MSigDB Curated collection
  eC2 <- imap(list(All = all_genes, Up = up_genes, Down = down_genes), \(gene_list, direction) {
    print(glue::glue("MSigDB - Curated - {direction} enrichment analysis for comparison {contrast}"))
    safe_enricher(gene_list, direction, C2_genes)
  })
  eC2_tb <- eC2 %>% 
    imap(\(x, idx) enrich2tb(x, idx)) %>% 
    list_rbind() %>% 
    dplyr::arrange(pvalue)
  
  # GO
  eGO <- imap(list(All = all_genes, Up = up_genes, Down = down_genes), \(gene_list, direction) {
    map(list(MF = "MF", BP = "BP", CC = "CC"), \(ont) {
      print(glue::glue("GO - {ont} - {direction} enrichment analysis for comparison {contrast}"))
      safe_enrichGO(gene_list, direction, ont)
      })
    }) 
  
  eGO_tb <- eGO %>%
    imap(\(x, idx) {
      x %>%
        map(\(y) {possibly(\(y, idx) enrich2tb(y, idx), otherwise = NULL)(y, idx)}) %>% 
        list_rbind()
    }) %>%
    list_rbind() %>%
    arrange(pvalue)
  
  # MSigDB Oncogenic collection
  eC6 <- imap(list(All = all_genes, Up = up_genes, Down = down_genes), \(gene_list, direction) {
    print(glue::glue("MSigDB - Oncogenic - {direction} enrichment analysis for comparison {contrast}"))
    safe_enricher(gene_list, direction, C6_genes)
  }) 
  
  eC6_tb <- eC6 %>% 
    imap(\(x, idx) enrich2tb(x, idx)) %>% 
    list_rbind() %>% 
    dplyr::arrange(pvalue)
  
  # Summary
  summ <- bind_rows(eH_tb %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                    eC2_tb %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                    eGO_tb %>% 
                      dplyr::mutate(database = "GO"),
                    eC6_tb %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & as.numeric(str_split_i(BgRatio, "/", 1)) < 400) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(database, ID, GO, Description, GeneDirection, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, pvalue_BH, Count, geneID)
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "Hallmarks - MSigDB" = eH_tb,
                   "Curated - MSigDB" = eC2_tb,
                   "GO" = eGO_tb,
                   "Oncogenic - MSigDB" = eC6_tb)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/enrichment_analysis/{prefix}_{contrast}_enrichment_analysis.xlsx"))
  
  return(list(H = eH,
              C2 = eC2,
              GO = eGO,
              C6 = eC6))
}
res_enrich <- df %>% 
  imap(\(x, idx) enrichment_function(x, idx, mode="adjusted"))
qs_save(res_enrich, file = glue::glue("pathway_enrichment_analysis/{prefix}_enrichResult.qs2"))

# Camera ----
safe_camera <- possibly(\(camera_index) {
  camera(y = cpm_camera,
         index = camera_index,
         design = design_matrix$design,
         contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::select(category, direction, pvalue, pvalue_BH, NGenes)
}, otherwise = NULL)
camera_function <- function(contrast, counts){
  if(!(dir.exists(glue::glue("{out.dir}/camera_analysis")))){
    dir.create(glue::glue("{out.dir}/camera_analysis"))
  }
  
  # Prepare Camera
  cpm_camera <- edgeR::cpm(counts,
                    normalized.lib.sizes = TRUE,
                    log = TRUE,
                    prior.count = 3)
  rownames(cpm_camera) <- counts$genes$entrezid
  
  # MSigDB Hallmarks collection
  print(glue::glue("MSigDB - Hallmarks - Camera analysis for comparison {contrast}"))
  cameraH <- safe_camera(idxH)
  
  # MSigDB Curated collection
  print(glue::glue("MSigDB - Curated - Camera analysis for comparison {contrast}"))
  cameraC2 <- safe_camera(idxC2)
  
  # GO
  print(glue::glue("MSigDB - GO - Camera analysis for comparison {contrast}"))
  cameraGO <- safe_camera(idxC5)
  
  # MSigDB Oncogenic collection
  print(glue::glue("MSigDB - Oncogenic - Camera analysis for comparison {contrast}"))
  cameraC6 <- safe_camera(idxC6)
  
  # Summary
  summ <- bind_rows(cameraH %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                      cameraC2 %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                    cameraGO %>% 
                      dplyr::mutate(database = "GO"),
                      cameraC6 %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & NGenes < 400) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(database, category, direction, pvalue, pvalue_BH, global_pvalue_BH, NGenes)
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "Hallmarks - MSigDB" = cameraH,
                   "Curated - MSigDB" = cameraC2,
                   "GO - MSigDB" = cameraGO,
                   "Oncogenic - MSigDB" = cameraC6)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/camera_analysis/{prefix}_{contrast}_camera_analysis.xlsx"))
  
  return(summ)
}
names(df) %>% 
  walk(\(x) camera_function(x, counts))

# GSEA ----
gsea_function <- function(de_results, contrast){
  if(!(dir.exists(glue::glue("{out.dir}/gsea_analysis")))){
    dir.create(glue::glue("{out.dir}/gsea_analysis"))
  }
  
  # Rank the genes for the GSEA
  ranks <- de_results %>% 
    tidyr::drop_na(entrez_id) %>% 
    dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
    dplyr::mutate(entrez_id = entrez_id,
                  ranks = ifelse("t" %in% colnames(de_results),
                                 t, #From fgsea, ranking using limma t-statistic:10.1101/060012
                                 log2fc * (-log10(pvalue))), #Rank significant genes: 10.1093/bioinformatics/btr671
                  .keep = "none") %>%  
    dplyr::arrange(ranks) %>%
    tibble::deframe()
  
  # KEGG
  print(glue::glue("KEGG GSEA analysis for comparison {contrast}"))
  gseaKEGG <- fgseaMultilevel(pathways = idxKEGG,
                              stats = ranks,
                              maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(category = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::left_join(kegg2paths, by = "category") %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(category, path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # Reactome
  print(glue::glue("Reactome GSEA analysis for comparison {contrast}"))
  gseaReactome <- fgseaMultilevel(pathways = idxReactome,
                                  stats = ranks,
                                  maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(category = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::left_join(reactome2paths, by = c("category" = "DB_ID")) %>% 
    dplyr::rename(path = path_name) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(category, path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # GO
  print(glue::glue("GO GSEA analysis for comparison {contrast}"))
  gseaGO <- fgseaMultilevel(pathways = C5_genes,
                            stats = ranks,
                            maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(path = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # MSigDB Hallmarks collection
  print(glue::glue("MSigDB GSEA analysis for comparison {contrast}"))
  gseaH <- fgseaMultilevel(pathways = H_genes,
                           stats = ranks,
                           maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(path = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # MSigDB Curated collection
  gseaC2 <- fgseaMultilevel(pathways = C2_genes,
                            stats = ranks,
                            maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(path = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # MSigDB Oncogenic collection
  gseaC6 <- fgseaMultilevel(pathways = C6_genes,
                            stats = ranks,
                            maxSize = 500) %>% 
    as_tibble() %>% 
    dplyr::rename(path = pathway,
                  pvalue = pval,
                  pvalue_BH = padj) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(path, ES, NES, pvalue, pvalue_BH, pvalue_bonferroni, log2err, size)
  
  # Summary
  summ <- bind_rows(gseaKEGG %>% 
                        dplyr::mutate(database = "KEGG"),
                      gseaReactome %>% 
                        dplyr::mutate(database = "Reactome"),
                      gseaGO %>% 
                        dplyr::mutate(database = "GO"),
                      gseaH %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                      gseaC2 %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                      gseaC6 %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH"),
                  global_pvalue_bonferroni = p.adjust(pvalue,
                                                      method = "bonferroni")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & size < 400) %>% 
    dplyr::arrange(global_pvalue_BH) %>% 
    dplyr::select(database, category, path, ES, NES, pvalue, pvalue_BH, global_pvalue_BH, pvalue_bonferroni, global_pvalue_bonferroni, log2err, size)
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "KEGG" = gseaKEGG,
                   "Reactome" = gseaReactome,
                   "GO - MSigDB" = gseaGO,
                   "Hallmarks - MSigDB" = gseaH,
                   "Curated - MSigDB" = gseaC2,
                   "Oncogenic - MSigDB" = gseaC6)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/gsea_analysis/{prefix}_{contrast}_gsea_analysis.xlsx"))
}
df %>% 
  iwalk(\(x, idx) gsea_function(x, idx))

## GSEA Plots ----
# plotEnrichment(H_genes[["HALLMARK_MYC_TARGETS_V1"]], ranks) +
#   labs(title = "HALLMARK_MYC_TARGETS_V1")
