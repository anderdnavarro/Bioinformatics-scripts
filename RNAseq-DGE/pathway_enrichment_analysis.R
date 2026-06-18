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
safe_gseGO <- possibly(\(gene_rank, ont, simplify.f=FALSE) {
  tmp <- gseGO(geneList = gene_rank,
        OrgDb = org.Hs.eg.db,
        ont = ont,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        method = "multilevel")
  tmp <- enrichplot::pairwise_termsim(tmp)
  
  if (simplify.f) {
    suppressWarnings(clusterProfiler::simplify(tmp, cutoff=0.7, by="p.adjust", select_fun=min)) %>%
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  } else {
    tmp %>% 
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
}, otherwise = NULL)
safe_gsea <- possibly(\(gene_rank, term) {
  GSEA(geneList = gene_rank,
       TERM2GENE = term,
       minGSSize = 10,
       maxGSSize = 500,
       pvalueCutoff = 1,
       pAdjustMethod = "BH",
       method = "multilevel") %>% 
    setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
}, otherwise = NULL)
gsea2tb <- function(gseaResult) {
  ont <- ifelse(gseaResult@setType == "unknown", NA, gseaResult@setType)
  gseaResult %>% 
    as_tibble() %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::mutate(GO = ont,
                  enrichSize = str_count(core_enrichment, "/") + 1) %>% 
    dplyr::rename(pvalue_BH = p.adjust) %>%
    dplyr::select(ID, GO, Description, enrichmentScore, NES, pvalue, pvalue_BH, log2err, setSize, enrichSize, leading_edge, core_enrichment)
}
gsea_function <- function(de_results, contrast){
  if(!(dir.exists(glue::glue("{out.dir}/gsea_analysis")))){
    dir.create(glue::glue("{out.dir}/gsea_analysis"))
  }
  
  # Rank the genes for the GSEA
  ranks <- de_results %>% 
    tidyr::drop_na(entrez_id) %>% 
    dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
    {
      if ("t" %in% colnames(.)) {
        dplyr::mutate(., entrez_id,
                      ranks = t, #From fgsea, ranking using limma t-statistic:10.1101/060012
                      .keep = "none")
      } else {
        dplyr::mutate(., entrez_id,
                      ranks = log2fc * (-log10(pvalue)), #Rank significant genes: 10.1093/bioinformatics/btr671
                      .keep = "none")
      }
    } %>% 
    dplyr::arrange(desc(ranks)) %>%
    tibble::deframe()
  
  # MSigDB Hallmarks collection
  print(glue::glue("MSigDB - Hallmarks - GSEA analysis for comparison {contrast}"))
  gseaH <- safe_gsea(ranks, H_genes)
  gseaH_tb <- gsea2tb(gseaH)
  
  # MSigDB Curated collection
  print(glue::glue("MSigDB - Curated - GSEA analysis for comparison {contrast}"))
  gseaC2 <- safe_gsea(ranks, C2_genes)
  gseaC2_tb <- gsea2tb(gseaC2)
  
  # GO
  gseaGO <- map(list(MF = "MF", BP = "BP", CC = "CC"), \(ont) {
    print(glue::glue("GO - {ont} - GSEA analysis for comparison {contrast}"))
    safe_gseGO(ranks, ont, simplify.f = FALSE) #Simplify GO terms takes a lot of time for GSEA (I think it's because we don't filter by significant terms before simplifying)
  }) 
  gseaGO_tb <- gseaGO %>% 
    map(\(x) gsea2tb(x)) %>% 
    list_rbind() %>% 
    dplyr::arrange(pvalue)
  
  # MSigDB Oncogenic collection
  print(glue::glue("MSigDB - Oncogenic - GSEA analysis for comparison {contrast}"))
  gseaC6 <- safe_gsea(ranks, C6_genes)
  gseaC6_tb <- gsea2tb(gseaC6)
  
  # Summary
  summ <- bind_rows(gseaH_tb %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                    gseaC2_tb %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                    gseaGO_tb %>% 
                      dplyr::mutate(database = "GO"),
                    gseaC6_tb %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & setSize < 400) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(database, ID, GO, Description, enrichmentScore, NES, pvalue, pvalue_BH, global_pvalue_BH, log2err, setSize, enrichSize, leading_edge, core_enrichment) 
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "Hallmarks - MSigDB" = gseaH_tb,
                   "Curated - MSigDB" = gseaC2_tb,
                   "GO" = gseaGO_tb,
                   "Oncogenic - MSigDB" = gseaC6_tb)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/gsea_analysis/{prefix}_{contrast}_gsea_analysis.xlsx"))
  
  return(list(H = gseaH,
              C2 = gseaC2,
              GO = gseaGO,
              C6 = gseaC6))
}
res_gsea <- df %>% 
  imap(\(x, idx) gsea_function(x, idx))
qs_save(res_gsea, file = glue::glue("pathway_enrichment_analysis/{prefix}_enrichGSEA.qs2"))

## GSEA Plots ----
# plotEnrichment(H_genes[["HALLMARK_MYC_TARGETS_V1"]], ranks) +
#   labs(title = "HALLMARK_MYC_TARGETS_V1")
