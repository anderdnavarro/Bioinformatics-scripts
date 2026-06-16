library(qs2)
library(openxlsx)
library(dplyr)
library(tidyr)
library(purrr)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(reactome.db)
library(ReactomePA)
library(goseq)
library(fgsea)
library(ggplot2)
library(ggrepel)

# Constants ----
lfcCutoff <- 0.5
pCutoff <- 0.05
prefix <- '20260616_default_protocol'
out.dir <- 'pathway_enrichment_analysis'
de.file <- 'default_results/20260209_default_protocol_all_contrasts_differential_expression_results_edgeR4.qs2'
counts.file <- 'default_results/dge_matrix.qs2'
design.file <- 'default_results/design_matrix.qs2'

# Open results ----
df <- qs_read(de.file)
counts <- qs_read(counts.file)
design_matrix <- qs_read(design.file)

# KEGG and Reactome datasets ----
kegg2paths <- read_delim('/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/KEGG/kegg_pathways.txt',
                         show_col_types = FALSE,
                         col_names = c("category", "path")) %>% 
  dplyr::mutate(path = str_remove(path, " - Homo sapiens \\(human\\)"))
kegg2genes <- read_delim('/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/KEGG/kegg_genes.txt',
                         show_col_types = FALSE,
                         col_names = c("category", "type", "genomic_coordinates", "genes")) %>% 
  dplyr::select(category, genes) %>% 
  dplyr::mutate(genes = str_split_i(genes, '[,;]', 1))
kegg2entrez <- as.list(org.Hs.egPATH2EG)
names(kegg2entrez) <- glue::glue("hsa{names(kegg2entrez)}")
idxKEGG <- ids2indices(gene.sets = kegg2entrez,
                       identifiers = counts$genes$entrezid)
reactome2paths <- as.data.frame(reactomePATHID2NAME)
reactome2entrez <- as.list(reactomePATHID2EXTID)
idxReactome <- ids2indices(gene.sets = reactome2entrez,
                           identifiers = counts$genes$entrezid)

# MSigDB datasets ----
H_genes <- readRDS("/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/MSigDB/Hs.h.all.v7.1.entrez.rds") #Hallmarks collection
idxH <- ids2indices(gene.sets = H_genes,
                    identifiers = counts$genes$entrezid)
C2_genes <- readRDS("/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/MSigDB/Hs.c2.all.v7.1.entrez.rds") #Curated collection
idxC2 <- ids2indices(gene.sets = C2_genes,
                     identifiers = counts$genes$entrezid)
C5_genes <- readRDS("/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/MSigDB/Hs.c5.all.v7.1.entrez.rds") #GO collection
idxC5 <- ids2indices(gene.sets = C5_genes,
                     identifiers = counts$genes$entrezid)
C6_genes <- readRDS("/Users/ander/Library/Mobile Documents/com~apple~CloudDocs/Protocolos/RNAseq/databases/MSigDB/Hs.c6.all.v7.1.entrez.rds") #Oncogenic collection
idxC6 <- ids2indices(gene.sets = C6_genes,
                     identifiers = counts$genes$entrezid)

# Enrichment analysis ----
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
    eGenes <- tb %>% 
      dplyr::mutate(entrez_id = entrez_id,
                    eGene = as.integer(pvalue < 0.05),
                    .keep = "none") %>% 
      tibble::deframe()
  } else {
    eGenes <- tb %>% 
      dplyr::mutate(entrez_id = entrez_id,
                    eGene = as.integer(pvalue_BH < 0.05),
                    .keep = "none") %>% 
      tibble::deframe()
  }
  
  # Genes length
  length_bias <- tb %>% 
    dplyr::select(entrez_id, length) %>% 
    tibble::deframe()
  
  # Fitting the Probability Weighting Function (PWF - goseq) 
  pwf <- nullp(DEgenes = eGenes,
               bias.data = length_bias,
               plot.fit  = FALSE)
  
  # KEGG
  print(glue::glue("KEGG enrichment analysis for comparison {contrast}"))
  eKEGG <- goseq(pwf,
                 genome = "hg38",
                 id = "knownGene",
                 test.cats = c("KEGG")) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni"),
                  category = glue::glue("hsa{category}")) %>% 
    dplyr::left_join(kegg2paths, by = "category") %>% 
    dplyr::select(category, path, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # Reactome
  print(glue::glue("Reactome enrichment analysis for comparison {contrast}"))
  eReactome <- goseq(pwf,
                     gene2cat = as.list(reactomeEXTID2PATHID)) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::left_join(reactome2paths, by = c("category" = "DB_ID")) %>% 
    dplyr::rename(path = path_name) %>% 
    dplyr::select(category, path, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # GO
  print(glue::glue("GO enrichment analysis for comparison {contrast}"))
  eGO <- goseq(pwf,
               genome = "hg38",
               id = "knownGene",
               test.cats = c("GO:MF", "GO:BP", "GO:CC")) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue,
                  path = term) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::select(category, path, ontology, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # MSigDB Hallmarks collection
  print(glue::glue("MSigDB enrichment analysis for comparison {contrast}"))
  eH_genes <- goseq(pwf,
                    gene2cat = H_genes) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::select(category, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # MSigDB Curated collection
  eC2_genes <- goseq(pwf,
                     gene2cat = C2_genes) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::select(category, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # MSigDB Oncogenic collection
  eC6_genes <- goseq(pwf,
                     gene2cat = C6_genes) %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = over_represented_pvalue) %>% 
    dplyr::mutate(pvalue_BH = p.adjust(pvalue, method = "BH"),
                  pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::select(category, pvalue, pvalue_BH, pvalue_bonferroni, numDEInCat, numInCat)
  
  # Summary
  summ <- bind_rows(eKEGG %>% 
                        dplyr::mutate(database = "KEGG"),
                      eReactome %>% 
                        dplyr::mutate(database = "Reactome"),
                      eGO %>% 
                        dplyr::mutate(database = "GO"),
                      eH_genes %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                      eC2_genes %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                      eC6_genes %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH"),
                  global_pvalue_bonferroni = p.adjust(pvalue,
                                                      method = "bonferroni")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & numInCat < 400) %>% 
    dplyr::arrange(global_pvalue_BH) %>% 
    dplyr::select(database, category, path, ontology, pvalue, pvalue_BH, global_pvalue_BH, pvalue_bonferroni, global_pvalue_bonferroni, numDEInCat, numInCat)
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "KEGG" = eKEGG,
                   "Reactome" = eReactome,
                   "GO" = eGO,
                   "Hallmarks - MSigDB" = eH_genes,
                   "Curated - MSigDB" = eC2_genes,
                   "Oncogenic - MSigDB" = eC6_genes)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/enrichment_analysis/{prefix}_{contrast}_enrichment_analysis.xlsx"))
}
df %>% 
  iwalk(\(x, idx) enrichment_function(x, idx, mode="adjusted"))

# Camera ----
camera_function <- function(contrast, counts){
  if(!(dir.exists(glue::glue("{out.dir}/camera_analysis")))){
    dir.create(glue::glue("{out.dir}/camera_analysis"))
  }
  
  # Prepare Camera
  cpm_camera <- cpm(counts,
                    normalized.lib.sizes = TRUE,
                    log = TRUE,
                    prior.count = 3)
  rownames(cpm_camera) <- counts$genes$entrezid
  
  # KEGG
  print(glue::glue("KEGG camera analysis for comparison {contrast}"))
  cameraKEGG <- camera(y = cpm_camera,
                       index = idxKEGG,
                       design = design_matrix$design,
                       contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::left_join(kegg2paths, by = "category") %>% 
    dplyr::select(category, path, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # Reactome
  print(glue::glue("Reactome camera analysis for comparison {contrast}"))
  cameraReactome <- camera(y = cpm_camera,
                           index = idxReactome,
                           design = design_matrix$design,
                           contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::left_join(reactome2paths, by = c("category" = "DB_ID")) %>% 
    dplyr::rename(path = path_name) %>% 
    dplyr::select(category, path, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # GO
  print(glue::glue("GO camera analysis for comparison {contrast}"))
  cameraGO <- camera(y = cpm_camera,
                     index = idxC5,
                     design = design_matrix$design,
                     contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::select(category, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # MSigDB Hallmarks collection
  print(glue::glue("MSigDB camera analysis for comparison {contrast}"))
  cameraH <- camera(y = cpm_camera,
                    index = idxH,
                    design = design_matrix$design,
                    contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::select(category, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # MSigDB Curated collection
  cameraC2 <- camera(y = cpm_camera,
                     index = idxC2,
                     design = design_matrix$design,
                     contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::select(category, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # MSigDB Oncogenic collection
  cameraC6 <- camera(y = cpm_camera,
                     index = idxC6,
                     design = design_matrix$design,
                     contrast = design_matrix$contr.matrix[,contrast]) %>% 
    tibble::rownames_to_column(var = "category") %>% 
    as_tibble() %>% 
    dplyr::rename(pvalue = PValue,
                  pvalue_BH = FDR,
                  direction = Direction) %>% 
    dplyr::mutate(pvalue_bonferroni = p.adjust(pvalue, method = "bonferroni")) %>% 
    dplyr::filter(pvalue < 0.05) %>% 
    dplyr::select(category, direction, pvalue, pvalue_BH, pvalue_bonferroni, NGenes)
  
  # Summary
  summ <- bind_rows(cameraKEGG %>% 
                        dplyr::mutate(database = "KEGG"),
                      cameraReactome %>% 
                        dplyr::mutate(database = "Reactome"),
                      cameraGO %>% 
                        dplyr::mutate(database = "GO"),
                      cameraH %>% 
                        dplyr::mutate(database = "Hallmarks - MSigDB"),
                      cameraC2 %>% 
                        dplyr::mutate(database = "Curated - MSigDB"),
                      cameraC6 %>% 
                        dplyr::mutate(database = "Oncogenic - MSigDB")) %>% 
    dplyr::mutate(global_pvalue_BH = p.adjust(pvalue,
                                              method = "BH"),
                  global_pvalue_bonferroni = p.adjust(pvalue,
                                                      method = "bonferroni")) %>% 
    dplyr::filter(global_pvalue_BH < 0.05 & NGenes < 400) %>% 
    dplyr::arrange(global_pvalue_BH) %>% 
    dplyr::select(database, category, path, direction, pvalue, pvalue_BH, global_pvalue_BH, pvalue_bonferroni, global_pvalue_bonferroni, NGenes)
  
  # Write the results in an excel file
  dfs_list <- list("Summary" = summ,
                   "KEGG" = cameraKEGG,
                   "Reactome" = cameraReactome,
                   "GO - MSigDB" = cameraGO,
                   "Hallmarks - MSigDB" = cameraH,
                   "Curated - MSigDB" = cameraC2,
                   "Oncogenic - MSigDB" = cameraC6)
  write.xlsx(dfs_list, file = glue::glue("{out.dir}/camera_analysis/{prefix}_{contrast}_camera_analysis.xlsx"))
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
