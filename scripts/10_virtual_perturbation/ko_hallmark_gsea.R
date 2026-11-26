rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

cran_pkgs <- c("data.table", "dplyr", "stringr", "purrr", "tibble")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("fgsea", "msigdbr")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(fgsea)
library(msigdbr)

base_dir <- Sys.getenv("OAGLY_KO_ROOT", unset = project_path("results_root", "virtual_ko"))
out_dir  <- file.path(base_dir, "HALLMARK_GSEA_EnrichmentMap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

minSize <- 10
maxSize <- 500
nperm   <- 10000
padj_cutoff <- 0.05

min_KO_freq <- 4

top_n_consensus <- 100

ko_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

ko_dirs <- ko_dirs[str_detect(basename(ko_dirs), "敲除")]

if (length(ko_dirs) == 0) {
  stop("未找到 KO 文件夹，请检查 base_dir。")
}

message("找到 KO 文件夹数: ", length(ko_dirs))

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")

hallmark_list <- msig_h %>%
  distinct(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop")

hallmark_pathways <- setNames(hallmark_list$genes, hallmark_list$gs_name)

gmt_path <- file.path(out_dir, "HALLMARK_msigdb_export.gmt")

write_gmt <- function(pathways, file) {
  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  for (nm in names(pathways)) {
    line <- paste(c(nm, "na", pathways[[nm]]), collapse = "\t")
    writeLines(line, con)
  }
}
write_gmt(hallmark_pathways, gmt_path)

extract_ko_gene <- function(folder_name) {

  x <- folder_name %>%
    str_remove("^\\d+_") %>%
    str_remove("敲除$")
  x
}

read_diffreg <- function(f, ko_gene) {
  dt <- fread(f)

  need_cols <- c("gene", "Z")
  miss_cols <- setdiff(need_cols, colnames(dt))
  if (length(miss_cols) > 0) {
    stop("文件缺少必要列: ", paste(miss_cols, collapse = ", "), " | 文件: ", f)
  }

  dt2 <- dt %>%
    transmute(
      gene = as.character(gene),
      Z = as.numeric(Z),
      p.adj = if ("p.adj" %in% colnames(dt)) as.numeric(.data[["p.adj"]]) else NA_real_,
      Significance = if ("Significance" %in% colnames(dt)) as.character(.data[["Significance"]]) else NA_character_
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(Z)) %>%
    group_by(gene) %>%
    summarise(Z = mean(Z, na.rm = TRUE), .groups = "drop")

  stats <- dt2$Z
  names(stats) <- dt2$gene

  stats <- sort(stats, decreasing = TRUE)

  return(list(table = dt2, stats = stats))
}

run_hallmark_fgsea <- function(stats_vec, ko_gene) {

  fg <- suppressWarnings(
    fgseaMultilevel(
      pathways = hallmark_pathways,
      stats = stats_vec,
      minSize = minSize,
      maxSize = maxSize,
      eps = 0
    )
  )

  fg <- as.data.table(fg)

  if (nrow(fg) == 0) {
    return(data.frame())
  }

  fg2 <- fg %>%
    as_tibble() %>%
    mutate(
      KO_gene = ko_gene,
      Description = pathway,
      GeneSetID = pathway,
      pvalue = pval,
      FDR = padj,
      Phenotype = ifelse(NES >= 0, 1, -1),
      leadingEdge_genes = map_chr(leadingEdge, ~ paste(.x, collapse = ",")),
      core_size = map_int(leadingEdge, length)
    ) %>%
    select(
      KO_gene, GeneSetID, Description,
      pvalue, FDR, ES, NES, size,
      leadingEdge_genes, core_size, Phenotype
    ) %>%
    arrange(FDR, desc(abs(NES)))

  return(fg2)
}

to_emap_generic <- function(df) {
  df %>%
    transmute(
      GeneSetID   = GeneSetID,
      Description = Description,
      pvalue      = pvalue,
      FDR         = FDR,
      Phenotype   = Phenotype,
      Genes       = leadingEdge_genes
    )
}

all_gsea_list <- list()

per_ko_dir <- file.path(out_dir, "per_KO_results")
dir.create(per_ko_dir, showWarnings = FALSE, recursive = TRUE)

for (d in ko_dirs) {
  folder_name <- basename(d)
  ko_gene <- extract_ko_gene(folder_name)

  diff_file <- file.path(d, "diffRegulation_full.txt")
  if (!file.exists(diff_file)) {
    warning("缺少 diffRegulation_full.txt，跳过: ", d)
    next
  }

  message("正在处理: ", folder_name, " | KO gene = ", ko_gene)

  x <- read_diffreg(diff_file, ko_gene)
  gsea_res <- run_hallmark_fgsea(x$stats, ko_gene)

  if (nrow(gsea_res) == 0) {
    warning("GSEA 无结果: ", folder_name)
    next
  }

  gsea_sig <- gsea_res %>%
    filter(!is.na(FDR), FDR <= padj_cutoff)

  ko_out_subdir <- file.path(per_ko_dir, folder_name)
  dir.create(ko_out_subdir, showWarnings = FALSE, recursive = TRUE)

  fwrite(gsea_res,
         file.path(ko_out_subdir, paste0(folder_name, "_HALLMARK_GSEA_all.csv")))

  fwrite(gsea_sig,
         file.path(ko_out_subdir, paste0(folder_name, "_HALLMARK_GSEA_sig.csv")))

  if (nrow(gsea_sig) > 0) {
    emap_ko <- to_emap_generic(gsea_sig)
    fwrite(emap_ko,
           file.path(ko_out_subdir, paste0(folder_name, "_EnrichmentMap_generic.txt")),
           sep = "\t", quote = FALSE)
  }

  all_gsea_list[[folder_name]] <- gsea_res
}

all_gsea <- bind_rows(all_gsea_list)

if (nrow(all_gsea) == 0) {
  stop("所有 KO 的 GSEA 结果均为空，请检查输入。")
}

all_gsea_sig <- all_gsea %>%
  filter(!is.na(FDR), FDR <= padj_cutoff)

fwrite(all_gsea, file.path(out_dir, "All_KO_HALLMARK_GSEA_all.csv"))
fwrite(all_gsea_sig, file.path(out_dir, "All_KO_HALLMARK_GSEA_sig.csv"))

consensus_summary <- all_gsea_sig %>%
  group_by(GeneSetID, Description) %>%
  summarise(
    n_KO = n_distinct(KO_gene),
    KO_list = paste(sort(unique(KO_gene)), collapse = ","),
    mean_FDR = mean(FDR, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    mean_absNES = mean(abs(NES), na.rm = TRUE),
    pos_count = sum(NES > 0, na.rm = TRUE),
    neg_count = sum(NES < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KO), mean_FDR, desc(mean_absNES))

fwrite(consensus_summary, file.path(out_dir, "Consensus_HALLMARK_summary.csv"))

consensus_keep <- consensus_summary %>%
  filter(n_KO >= min_KO_freq) %>%
  slice_head(n = top_n_consensus)

message("进入共识 EnrichmentMap 的通路数: ", nrow(consensus_keep))

consensus_emap <- all_gsea_sig %>%
  semi_join(consensus_keep, by = c("GeneSetID", "Description")) %>%
  group_by(GeneSetID, Description) %>%
  summarise(
    pvalue = min(pvalue, na.rm = TRUE),
    FDR = min(FDR, na.rm = TRUE),

    Phenotype = ifelse(sum(NES > 0, na.rm = TRUE) >= sum(NES < 0, na.rm = TRUE), 1, -1),
    Genes = paste(sort(unique(unlist(strsplit(paste(leadingEdge_genes, collapse = ","), ",")))), collapse = ","),
    n_KO = n_distinct(KO_gene),
    KO_list = paste(sort(unique(KO_gene)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KO), FDR)

consensus_emap_export <- consensus_emap %>%
  select(GeneSetID, Description, pvalue, FDR, Phenotype, Genes)

fwrite(consensus_emap_export,
       file.path(out_dir, "EnrichmentMap_HALLMARK_consensus_generic.txt"),
       sep = "\t", quote = FALSE)

fwrite(consensus_emap,
       file.path(out_dir, "EnrichmentMap_HALLMARK_consensus_with_meta.csv"))

ko_sig_count <- all_gsea_sig %>%
  group_by(KO_gene) %>%
  summarise(
    n_sig_pathways = n(),
    top_pathways = paste(head(Description[order(FDR, -abs(NES))], 5), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_pathways))

fwrite(ko_sig_count, file.path(out_dir, "KO_sig_pathway_counts.csv"))

ko_pathway_matrix <- all_gsea %>%
  mutate(score = ifelse(is.na(FDR), 0, ifelse(FDR <= padj_cutoff, sign(NES) * -log10(FDR + 1e-300), 0))) %>%
  select(KO_gene, Description, score) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = Description, values_from = score, values_fill = 0)

fwrite(as.data.table(ko_pathway_matrix), file.path(out_dir, "KO_x_HALLMARK_score_matrix.csv"))

cat(
  "================ Cytoscape / EnrichmentMap 使用说明 ================\n",
  "1. 打开 Cytoscape\n",
  "2. Apps -> App Manager，安装 EnrichmentMap\n",
  "3. Apps -> EnrichmentMap\n",
  "4. 选择 Generic enrichment result\n",
  "5. 导入以下文件：\n",
  "   - EnrichmentMap_HALLMARK_consensus_generic.txt\n",
  "   - GMT 文件: ", gmt_path, "\n",
  "6. 初始阈值建议：\n",
  "   - FDR cutoff = 0.05\n",
  "   - p-value cutoff = 0.05\n",
  "   - overlap cutoff 先用中等\n",
  "7. 观察模块是否聚成：\n",
  "   - Glycolysis / metabolism\n",
  "   - ECM / matrix remodeling\n",
  "   - Hypertrophy / differentiation\n",
  "   - Stress / inflammatory response\n",
  "===================================================================\n",
  sep = ""
)

message("全部完成。结果输出到: ", out_dir)
