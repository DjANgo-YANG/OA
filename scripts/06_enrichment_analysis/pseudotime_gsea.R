rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

INPUT_RDS <- project_path("data_root", "auxiliary/oa/slingshot/rds/OA_Seurat_with_slingshot.rds")
OUT_DIR   <- project_path("results_root", "pseudotime_enrichment")
TIME_COL  <- "pseudotime_L5"
SPECIES   <- "Homo sapiens"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

packages <- c("Seurat", "clusterProfiler", "enrichplot", "msigdbr", "ggplot2", "dplyr", "tibble")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

message("正在准备 MSigDB 基因集...")

m_H <- msigdbr(species = SPECIES, category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

m_C5 <- msigdbr(species = SPECIES, category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

m_t2g <- rbind(m_H, m_C5)

message("正在读取 Seurat 对象...")
seu <- readRDS(INPUT_RDS)

if (!TIME_COL %in% colnames(seu@meta.data)) stop(paste("找不到列:", TIME_COL))

message("正在计算基因与拟时序的相关性 (这可能需要几分钟)...")

expr_mat <- GetAssayData(seu, slot = "data", assay = "RNA")
keep_genes <- rowSums(expr_mat > 0) > 10
expr_mat <- expr_mat[keep_genes, ]

pt_info <- seu@meta.data[[TIME_COL]]
names(pt_info) <- rownames(seu@meta.data)
valid_cells <- names(pt_info)[!is.na(pt_info)]

expr_use <- as.matrix(expr_mat[, valid_cells])
pt_use   <- pt_info[valid_cells]

gene_cor <- cor(t(expr_use), pt_use, method = "spearman")

gene_list <- gene_cor[, 1]
gene_list <- sort(gene_list, decreasing = TRUE)

message("排序列表构建完成。Top 5 正相关基因 (MMP13+端):")
print(head(gene_list, 5))

write.csv(data.frame(Gene=names(gene_list), Correlation=gene_list),
          file.path(OUT_DIR, "Gene_Pseudotime_Correlation_Rank.csv"), row.names = FALSE)

message("正在运行 GSEA (可能需要几分钟)...")

set.seed(1234)

gsea_res <- GSEA(gene_list,
                 TERM2GENE = m_t2g,
                 pvalueCutoff = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 0,
                 seed = TRUE)

write.csv(gsea_res@result, file.path(OUT_DIR, "GSEA_Full_Results.csv"))

gsea_sig <- filter(gsea_res, p.adjust < 0.05)
message(paste("共发现", nrow(gsea_sig), "条显著富集通路"))

write.csv(gsea_sig, file.path(OUT_DIR, "GSEA_Significant_Results.csv"))

keywords <- c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TNFA", "HYPOXIA", "MATRIX", "COLLAGEN", "ECM")

target_paths <- grep(paste(keywords, collapse="|"), gsea_sig$ID, value=TRUE, ignore.case=TRUE)

if (length(target_paths) > 0) {
  message("发现以下目标通路，正在绘制 GseaPlot:")
  print(head(target_paths))

  for (pathID in head(target_paths, 10)) {
    p <- gseaplot2(gsea_res, geneSetID = pathID, title = pathID, pvalue_table = TRUE)

    fname <- gsub("[:/]", "_", pathID)
    ggsave(file.path(OUT_DIR, paste0("GSEA_Plot_", fname, ".png")), p, width = 8, height = 6)
    ggsave(file.path(OUT_DIR, paste0("GSEA_Plot_", fname, ".pdf")), p, width = 8, height = 6)
  }
} else {
  warning("在显著结果中未找到包含 EMT/Matrix/Inflammation 的通路。")
}

top_up <- gsea_sig %>% arrange(desc(NES)) %>% head(10) %>% pull(ID)
top_down <- gsea_sig %>% arrange(NES) %>% head(10) %>% pull(ID)
top_paths <- c(top_up, top_down)

if (length(top_paths) > 0) {
  p_dot <- dotplot(gsea_res, showCategory = top_paths, split=".sign") +
    facet_grid(.~.sign) +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle(paste0("Top GSEA Pathways along ", TIME_COL))

  ggsave(file.path(OUT_DIR, "GSEA_Top_Dotplot.png"), p_dot, width = 12, height = 10)
  ggsave(file.path(OUT_DIR, "GSEA_Top_Dotplot.pdf"), p_dot, width = 12, height = 10)
}

if ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" %in% gsea_sig$ID) {
  emt_nes <- gsea_sig[gsea_sig$ID == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "NES"]
  if(emt_nes > 0) {
    message("\n>>> 【验证成功】 >>>")
    message("检测到 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION' 显著富集且 NES > 0。")
    message("这表明随着拟时序增加 (MMP13- -> MMP13+)，EMT 通路被激活。")
    message("这与 SCENIC 发现的 ZEB1 (EMT 关键 TF) 活性升高完美吻合！")
  }
}

message("=========================================")
message("GSEA 分析完成！结果已保存至: ", OUT_DIR)
