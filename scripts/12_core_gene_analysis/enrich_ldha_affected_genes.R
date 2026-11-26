rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
})

IN_FILE  <- project_path("data_root", "figure_inputs/ldha_ko/diffRegulation_full.txt")
OUT_DIR  <- file.path(dirname(IN_FILE), "Enrichment_LDHA_VirtualKO")
SPECIES  <- "human"
P_ADJ_CUT <- 0.05
TOPN_ORA  <- 20
TOPN_GSEA <- 20

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

pick_col <- function(df, candidates) {
  cn <- tolower(colnames(df))
  hit <- intersect(tolower(candidates), cn)
  if (length(hit) == 0) return(NA_character_)
  colnames(df)[match(hit[1], cn)]
}

guess_org <- function(species){
  if (tolower(species) %in% c("human","hs","hsa","homo sapiens")) {
    list(org_pkg="org.Hs.eg.db", kegg="hsa")
  } else if (tolower(species) %in% c("mouse","mm","mmu","mus musculus")) {
    list(org_pkg="org.Mm.eg.db", kegg="mmu")
  } else {
    stop("SPECIES must be 'human' or 'mouse'.")
  }
}

safe_pdf <- function(file, width=9, height=7) {
  grDevices::pdf(file, width=width, height=height)
}

save_gg <- function(p, filename, w=9, h=7) {
  ggsave(filename = filename, plot = p, width = w, height = h, dpi = 300)
}

cat("Reading:", IN_FILE, "\n")
df <- data.table::fread(IN_FILE, data.table = FALSE)
stopifnot(nrow(df) > 0)

gene_col <- pick_col(df, c("gene","genes","symbol","feature","features","genename","gene_name"))
if (is.na(gene_col)) stop("Cannot find gene column. Please ensure a column named gene/symbol/feature exists.")

rank_col <- pick_col(df, c("log2fc","logfc","avg_log2fc","avg_logfc","effect","score","z","zscore","stat","t","lr","w"))
if (is.na(rank_col)) {
  stop("Cannot find rank column (logFC/score/z/stat). Please ensure one exists for GSEA ranking.")
}

padj_col <- pick_col(df, c("p_val_adj","p_adj","padj","fdr","qval","q_value","p.adjust","p.adj","pval_adj"))

if (is.na(padj_col)) {
  warning("No adjusted p-value column found. ORA will use top genes by |rank| as a proxy.")
}

df2 <- df %>%
  dplyr::select(gene = all_of(gene_col),
                rank = all_of(rank_col),
                dplyr::everything()) %>%
  dplyr::mutate(gene = as.character(gene),
                rank = as.numeric(rank)) %>%
  dplyr::filter(!is.na(gene) & gene != "" & !is.na(rank)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

cat("Rows after cleaning:", nrow(df2), "\n")
write.csv(df2, file = "input_cleaned.csv", row.names = FALSE)

org <- guess_org(SPECIES)
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
})
if (org$org_pkg == "org.Mm.eg.db") suppressPackageStartupMessages(library(org.Mm.eg.db))

orgdb <- get(org$org_pkg)

map <- suppressWarnings(
  bitr(df2$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
)

if (is.null(map) || nrow(map) == 0) {
  stop("No genes could be mapped from SYMBOL to ENTREZID. Check SPECIES or gene symbols.")
}

dfm <- df2 %>%
  inner_join(map, by = c("gene" = "SYMBOL")) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

cat("Mapped genes:", nrow(dfm), "\n")
write.csv(dfm, file = "mapped_genes.csv", row.names = FALSE)

if (!is.na(padj_col)) {
  dfm$padj <- as.numeric(dfm[[padj_col]])
  ora_genes <- dfm %>% filter(!is.na(padj) & padj < P_ADJ_CUT) %>% pull(ENTREZID)
  if (length(ora_genes) < 10) {
    warning("Significant genes < 10; ORA may be unstable. Falling back to top-|rank| genes (top 300).")
    ora_genes <- dfm %>% arrange(desc(abs(rank))) %>% head(300) %>% pull(ENTREZID)
  }
} else {
  ora_genes <- dfm %>% arrange(desc(abs(rank))) %>% head(300) %>% pull(ENTREZID)
}

gene_rank <- dfm$rank
names(gene_rank) <- dfm$ENTREZID
gene_rank <- sort(gene_rank, decreasing = TRUE)

write.table(ora_genes, file = "ORA_gene_Entrez.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Running ORA (GO/KEGG)...\n")
ego <- enrichGO(
  gene          = ora_genes,
  OrgDb         = orgdb,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ekk <- enrichKEGG(
  gene          = ora_genes,
  organism      = org$kegg,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) write.csv(as.data.frame(ego), "ORA_GO_BP.csv", row.names = FALSE)
if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) write.csv(as.data.frame(ekk), "ORA_KEGG.csv", row.names = FALSE)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  p1 <- dotplot(ego, showCategory = TOPN_ORA) + ggtitle("ORA: GO BP")
  save_gg(p1, "ORA_GO_BP_dotplot.png")
  safe_pdf("ORA_GO_BP_dotplot.pdf"); print(p1); dev.off()

  p2 <- barplot(ego, showCategory = TOPN_ORA) + ggtitle("ORA: GO BP")
  save_gg(p2, "ORA_GO_BP_barplot.png")
  safe_pdf("ORA_GO_BP_barplot.pdf"); print(p2); dev.off()
}

if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
  p3 <- dotplot(ekk, showCategory = TOPN_ORA) + ggtitle("ORA: KEGG")
  save_gg(p3, "ORA_KEGG_dotplot.png")
  safe_pdf("ORA_KEGG_dotplot.pdf"); print(p3); dev.off()

  p4 <- barplot(ekk, showCategory = TOPN_ORA) + ggtitle("ORA: KEGG")
  save_gg(p4, "ORA_KEGG_barplot.png")
  safe_pdf("ORA_KEGG_barplot.pdf"); print(p4); dev.off()
}

cat("Running GSEA (GO/KEGG)...\n")
ggo <- gseGO(
  geneList      = gene_rank,
  OrgDb         = orgdb,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 1e-10,
  verbose       = FALSE
)

gkk <- gseKEGG(
  geneList      = gene_rank,
  organism      = org$kegg,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 1e-10,
  verbose       = FALSE
)

if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) write.csv(as.data.frame(ggo), "GSEA_GO_BP.csv", row.names = FALSE)
if (!is.null(gkk) && nrow(as.data.frame(gkk)) > 0) write.csv(as.data.frame(gkk), "GSEA_KEGG.csv", row.names = FALSE)

if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) {
  p5 <- ridgeplot(ggo, showCategory = TOPN_GSEA) + ggtitle("GSEA: GO BP (ridgeplot)")
  save_gg(p5, "GSEA_GO_BP_ridgeplot.png", w=10, h=7)
  safe_pdf("GSEA_GO_BP_ridgeplot.pdf", width=10, height=7); print(p5); dev.off()

  p6 <- dotplot(ggo, showCategory = TOPN_GSEA, split = ".sign") + facet_grid(.~.sign) + ggtitle("GSEA: GO BP")
  save_gg(p6, "GSEA_GO_BP_dotplot.png", w=10, h=7)
  safe_pdf("GSEA_GO_BP_dotplot.pdf", width=10, height=7); print(p6); dev.off()
}

if (!is.null(gkk) && nrow(as.data.frame(gkk)) > 0) {
  p7 <- ridgeplot(gkk, showCategory = TOPN_GSEA) + ggtitle("GSEA: KEGG (ridgeplot)")
  save_gg(p7, "GSEA_KEGG_ridgeplot.png", w=10, h=7)
  safe_pdf("GSEA_KEGG_ridgeplot.pdf", width=10, height=7); print(p7); dev.off()

  p8 <- dotplot(gkk, showCategory = TOPN_GSEA, split = ".sign") + facet_grid(.~.sign) + ggtitle("GSEA: KEGG")
  save_gg(p8, "GSEA_KEGG_dotplot.png", w=10, h=7)
  safe_pdf("GSEA_KEGG_dotplot.pdf", width=10, height=7); print(p8); dev.off()
}

gly_id <- paste0(org$kegg, "00010")
tca_id <- paste0(org$kegg, "00020")
oxp_id <- paste0(org$kegg, "00190")

focus_ids <- c(gly_id, tca_id, oxp_id)

plot_gsea_if_exists <- function(gsea_obj, pid, label) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  tab <- as.data.frame(gsea_obj)
  if (!pid %in% tab$ID) return(NULL)
  p <- gseaplot2(gsea_obj, geneSetID = pid, title = label)
  fn_png <- paste0("FOCUS_", gsub("[^A-Za-z0-9]+","_",label), "_", pid, ".png")
  fn_pdf <- paste0("FOCUS_", gsub("[^A-Za-z0-9]+","_",label), "_", pid, ".pdf")
  save_gg(p, fn_png, w=9, h=6)
  safe_pdf(fn_pdf, width=9, height=6); print(p); dev.off()
}

plot_gsea_if_exists(gkk, gly_id, "KEGG Glycolysis / Gluconeogenesis (LDHA virtual KO)")
plot_gsea_if_exists(gkk, tca_id, "KEGG TCA cycle (LDHA virtual KO)")
plot_gsea_if_exists(gkk, oxp_id, "KEGG Oxidative phosphorylation (LDHA virtual KO)")

if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
  ekk_tab <- as.data.frame(ekk)
  focus_ora <- ekk_tab %>% filter(ID %in% focus_ids)
  write.csv(focus_ora, "FOCUS_ORA_KEGG_gly_TCA_OXPHOS.csv", row.names = FALSE)
}

if (!is.null(gkk) && nrow(as.data.frame(gkk)) > 0) {
  gkk_tab <- as.data.frame(gkk)
  focus_gsea <- gkk_tab %>% filter(ID %in% focus_ids) %>% arrange(p.adjust)
  write.csv(focus_gsea, "FOCUS_GSEA_KEGG_gly_TCA_OXPHOS.csv", row.names = FALSE)
}

cat("Done. Outputs written to:", OUT_DIR, "\n")
