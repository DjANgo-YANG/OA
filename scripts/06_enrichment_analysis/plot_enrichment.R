rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(clusterProfiler)
})

setwd(project_path("results_root", "auxiliary_enrichment"))

CFG <- list(
  outdir       = project_path("results_root", "auxiliary_enrichment/RESULTS_ENRICH"),
  file_prefix  = "OPLL2",
  date_tag     = format(Sys.Date(), "%Y%m%d"),

  size = list(
    Hallmark_bar = c(w=10, h=12),
    ORA_GO       = c(w=10, h=12),
    ORA_KEGG     = c(w=10, h=12),
    GSEA_GO      = c(w=10, h=18),
    GSEA_KEGG    = c(w=10, h=18),
    Hallmark_hm  = c(w=12, h=12),
    Hallmark_bb  = c(w=12, h=12)
  ),

  pdf = list(engine="cairo", pointsize=10, useDingbats=FALSE),

  n_hallmark_bar = 20,
  n_dotplot      = 20,
  summary_top_pathways = 30,

  contrast_label = "OPLL vs PLL"
)

safe_name <- function(x) gsub("[\\\\/:*?\"<>|\\s]+", "_", x)

make_pdf <- function(dir, stem) {
  file.path(dir, paste0(CFG$file_prefix, "_", CFG$date_tag, "_", stem, ".pdf"))
}

open_pdf <- function(file, w, h) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  w <- suppressWarnings(as.numeric(w))[1]
  h <- suppressWarnings(as.numeric(h))[1]
  if (!is.finite(w) || w <= 0) stop("Invalid PDF width: ", w)
  if (!is.finite(h) || h <= 0) stop("Invalid PDF height: ", h)

  if (CFG$pdf$engine == "cairo") {
    grDevices::cairo_pdf(file = file, width = w, height = h)
    par(ps = CFG$pdf$pointsize)
  } else {
    grDevices::pdf(file = file, width = w, height = h,
                   pointsize = CFG$pdf$pointsize,
                   useDingbats = CFG$pdf$useDingbats)
  }
}
close_pdf <- function() try(grDevices::dev.off(), silent = TRUE)

ct_dirs <- list.dirs(CFG$outdir, full.names = TRUE, recursive = FALSE)
if (length(ct_dirs) == 0) stop("outdir 下没有任何子目录：", CFG$outdir)

for (d in ct_dirs) {
  ct <- basename(d)
  ct_safe <- safe_name(ct)
  message("[PLOT] ", ct)

  fg_file <- file.path(d, "GSEA_Hallmark_fgsea.tsv")
  if (file.exists(fg_file)) {
    fg <- fread(fg_file)
    if (all(c("pathway","NES","padj") %in% colnames(fg))) {
      top_fg <- fg %>%
        as.data.frame() %>%
        dplyr::filter(!is.na(padj), !is.na(NES)) %>%
        dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
        dplyr::slice_head(n = CFG$n_hallmark_bar)

      outp <- make_pdf(d, paste0(ct_safe, "_GSEA_Hallmark_top", CFG$n_hallmark_bar))
      open_pdf(outp, unname(CFG$size$Hallmark_bar["w"]), unname(CFG$size$Hallmark_bar["h"]))
      if (nrow(top_fg) > 0) {
        tt <- ifelse(CFG$contrast_label=="", paste0(ct, ": Hallmark GSEA"), paste0(ct, ": Hallmark GSEA (", CFG$contrast_label, ")"))
        print(
          ggplot(top_fg, aes(x=reorder(pathway, NES), y=NES)) +
            geom_col() + coord_flip() + theme_classic() +
            labs(title=tt, x=NULL, y="NES")
        )
      } else {
        plot.new(); title("No Hallmark GSEA terms.")
      }
      close_pdf()
    }
  }

  go_rds <- file.path(d, "ORA_GO_BP.rds")
  if (file.exists(go_rds)) {
    ego <- readRDS(go_rds)
    ego_df <- try(as.data.frame(ego), silent = TRUE)
    if (!inherits(ego_df,"try-error") && nrow(ego_df) > 0) {
      outp <- make_pdf(d, paste0(ct_safe, "_ORA_GO_BP_dotplot"))
      open_pdf(outp, unname(CFG$size$ORA_GO["w"]), unname(CFG$size$ORA_GO["h"]))
      print(dotplot(ego, showCategory = CFG$n_dotplot) + ggtitle(paste0(ct, ": GO BP ORA")))
      close_pdf()
    }
  }

  kk_rds <- file.path(d, "ORA_KEGG.rds")
  if (file.exists(kk_rds)) {
    ekk <- readRDS(kk_rds)
    ekk_df <- try(as.data.frame(ekk), silent = TRUE)
    if (!inherits(ekk_df,"try-error") && nrow(ekk_df) > 0) {
      outp <- make_pdf(d, paste0(ct_safe, "_ORA_KEGG_dotplot"))
      open_pdf(outp, unname(CFG$size$ORA_KEGG["w"]), unname(CFG$size$ORA_KEGG["h"]))
      print(dotplot(ekk, showCategory = CFG$n_dotplot) + ggtitle(paste0(ct, ": KEGG ORA")))
      close_pdf()
    }
  }

  ggo_rds <- file.path(d, "GSEA_GO_BP_gseGO.rds")
  if (file.exists(ggo_rds)) {
    ggo <- readRDS(ggo_rds)
    ggo_df <- try(as.data.frame(ggo), silent = TRUE)
    if (!inherits(ggo_df,"try-error") && nrow(ggo_df) > 0) {
      outp <- make_pdf(d, paste0(ct_safe, "_GSEA_GO_BP_dotplot"))
      open_pdf(outp, unname(CFG$size$GSEA_GO["w"]), unname(CFG$size$GSEA_GO["h"]))
      print(dotplot(ggo, showCategory = CFG$n_dotplot, split = ".sign") +
              facet_grid(.~.sign) + ggtitle(paste0(ct, ": GO BP GSEA (gseGO)")))
      close_pdf()
    }
  }

  gkk_rds <- file.path(d, "GSEA_KEGG_gseKEGG.rds")
  if (file.exists(gkk_rds)) {
    gkk <- readRDS(gkk_rds)
    gkk_df <- try(as.data.frame(gkk), silent = TRUE)
    if (!inherits(gkk_df,"try-error") && nrow(gkk_df) > 0) {
      outp <- make_pdf(d, paste0(ct_safe, "_GSEA_KEGG_dotplot"))
      open_pdf(outp, unname(CFG$size$GSEA_KEGG["w"]), unname(CFG$size$GSEA_KEGG["h"]))
      print(dotplot(gkk, showCategory = CFG$n_dotplot, split = ".sign") +
              facet_grid(.~.sign) + ggtitle(paste0(ct, ": KEGG GSEA (gseKEGG)")))
      close_pdf()
    }
  }
}

sum_file <- list.files(CFG$outdir, pattern = "SUMMARY_GSEA_Hallmark_all_celltypes.tsv$", full.names = TRUE)
if (length(sum_file) == 1) {
  gsea_all <- fread(sum_file)
  if (all(c("cell_type","pathway","NES","padj") %in% colnames(gsea_all))) {

    gsea_pick <- gsea_all %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(padj), !is.na(NES)) %>%
      dplyr::group_by(pathway) %>%
      dplyr::summarise(min_padj=min(padj, na.rm=TRUE), max_absNES=max(abs(NES), na.rm=TRUE), .groups="drop") %>%
      dplyr::arrange(min_padj, dplyr::desc(max_absNES)) %>%
      dplyr::slice_head(n = CFG$summary_top_pathways)

    pathways_show <- gsea_pick$pathway
    df_show <- gsea_all %>%
      as.data.frame() %>%
      dplyr::filter(pathway %in% pathways_show) %>%
      dplyr::mutate(log10padj = -log10(pmax(padj, 1e-300)))

    if (nrow(df_show) > 0) {

      out1 <- make_pdf(CFG$outdir, "SUMMARY_Hallmark_NES_heatmap")
      open_pdf(out1, unname(CFG$size$Hallmark_hm["w"]), max(unname(CFG$size$Hallmark_hm["h"]), 0.25*length(pathways_show)+2))
      print(
        ggplot(df_show, aes(x=cell_type, y=pathway, fill=NES)) +
          geom_tile() + theme_bw() +
          theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), panel.grid=element_blank()) +
          labs(title="Hallmark GSEA: NES heatmap (across cell types)", x=NULL, y=NULL)
      )
      close_pdf()

      out2 <- make_pdf(CFG$outdir, "SUMMARY_Hallmark_NES_bubble")
      open_pdf(out2, unname(CFG$size$Hallmark_bb["w"]), max(unname(CFG$size$Hallmark_bb["h"]), 0.25*length(pathways_show)+2))
      print(
        ggplot(df_show, aes(x=cell_type, y=pathway)) +
          geom_point(aes(size=log10padj, color=NES), alpha=0.85) +
          theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
          labs(title="Hallmark GSEA: bubble (size=-log10(padj), color=NES)", x=NULL, y=NULL)
      )
      close_pdf()
    }
  }
} else {
  message("[SKIP] No SUMMARY_GSEA_Hallmark_all_celltypes.tsv found.")
}

message("\n[DONE] Plot-only script finished. PDFs regenerated without re-analysis.")
