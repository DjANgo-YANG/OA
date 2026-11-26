rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(edgeR)
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
})

setwd(project_path("results_root", "auxiliary_enrichment"))

CFG_MAIN <- list(
  ann_rds      = project_path("data_root", "auxiliary/trajectory/annotation/ann.rds"),
  outdir       = project_path("results_root", "auxiliary_enrichment/RESULTS_ENRICH"),
  celltype_col = "cell_type",
  group_col    = "group",
  patient_col  = "patient",

  target_celltypes = NULL,
  force_groups = NULL,

  fdr_cutoff   = 0.05,
  logfc_cutoff = 0.25,

  file_prefix  = "OPLL2",
  date_tag     = format(Sys.Date(), "%Y%m%d"),

  pdf = list(engine = c("cairo", "pdf")[1], pointsize = 10, useDingbats = FALSE),

  size = list(
    GSEA_top20 = c(w = 10, h = 6),
    ORA_GO     = c(w = 10, h = 7),
    ORA_KEGG   = c(w = 10, h = 7),

    GSEA_GO    = c(w = 10, h = 7),
    GSEA_KEGG  = c(w = 10, h = 7)
  )
)

CFG_EXTRA <- list(
  summary_top_pathways = 30,
  curves_top_up   = 2,
  curves_top_down = 2,
  min_gs_genes_in_rank = 10
)

CFG_LINK <- list(

  opll_up_logfc_sign = "pos",
  hallmark_top_n = 5,
  go_top_n = 10,
  kegg_top_n = 10
)

CFG_GSEA2 <- list(

  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  padj_show = 0.05,
  showCategory = 20
)

dir.create(CFG_MAIN$outdir, showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[\\\\/:*?\"<>|\\s]+", "_", x)

make_pdf_path <- function(dir, stem) {
  stem2 <- paste0(CFG_MAIN$file_prefix, "_", CFG_MAIN$date_tag, "_", stem)
  file.path(dir, paste0(stem2, ".pdf"))
}

with_pdf <- function(file, w, h, expr) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  if (CFG_MAIN$pdf$engine == "cairo") {
    grDevices::cairo_pdf(file = file, width = w, height = h, pointsize = CFG_MAIN$pdf$pointsize)
  } else {
    grDevices::pdf(file = file, width = w, height = h,
                   pointsize = CFG_MAIN$pdf$pointsize, useDingbats = CFG_MAIN$pdf$useDingbats)
  }
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  tryCatch(force(expr), error = function(e) { plot.new(); title(paste("Plot error:", conditionMessage(e))) })
  invisible(TRUE)
}

pick_two_groups <- function(x, force_groups = NULL) {
  x <- as.character(x)
  tab <- sort(table(x), decreasing = TRUE)
  lev <- names(tab)
  if (!is.null(force_groups)) {
    if (!all(force_groups %in% lev)) stop("force_groups 不在数据里：", paste(force_groups, collapse = ", "))
    return(force_groups)
  }
  if (all(c("PLL","OPLL") %in% lev)) return(c("PLL","OPLL"))
  if (all(c("HC","OA") %in% lev)) return(c("HC","OA"))
  if (length(lev) < 2) stop("group 水平不足2个：", paste(lev, collapse = ", "))
  lev[1:2]
}

ensure_patient_col <- function(obj, patient_col = "patient") {
  md <- obj@meta.data
  if (patient_col %in% colnames(md)) return(obj)
  candidates <- c("patient","Patient","patient_id","donor","donor_id","sample","sample_id","orig.ident")
  hit <- intersect(candidates, colnames(md))
  if (length(hit) > 0) {
    src <- hit[1]
    obj@meta.data[[patient_col]] <- as.character(obj@meta.data[[src]])
    message("[INFO] meta.data 缺少 '", patient_col, "'，已用列 '", src, "' 复制生成。")
    return(obj)
  }
  cn <- colnames(obj)
  if (any(grepl("_", cn))) obj@meta.data[[patient_col]] <- sub("_.*","",cn) else obj@meta.data[[patient_col]] <- "P1"
  message("[WARN] meta.data 缺少 '", patient_col, "'，已用兜底策略生成（请尽量提供 sample/patient 列）。")
  obj
}

run_edger_de <- function(pb_counts, meta_samp, group_levels) {
  meta_samp$group <- factor(meta_samp$group, levels = group_levels)
  y <- DGEList(counts = pb_counts, group = meta_samp$group)
  keep <- filterByExpr(y, group = meta_samp$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ group, data = meta_samp)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)
  tt
}

make_rank <- function(de_table) {
  p <- de_table$PValue
  p[p <= 1e-300] <- 1e-300
  score <- -log10(p) * sign(de_table$logFC)
  names(score) <- de_table$gene
  score <- sort(score, decreasing = TRUE)
  score[is.finite(score)]
}

map_to_entrez <- function(genes_sym, universe_sym) {
  g1 <- bitr(genes_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  u1 <- bitr(universe_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  list(genes = unique(g1$ENTREZID), universe = unique(u1$ENTREZID))
}

add_scalar_col <- function(df, col, value) {
  df[[col]] <- rep(value, nrow(df))
  df
}

map_symbol_to_entrez <- function(symbols) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(symbols), columns="ENTREZID", keytype="SYMBOL")
  m <- m[!is.na(m$ENTREZID), c("SYMBOL","ENTREZID")]
  m[!duplicated(m),]
}
entrez_to_symbol <- function(entrez_ids) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(as.character(entrez_ids)), columns="SYMBOL", keytype="ENTREZID")
  unique(m$SYMBOL[!is.na(m$SYMBOL) & m$SYMBOL!=""])
}
parse_geneID <- function(x) unique(unlist(strsplit(as.character(x), "/", fixed=TRUE)))

pick_opll_up <- function(de_df) {
  if (CFG_LINK$opll_up_logfc_sign == "pos") {
    de_df %>% dplyr::filter(FDR < CFG_MAIN$fdr_cutoff, logFC >= CFG_MAIN$logfc_cutoff)
  } else {
    de_df %>% dplyr::filter(FDR < CFG_MAIN$fdr_cutoff, logFC <= -CFG_MAIN$logfc_cutoff)
  }
}

make_rank_entrez <- function(ranks_sym) {
  if (is.null(names(ranks_sym)) || length(ranks_sym) == 0) return(NULL)
  sym <- names(ranks_sym)
  mp <- map_symbol_to_entrez(sym)
  if (nrow(mp) == 0) return(NULL)
  df <- data.frame(SYMBOL=sym, score=as.numeric(ranks_sym), stringsAsFactors = FALSE)
  df2 <- merge(df, mp, by="SYMBOL", all=FALSE)
  if (nrow(df2) == 0) return(NULL)

  agg <- tapply(df2$score, df2$ENTREZID, function(v) v[which.max(abs(v))])
  ranks_e <- sort(as.numeric(agg), decreasing = TRUE)
  names(ranks_e) <- names(agg)[order(as.numeric(agg), decreasing = TRUE)]
  ranks_e <- ranks_e[is.finite(ranks_e)]
  ranks_e
}

run_link_report_opllup_one_ct <- function(ct, ct_dir, de, fg, hallmark_list) {
  fg_df <- as.data.frame(fg)
  fg_opll <- fg_df %>%
    dplyr::filter(!is.na(padj), !is.na(NES), padj < 0.05, NES > 0) %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    dplyr::slice_head(n = CFG_LINK$hallmark_top_n)
  if (nrow(fg_opll) == 0) return(invisible(FALSE))

  hallmark_le <- lapply(seq_len(nrow(fg_opll)), function(i) as.character(fg_opll$leadingEdge[[i]]))
  names(hallmark_le) <- fg_opll$pathway

  de_up <- pick_opll_up(de)
  if (nrow(de_up) < 10) return(invisible(FALSE))

  universe_sym <- unique(as.character(de$gene))
  sig_sym <- unique(as.character(de_up$gene))

  map_sig <- map_symbol_to_entrez(sig_sym)
  map_uni <- map_symbol_to_entrez(universe_sym)
  sig_entrez <- unique(map_sig$ENTREZID)
  uni_entrez <- unique(map_uni$ENTREZID)
  if (length(sig_entrez) < 10 || length(uni_entrez) < 100) return(invisible(FALSE))

  ego <- try(enrichGO(gene=sig_entrez, universe=uni_entrez, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                      ont="BP", pAdjustMethod="BH", qvalueCutoff=0.2, readable=FALSE), silent=TRUE)
  ekk <- try(enrichKEGG(gene=sig_entrez, universe=uni_entrez, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.2), silent=TRUE)

  ego_df <- if (!inherits(ego,"try-error")) as.data.frame(ego) else data.frame()
  ekk_df <- if (!inherits(ekk,"try-error")) as.data.frame(ekk) else data.frame()

  go_terms <- character(); go_hit_sym <- character()
  if (nrow(ego_df) > 0) {
    ego_top <- ego_df %>% dplyr::arrange(p.adjust, dplyr::desc(Count)) %>% dplyr::slice_head(n = CFG_LINK$go_top_n)
    go_terms <- ego_top$Description
    go_hit_sym <- entrez_to_symbol(parse_geneID(ego_top$geneID))
  }
  kegg_terms <- character(); kegg_hit_sym <- character()
  if (nrow(ekk_df) > 0) {
    ekk_top <- ekk_df %>% dplyr::arrange(p.adjust, dplyr::desc(Count)) %>% dplyr::slice_head(n = CFG_LINK$kegg_top_n)
    kegg_terms <- ekk_top$Description
    kegg_hit_sym <- entrez_to_symbol(parse_geneID(ekk_top$geneID))
  }

  de_core <- de %>% dplyr::transmute(gene=as.character(gene), logFC=as.numeric(logFC), FDR=as.numeric(FDR), PValue=as.numeric(PValue))

  summary_rows <- list(); core_rows <- list()
  for (pw in names(hallmark_le)) {
    le_sym <- unique(hallmark_le[[pw]])
    inter_go <- intersect(le_sym, go_hit_sym)
    inter_k  <- intersect(le_sym, kegg_hit_sym)
    inter_all <- intersect(inter_go, inter_k)

    core_sym <- inter_all; mode <- "LE∩GO∩KEGG"
    if (length(core_sym) == 0) {
      core_sym <- unique(c(inter_go, inter_k)); mode <- "LE∩GO + LE∩KEGG"
      if (length(core_sym) == 0) { core_sym <- le_sym; mode <- "LE only" }
    }

    core_tbl <- de_core %>% dplyr::filter(gene %in% core_sym) %>% dplyr::arrange(FDR, dplyr::desc(abs(logFC)))
    if (nrow(core_tbl) > 0) {
      core_tbl <- core_tbl %>%
        dplyr::mutate(cell_type=ct, hallmark_pathway=pw, link_mode=mode,
                      in_leadingEdge = gene %in% le_sym,
                      in_GO_hits = gene %in% go_hit_sym,
                      in_KEGG_hits = gene %in% kegg_hit_sym) %>%
        dplyr::select(cell_type, hallmark_pathway, link_mode, gene, logFC, FDR, PValue, in_leadingEdge, in_GO_hits, in_KEGG_hits)
      core_rows[[pw]] <- core_tbl
    }

    fg_row <- fg_opll %>% dplyr::filter(pathway==pw) %>% dplyr::slice_head(n=1)
    summary_rows[[pw]] <- data.frame(
      cell_type=ct, hallmark_pathway=pw, hallmark_padj=fg_row$padj, hallmark_NES=fg_row$NES,
      n_leadingEdge=length(le_sym),
      n_LE_GO=length(inter_go), n_LE_KEGG=length(inter_k), n_LE_GO_KEGG=length(inter_all),
      link_mode_used=mode,
      core_genes=paste(core_sym, collapse=","),
      top_GO_terms=paste(go_terms, collapse=" | "),
      top_KEGG_terms=paste(kegg_terms, collapse=" | "),
      stringsAsFactors=FALSE
    )
  }

  sum_df <- rbindlist(summary_rows, fill=TRUE)
  core_df <- if (length(core_rows)>0) rbindlist(core_rows, fill=TRUE) else data.table()

  out_sum  <- file.path(ct_dir, paste0(CFG_MAIN$file_prefix, "_", CFG_MAIN$date_tag, "_", safe_name(ct), "_LINK_REPORT_OPLLup_summary.tsv"))
  out_core <- file.path(ct_dir, paste0(CFG_MAIN$file_prefix, "_", CFG_MAIN$date_tag, "_", safe_name(ct), "_LINK_REPORT_OPLLup_coregenes.tsv"))
  fwrite(sum_df, out_sum, sep="\t")
  if (nrow(core_df)>0) fwrite(core_df, out_core, sep="\t")
  invisible(TRUE)
}

obj <- readRDS(CFG_MAIN$ann_rds)
stopifnot(inherits(obj, "Seurat"))

md0 <- obj@meta.data
need_base <- c(CFG_MAIN$celltype_col, CFG_MAIN$group_col)
if (!all(need_base %in% colnames(md0))) stop("meta.data 缺少列：", paste(setdiff(need_base, colnames(md0)), collapse=", "))

obj <- ensure_patient_col(obj, CFG_MAIN$patient_col)
md <- obj@meta.data

groups_use <- pick_two_groups(md[[CFG_MAIN$group_col]], CFG_MAIN$force_groups)
message("Using groups: ", paste(groups_use, collapse = " vs "))

cells_keep <- rownames(md)[md[[CFG_MAIN$group_col]] %in% groups_use]
obj <- subset(obj, cells = cells_keep)
md <- obj@meta.data

DefaultAssay(obj) <- "RNA"
counts <- GetAssayData(obj, slot="counts")

msig_h <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, gene_symbol)
hallmark_list <- split(msig_h$gene_symbol, msig_h$gs_name)

ct_all <- sort(unique(as.character(md[[CFG_MAIN$celltype_col]])))
if (!is.null(CFG_MAIN$target_celltypes)) ct_all <- CFG_MAIN$target_celltypes

summary_gsea <- list()
summary_ora  <- list()

for (ct in ct_all) {
  message("\n===== Cell type: ", ct, " =====")
  ct_safe <- safe_name(ct)
  ct_dir  <- file.path(CFG_MAIN$outdir, ct_safe)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  cells_ct <- rownames(md)[as.character(md[[CFG_MAIN$celltype_col]]) == ct]
  if (length(cells_ct) < 50) { message("Skip (too few cells): ", length(cells_ct)); next }

  md_ct <- md[cells_ct, , drop=FALSE]
  counts_ct <- counts[, cells_ct, drop=FALSE]

  md_ct$sample_id <- paste0(md_ct[[CFG_MAIN$patient_col]], "__", md_ct[[CFG_MAIN$group_col]])
  samp_tab <- md_ct %>%
    dplyr::select(sample_id, group=dplyr::all_of(CFG_MAIN$group_col), patient=dplyr::all_of(CFG_MAIN$patient_col)) %>%
    dplyr::distinct()

  sample_ids <- unique(md_ct$sample_id)
  pb <- vapply(sample_ids, function(sid) Matrix::rowSums(counts_ct[, md_ct$sample_id == sid, drop=FALSE]),
               FUN.VALUE = numeric(nrow(counts_ct)))
  pb <- as.matrix(pb); rownames(pb) <- rownames(counts_ct); colnames(pb) <- sample_ids

  meta_samp <- samp_tab %>% dplyr::filter(sample_id %in% colnames(pb)) %>% dplyr::arrange(match(sample_id, colnames(pb)))
  ns <- table(meta_samp$group)
  if (any(ns < 2)) { message("Skip (replicates <2): ", paste(names(ns), ns, sep="=", collapse=", ")); next }

  de <- run_edger_de(pb, meta_samp, groups_use)
  fwrite(de, file.path(ct_dir, "DE_edgeR_pseudobulk.tsv"), sep="\t")

  ranks <- make_rank(de)
  hallmark_list2 <- lapply(hallmark_list, function(gs) intersect(gs, names(ranks)))
  hallmark_list2 <- hallmark_list2[lengths(hallmark_list2) >= 10]

  fg <- fgsea(pathways=hallmark_list2, stats=ranks, minSize=10, maxSize=500, nperm=20000)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  fwrite(fg, file.path(ct_dir, "GSEA_Hallmark_fgsea.tsv"), sep="\t")

  top_fg <- fg %>% as.data.frame() %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n=20)

  gsea_pdf <- make_pdf_path(ct_dir, paste0(ct_safe, "_GSEA_Hallmark_top20"))
  with_pdf(gsea_pdf, CFG_MAIN$size$GSEA_top20["w"], CFG_MAIN$size$GSEA_top20["h"], {
    if (nrow(top_fg)>0) {
      print(ggplot(top_fg, aes(x=reorder(pathway,NES), y=NES)) +
              geom_col() + coord_flip() + theme_classic() +
              labs(title=paste0(ct, ": Hallmark GSEA (", groups_use[2], " vs ", groups_use[1], ")"), x=NULL, y="NES"))
    } else { plot.new(); title("No Hallmark GSEA terms.") }
  })

  ranks_entrez <- make_rank_entrez(ranks)
  if (!is.null(ranks_entrez) && length(ranks_entrez) >= 200) {

    ggo <- try(
      gseGO(
        geneList      = ranks_entrez,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        minGSSize     = CFG_GSEA2$minGSSize,
        maxGSSize     = CFG_GSEA2$maxGSSize,
        pvalueCutoff  = CFG_GSEA2$pvalueCutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(ggo, "try-error")) {
      ggo_df <- as.data.frame(ggo)
      saveRDS(ggo, file.path(ct_dir, "GSEA_GO_BP_gseGO.rds"))
      fwrite(ggo_df, file.path(ct_dir, "GSEA_GO_BP_gseGO_table.tsv"), sep="\t")

      if (nrow(ggo_df) > 0) {
        go_gsea_pdf <- make_pdf_path(ct_dir, paste0(ct_safe, "_GSEA_GO_BP_dotplot"))
        with_pdf(go_gsea_pdf, CFG_MAIN$size$GSEA_GO["w"], CFG_MAIN$size$GSEA_GO["h"], {
          print(
            dotplot(ggo, showCategory = CFG_GSEA2$showCategory, split = ".sign") +
              facet_grid(.~.sign) +
              ggtitle(paste0(ct, ": GO BP GSEA (gseGO)"))
          )
        })
      }
    }

    gkk <- try(
      gseKEGG(
        geneList      = ranks_entrez,
        organism      = "hsa",
        minGSSize     = CFG_GSEA2$minGSSize,
        maxGSSize     = CFG_GSEA2$maxGSSize,
        pvalueCutoff  = CFG_GSEA2$pvalueCutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE
      ),
      silent = TRUE
    )

    if (!inherits(gkk, "try-error")) {
      gkk_df <- as.data.frame(gkk)
      saveRDS(gkk, file.path(ct_dir, "GSEA_KEGG_gseKEGG.rds"))
      fwrite(gkk_df, file.path(ct_dir, "GSEA_KEGG_gseKEGG_table.tsv"), sep="\t")

      if (nrow(gkk_df) > 0) {
        kk_gsea_pdf <- make_pdf_path(ct_dir, paste0(ct_safe, "_GSEA_KEGG_dotplot"))
        with_pdf(kk_gsea_pdf, CFG_MAIN$size$GSEA_KEGG["w"], CFG_MAIN$size$GSEA_KEGG["h"], {
          print(
            dotplot(gkk, showCategory = CFG_GSEA2$showCategory, split = ".sign") +
              facet_grid(.~.sign) +
              ggtitle(paste0(ct, ": KEGG GSEA (gseKEGG)"))
          )
        })
      }
    }
  } else {
    message("GSEA(GO/KEGG) skipped: ranks_entrez too small after mapping (need ~>=200).")
  }

  sig <- de %>% dplyr::filter(FDR < CFG_MAIN$fdr_cutoff, abs(logFC) >= CFG_MAIN$logfc_cutoff) %>%
    dplyr::pull(gene) %>% unique()
  universe <- rownames(pb)[Matrix::rowSums(pb)>0]

  if (length(sig) >= 10) {
    mapped <- try(map_to_entrez(sig, universe), silent=TRUE)
    if (!inherits(mapped,"try-error") && length(mapped$genes)>=10 && length(mapped$universe)>=100) {
      ego <- enrichGO(gene=mapped$genes, universe=mapped$universe, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                      ont="BP", pAdjustMethod="BH", qvalueCutoff=0.2, readable=TRUE)
      ekk <- enrichKEGG(gene=mapped$genes, universe=mapped$universe, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.2)

      saveRDS(ego, file.path(ct_dir, "ORA_GO_BP.rds"))
      saveRDS(ekk, file.path(ct_dir, "ORA_KEGG.rds"))

      ego_df <- as.data.frame(ego); ekk_df <- as.data.frame(ekk)
      fwrite(ego_df, file.path(ct_dir, "ORA_GO_BP_table.tsv"), sep="\t")
      fwrite(ekk_df, file.path(ct_dir, "ORA_KEGG_table.tsv"), sep="\t")

      go_pdf <- make_pdf_path(ct_dir, paste0(ct_safe, "_ORA_GO_BP_dotplot"))
      with_pdf(go_pdf, CFG_MAIN$size$ORA_GO["w"], CFG_MAIN$size$ORA_GO["h"], {
        if (nrow(ego_df)>0) print(dotplot(ego, showCategory=20) + ggtitle(paste0(ct, ": GO BP ORA"))) else { plot.new(); title("No GO BP.") }
      })

      kegg_pdf <- make_pdf_path(ct_dir, paste0(ct_safe, "_ORA_KEGG_dotplot"))
      with_pdf(kegg_pdf, CFG_MAIN$size$ORA_KEGG["w"], CFG_MAIN$size$ORA_KEGG["h"], {
        if (nrow(ekk_df)>0) print(dotplot(ekk, showCategory=20) + ggtitle(paste0(ct, ": KEGG ORA"))) else { plot.new(); title("No KEGG.") }
      })

      if (nrow(ego_df)>0) { ego_df <- add_scalar_col(ego_df, "cell_type", ct); ego_df <- add_scalar_col(ego_df, "db", "GO_BP"); summary_ora[[paste0(ct,"_GO")]] <- ego_df }
      if (nrow(ekk_df)>0) { ekk_df <- add_scalar_col(ekk_df, "cell_type", ct); ekk_df <- add_scalar_col(ekk_df, "db", "KEGG"); summary_ora[[paste0(ct,"_KEGG")]] <- ekk_df }
    }
  } else {
    message("ORA skipped: significant genes < 10.")
  }

  fg2 <- fg %>% as.data.frame()
  fg2 <- add_scalar_col(fg2, "cell_type", ct)
  summary_gsea[[ct]] <- fg2

  run_link_report_opllup_one_ct(ct, ct_dir, de, fg, hallmark_list)
}

if (length(summary_gsea)>0) {
  gsea_all <- rbindlist(summary_gsea, fill=TRUE)
  fwrite(gsea_all, file.path(CFG_MAIN$outdir, paste0(CFG_MAIN$file_prefix,"_",CFG_MAIN$date_tag,"_SUMMARY_GSEA_Hallmark_all_celltypes.tsv")), sep="\t")
}
if (length(summary_ora)>0) {
  ora_all <- rbindlist(summary_ora, fill=TRUE)
  fwrite(ora_all, file.path(CFG_MAIN$outdir, paste0(CFG_MAIN$file_prefix,"_",CFG_MAIN$date_tag,"_SUMMARY_ORA_GO_KEGG_all_celltypes.tsv")), sep="\t")
}

message("\n[OK] Main enrichment DONE: ", CFG_MAIN$outdir)

open_pdf <- function(file, w, h) {
  dir.create(dirname(file), showWarnings=FALSE, recursive=TRUE)
  if (CFG_MAIN$pdf$engine=="cairo") grDevices::cairo_pdf(file, w, h, pointsize=CFG_MAIN$pdf$pointsize)
  else grDevices::pdf(file, w, h, pointsize=CFG_MAIN$pdf$pointsize, useDingbats=CFG_MAIN$pdf$useDingbats)
}
close_pdf <- function() try(grDevices::dev.off(), silent=TRUE)
make_pdf <- function(stem) file.path(CFG_MAIN$outdir, paste0(CFG_MAIN$file_prefix,"_",CFG_MAIN$date_tag,"_",stem,".pdf"))

ct_dirs2 <- list.dirs(CFG_MAIN$outdir, full.names=TRUE, recursive=FALSE)
ct_dirs2 <- ct_dirs2[file.exists(file.path(ct_dirs2, "GSEA_Hallmark_fgsea.tsv"))]

collect_gsea2 <- list()
for (d in ct_dirs2) {
  ct <- basename(d)
  x <- fread(file.path(d, "GSEA_Hallmark_fgsea.tsv"))
  if (!all(c("pathway","NES","padj") %in% colnames(x))) next
  x$cell_type <- ct
  collect_gsea2[[ct]] <- x
}
gsea_all2 <- rbindlist(collect_gsea2, fill=TRUE)

gsea_pick <- gsea_all2 %>% as.data.frame() %>%
  dplyr::filter(!is.na(padj), !is.na(NES)) %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(min_padj=min(padj), max_absNES=max(abs(NES)), .groups="drop") %>%
  dplyr::arrange(min_padj, dplyr::desc(max_absNES)) %>%
  dplyr::slice_head(n=CFG_EXTRA$summary_top_pathways)

pathways_show <- gsea_pick$pathway
df_show <- gsea_all2 %>% as.data.frame() %>%
  dplyr::filter(pathway %in% pathways_show) %>%
  dplyr::mutate(log10padj = -log10(pmax(padj, 1e-300)))

if (nrow(df_show)>0) {
  pdf1 <- make_pdf("SUMMARY_Hallmark_NES_heatmap")
  open_pdf(pdf1, 12, max(6, 0.25*length(pathways_show)+2))
  print(
    ggplot(df_show, aes(x=cell_type, y=pathway, fill=NES)) +
      geom_tile() + theme_bw() +
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), panel.grid=element_blank()) +
      labs(title="Hallmark GSEA: NES heatmap", x=NULL, y=NULL)
  )
  close_pdf()

  pdf2 <- make_pdf("SUMMARY_Hallmark_NES_bubble")
  open_pdf(pdf2, 12, max(6, 0.25*length(pathways_show)+2))
  print(
    ggplot(df_show, aes(x=cell_type, y=pathway)) +
      geom_point(aes(size=log10padj, color=NES), alpha=0.85) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
      labs(title="Hallmark GSEA: bubble (size=-log10padj, color=NES)", x=NULL, y=NULL)
  )
  close_pdf()
}

message("\n[OK] Extra plots DONE.")
