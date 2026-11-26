# OAgly analysis code

R scripts for the osteoarthritis (OA) project. The single-cell workflow starts with existing per-sample 10x count matrices. The repository also includes bulk RNA-seq validation, glycolysis gene-set scoring, and in silico perturbation analyses.

## Repository contents

| Analysis | Start here |
| --- | --- |
| Single-cell setup, doublet removal, QC and annotation | [Seurat setup](scripts/01_seurat_integration/build_seurat.R), [doublet removal](scripts/02_doublet_removal/remove_doublets.R), [QC and clustering](scripts/03_preprocessing/qc_and_clustering.R), [annotation](scripts/04_cell_annotation/annotate_clusters.R) |
| Trajectory inference | [Slingshot](scripts/05_trajectory_inference/slingshot/fit_slingshot.R), [Monocle 3](scripts/05_trajectory_inference/monocle3.R), [CytoTRACE2](scripts/05_trajectory_inference/cytotrace2.R) |
| Bulk and pseudobulk analyses | [Bulk Hallmark GSEA](scripts/08_bulk_validation/bulk_hallmark_gsea.R), [tensor decomposition](scripts/08_bulk_validation/cp_tensor_decomposition.R), [pseudobulk](scripts/07_pseudobulk/pseudobulk.R) |
| Glycolysis gene-set scoring | [Whole-dataset leading-edge scoring](scripts/09_gene_set_scoring/score_glycolysis_leading_edge.R), [pseudotime gene-set analysis](scripts/05_trajectory_inference/slingshot/plot_gene_set_trends.R) |
| Virtual perturbation and LDHA | [Virtual KO](scripts/10_virtual_perturbation/run_ko.R), [virtual OE v2](scripts/10_virtual_perturbation/overexpression/run_oe.R), [cross-dataset overlap](scripts/10_virtual_perturbation/intersect_candidate_genes.R), [23-gene ROC](scripts/08_bulk_validation/roc_23_genes.R), [LDHA Hallmark GSEA by absolute Z](scripts/10_virtual_perturbation/ldha_absz_gsea.R) |

See [workflow and inputs](docs/WORKFLOW.md) for script order, intermediate objects, and the figure-related branches. The [Fig. 3.3 KO workflow](docs/WORKFLOW.md#fig-33-ko-workflow) gives commands for the 23-target batch and its summary plots. [Script manifest](script_manifest.csv) lists the 54 analysis scripts with their repository paths and SHA-256 checksums.

## Supplementary material

[23-gene ROC results](supplementary/Supplementary_Table_ROC_23_genes.csv) lists gene-level AUC and 95% confidence intervals, ranked by descending AUC. The [Fig. 3.4 workflow](docs/WORKFLOW.md#fig-34-roc-and-ldha-enrichment) gives the analysis inputs and output files.

## Getting started

1. Install R and the packages for the analyses you plan to run. From the repository root, run the commands below in R to install the listed CRAN and Bioconductor dependencies. For specialized packages such as DoubletFinder, Monocle 3 and scTenifoldKnk, follow their upstream installation instructions. Check the package imports in each script you use.

   ```r
   source("requirements.R")
   install_oagly_dependencies(install_specialized = TRUE)
   ```

2. Place the inputs under `data/` or set `OAGLY_DATA_ROOT` to your input directory. The Seurat entry script expects one subdirectory per sample under `data/single_cell/10x`, each readable by `Seurat::Read10X()`. Additional path options are listed in [the example configuration](config/paths.example.R). Keep sample identifiers consistent between directory names and downstream metadata.
3. Run individual scripts from the repository root. The first script writes its merged object under `results/single_cell/merged`.

   ```sh
   Rscript scripts/01_seurat_integration/build_seurat.R
   ```

The scripts use a mix of locally generated objects and existing inputs, including bulk count or TPM tables, sample metadata, gene sets and saved perturbation results. These files are supplied separately. Follow the [RDS handoffs](docs/WORKFLOW.md#rds-handoffs) when running the single-cell steps in sequence. See [portability notes](docs/PORTABILITY.md) for current verification and environment details.
