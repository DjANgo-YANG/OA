# Analysis workflow

## Single-cell workflow

| Step | Script | Input | Output |
| --- | --- | --- | --- |
| 1. Create and merge Seurat objects | [build_seurat.R](../scripts/01_seurat_integration/build_seurat.R) | Existing per-sample 10x count matrices in `data/single_cell/10x` | `1.rds` |
| 2. Detect and remove doublets | [remove_doublets.R](../scripts/02_doublet_removal/remove_doublets.R) | `1.rds` | `nodouble.rds` |
| 3. Run QC, normalization, batch correction and clustering | [qc_and_clustering.R](../scripts/03_preprocessing/qc_and_clustering.R) | `nodouble.rds` | `3.rds`, `4.rds` |
| 4. Annotate clusters | [annotate_clusters.R](../scripts/04_cell_annotation/annotate_clusters.R) | `4.rds` | `seurat_object_celltype_0.6.rds` |

### RDS handoffs

These are the default repository-relative paths. Each step creates its output directory and the next step reads the preceding output.

| File | Written to | Next script reads from |
| --- | --- | --- |
| `1.rds` | `results/single_cell/merged/1.rds` | `results/single_cell/merged/1.rds` |
| `nodouble.rds` | `results/single_cell/doublets/nodouble.rds` | `results/single_cell/doublets/nodouble.rds` |
| `4.rds` | `results/single_cell/qc/4.rds` | `results/single_cell/qc/4.rds` |

## Downstream analyses

| Analysis | Scripts | Inputs and outputs |
| --- | --- | --- |
| Trajectory | [Slingshot branch](../scripts/05_trajectory_inference/slingshot), [Monocle 3](../scripts/05_trajectory_inference/monocle3.R), [CytoTRACE2](../scripts/05_trajectory_inference/cytotrace2.R) | Annotated Seurat objects to inferred lineages, pseudotime or developmental-potential estimates |
| Enrichment | [Cell type and pseudotime enrichment](../scripts/06_enrichment_analysis) | Expression, metadata and gene sets to enrichment tables and plots. This folder also contains exploratory virtual perturbation scripts. |
| Pseudobulk | [Pseudobulk scripts](../scripts/07_pseudobulk) | Single-cell counts and metadata to aggregated expression and differential expression |
| Bulk validation | [Bulk validation scripts](../scripts/08_bulk_validation) | Bulk counts or TPM and metadata to Hallmark GSEA, CP decomposition, deconvolution and ROC outputs |
| Gene-set scoring | [Scoring scripts](../scripts/09_gene_set_scoring) | Annotated Seurat objects and gene lists to per-cell scores |
| Virtual perturbation | [KO/OE scripts](../scripts/10_virtual_perturbation) | Seurat objects and target genes to differential regulation and pathway summaries. Several plotting scripts read saved results. |
| Regulatory and core-gene analyses | [TF and target analysis](../scripts/11_regulatory_analysis), [core-gene analysis](../scripts/12_core_gene_analysis) | Regulon activity, pseudotime, bulk expression or KO results to candidate regulators and gene-level enrichment |

## Fig. 3.3 KO workflow

Run these commands from the repository root after generating the scored Seurat object at `results/glycolysis_scoring/seurat_with_LeadingEdge_scores.rds`.

```sh
Rscript "scripts/10_virtual_perturbation/run_ko.R" --all --plan
Rscript "scripts/10_virtual_perturbation/run_ko.R" --all
Rscript "scripts/10_virtual_perturbation/ko_hallmark_gsea.R"
Rscript "scripts/10_virtual_perturbation/cluster_ko_pathways.R"
Rscript "scripts/10_virtual_perturbation/plot_ko_sankey.R"
```

For one target, replace `--all` with its gene symbol, for example `LDHA`. The KO script writes a result RDS, `diffRegulation_full.txt`, and per-target plots to numbered directories such as `15_LDHA敲除` under `results/virtual_ko`. All four scripts use the same root. Set `OAGLY_KO_ROOT` to an absolute path to use a different root. The KO script stops if a selected target directory already exists. After a partial batch, run the missing targets individually or use a fresh root.

The Hallmark script reads the per-target difference tables and writes per-target enrichments under `HALLMARK_GSEA_EnrichmentMap/per_KO_results`. The PCA script reads those enrichments for clustering, PCA, UMAP and the heatmap. The Sankey script reads both the difference tables and the enrichments. Use a separate KO root to regenerate summaries from previously saved inputs.

## Figure-related scripts

| Figure | Script entry points and inputs |
| --- | --- |
| Fig. 3.2, glycolysis scoring and trajectory | [Whole-dataset scoring](../scripts/09_gene_set_scoring/score_glycolysis_leading_edge.R) uses the 80-gene leading edge from the overall OA–HC Hallmark result. [Slingshot](../scripts/05_trajectory_inference/slingshot/fit_slingshot.R) and [pseudotime gene-set analysis](../scripts/05_trajectory_inference/slingshot/plot_gene_set_trends.R) cover the preHTC-to-HTC branch. The Comp2 leading edge supplies the separate 23-gene panel. |
| Fig. 3.3, virtual perturbation | The current manuscript specifies [KO](../scripts/10_virtual_perturbation/run_ko.R) with up to 10,000 highly variable genes, 10 networks and up to 800 sampled cells per network, and [OE v2](../scripts/10_virtual_perturbation/overexpression/run_oe.R) with up to 5,000 highly variable genes, 5 networks and up to 500 sampled cells per network. Both scripts add the target gene if it is outside the selected set. The [earlier OE implementation](../scripts/10_virtual_perturbation/overexpression/run_oe_legacy.R) is retained for reference. The archived COG2 KO table has 5,001 records. Check its run provenance and the saved inputs used for each panel before comparing results. |
| Fig. 3.4, overlap, ROC and LDHA | [Three-way overlap](../scripts/10_virtual_perturbation/intersect_candidate_genes.R), [23-gene ROC](../scripts/08_bulk_validation/roc_23_genes.R) and [LDHA Hallmark GSEA by absolute Z](../scripts/10_virtual_perturbation/ldha_absz_gsea.R). See the commands and inputs below. |

## Fig. 3.4 ROC and LDHA enrichment

**Panel C.** Run the [23-gene ROC script](../scripts/08_bulk_validation/roc_23_genes.R) on `data/bulk/tpm.csv` and the matching `sample_info.csv`. The analysis used 21 HC and 103 OA bulk samples and log₂(TPM + 1) expression, with HC and OA as the control and case levels and pROC direction `"<"`. The [supplementary ROC table](../supplementary/Supplementary_Table_ROC_23_genes.csv) reports AUC and 95% DeLong confidence intervals in descending AUC order. The panel displays the LDHA curve.

**Panel E.** Run `Rscript "scripts/10_virtual_perturbation/ldha_absz_gsea.R"` from the repository root. The script reads `15_LDHA敲除/diffRegulation_full.txt` under the Fig. 3.3 KO root, retains LDHA, keeps the largest signed Z for repeated gene symbols, then ranks by absolute Z. Hallmark enrichment uses `fgseaSimple` with 10,000 permutations, positive-tail scoring, gene-set sizes 10–500 and seed 123. The script selects 15 pathways by nominal P value and displays them in NES order.

Set `OAGLY_KO_ROOT` to the archived KO root when using the original LDHA differential-regulation table. Set `OAGLY_HALLMARK_GMT` to the archived `HALLMARK_GSEA_EnrichmentMap/HALLMARK_msigdb_export.gmt` file for a repeat run with the archived gene sets. Without a GMT override the script uses the Hallmark release available through `msigdbr`. Outputs include the ranked gene list, all pathway results, the top-15 table and PDF in `LDHA_HALLMARK_absZ` under the KO root. To write elsewhere, set `OAGLY_LDHA_OUTDIR`. The archived inputs reproduced the 15 displayed pathways and their enrichment scores. The resampled P values and normalized enrichment scores can differ slightly from the archived table.

[Script manifest](../script_manifest.csv) records each release script's repository path, byte size and SHA-256 checksum.
