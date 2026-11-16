# Reproducing the TCGA analysis

Step-by-step instructions to download the processed TCGA-LUAD/LUSC inputs, run ExposureGAT with the
published settings, and regenerate the manuscript tables and figures. Run every command from the
repository root (`ExposureGAT/`).

Processed input data: **https://doi.org/10.7910/DVN/ABBWZG** (Harvard Dataverse; file `ExposureGAT_TCGA_inputs.zip`)

---

## 1. Get the code and install requirements

```bash
git clone https://github.com/abhinav2kaushik/ExposureGAT.git
cd ExposureGAT

python3 -m pip install -r requirements.txt          # includes Snakemake
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite","ggpubr"), repos = "https://cloud.r-project.org")'
```

## 2. Download the input data

Download `ExposureGAT_TCGA_inputs.zip` from https://doi.org/10.7910/DVN/ABBWZG and place it in the
repository root:

```text
ExposureGAT/ExposureGAT_TCGA_inputs.zip
```

That is all: the workflow unpacks it into `data/tcga/` on the first run. (If you prefer to unpack it
yourself, put the files so that `data/tcga/metadata_samples.csv` exists.)

## 3. MSigDB gene sets

The gene sets are included in the repository: the full MSigDB v2026.1 file with **Entrez IDs**
(all collections), zipped because of GitHub's file-size limit:

```text
ExposureGAT/MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip
```

Nothing needs to be downloaded or unzipped; the gene-signature step reads the zipped file directly.
The 49 signatures used in the paper (smoking-, nicotine- and lung cancer-related; 44 from C2 and 5 from
C5) are listed in `resources/signature_catalog.csv`. Source: https://www.gsea-msigdb.org.

## 4. Run the analysis

Check what will run (dry run), then run it:

```bash
snakemake -n cohort
snakemake --cores 8 cohort
```

The workflow uses the published settings in `config/config.yaml` (600 epochs, learning rate 0.003,
hidden dimension 8, 10 initializations, 20 patient subsets, 20 exposure permutations, seed 42, weighted
node score with gamma 0.10 and top 10%). Model fitting takes several hours on a CPU. To keep it running
after closing the terminal:

```bash
nohup snakemake --cores 8 cohort > tcga_run.log 2>&1 &
tail -f tcga_run.log
```

If the run is interrupted, run the same command again: Snakemake resumes from the last completed
step. Each step writes a log to `results/tcga/logs/`.

## 5. Where each manuscript result comes from

All files are in `results/tcga/`. Figures are in `figures/` as PDF and PNG, with their data as CSV.

| Manuscript item | File |
|---|---|
| Weighted node scores, top 10% genes (main output) | `weighted_node_score/weighted_node_scores.csv` |
| Figure 3A (fitted vs recorded exposure rank) | `figures/Figure_3A.pdf` |
| Figure 3B (top 30 genes, signed contributions) | `figures/Figure_3B.pdf`; data `weighted_node_score/top_genes_signed_contributions.csv` |
| Figure 3C (exposure alignment of top 5/10/20%) | `figures/Figure_3C.pdf`; data `weighted_node_score/top_gene_alignment_bins.csv` |
| Supplementary Figure 6 (training losses) | `figures/Supplementary_Figure_06.pdf` |
| Supplementary Figure 7 (node and edge priority distributions) | `figures/Supplementary_Figure_07.pdf` |
| Supplementary Figure 8 (edge vs node priority) | `figures/Supplementary_Figure_08.pdf` |
| Supplementary Figure 9 (weighted scores of signature genes) | `figures/Supplementary_Figure_09A.pdf`, `09B.pdf` |
| Supplementary Figure 10 (ranking agreement, real vs shuffled labels) | `figures/Supplementary_Figure_10.pdf` |
| Supplementary Figure 11 (aligned-gene counts; fixed top 5% agreement) | `figures/Supplementary_Figure_11A.pdf`, `11B.pdf` |
| Supplementary Figures 12–13 (signature representation) | `figures/Supplementary_Figure_12.pdf`, `13.pdf` |
| Gene-signature results (all 49 sets) | `gene_signatures/gene_signature_table.csv` |
| Real- vs shuffled-label statistics (Holm-adjusted) | `label_controls/statistical_comparisons.csv` |
| Held-out agreement and ranking stability | `reproducibility/Held_out_summary.csv` |
| Direction of tumor–normal change of the top 10% genes | `genes/direction_of_change_top_genes.csv` |
| Ranks of the smoking-anchor genes and CpGs | `genes/marker_ranks.csv` |

Figure 3B labels genes by Entrez ID; to show gene symbols, set `gene_symbols` in `config/config.yaml`
to a CSV with columns `entrez_id,symbol`.

## 6. Regenerate tables and figures only (no refitting)

The fitted models in `results/tcga/model/` are reused automatically. To redraw figures or recompute a
post-fit step, force only that step, for example:

```bash
snakemake --cores 4 cohort --forcerun cohort_figures label_controls gene_signatures
```

To start again from scratch: `rm -rf results/tcga`.

## 7. Optional: training-length check

Reproduces the choice of 600 epochs by comparing real-label and shuffled-label fits up to 800 epochs
(about 30 minutes):

```bash
python3 src/diagnose_exposure.py --data_dir data/tcga --out_dir results/tcga_training_length \
  --stage opt --lrs 0.003 --epochs 800 --checkpoints 400 500 600 800 --threads 4
```

---

Settings and output descriptions: [README.md](README.md) and [FileDescription.md](FileDescription.md).
