# Reproducing the TCGA analysis

Step-by-step commands to download the processed TCGA-LUAD/LUSC inputs, run ExposureGAT with the published
settings, and regenerate all tables and figures. Run every command from the repository root (`ExposureGAT/`).

Processed input data: **https://doi.org/10.7910/DVN/ABBWZG** (Harvard Dataverse; file `ExposureGAT_TCGA_inputs.zip`)

---

## 1. Get the code and install requirements

```bash
git clone https://github.com/abhinav2kaushik/ExposureGAT.git
cd ExposureGAT

python3 -m pip install -r requirements.txt
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite","ggpubr"), repos = "https://cloud.r-project.org")'
```

## 2. Download the input data

The inputs must be unpacked into `tmp/tcga_data/`, so that this file exists:

```text
ExposureGAT/tmp/tcga_data/metadata_samples.csv
```

**Option A: automatic.** `exposureGAT.py` downloads and unpacks the archive itself when given the direct
download link of `ExposureGAT_TCGA_inputs.zip` from the Dataverse record (this also starts the analysis, step 4):

```bash
python3 exposureGAT.py tcga --url "<direct link to ExposureGAT_TCGA_inputs.zip>"
```

**Option B: manual.** Download `ExposureGAT_TCGA_inputs.zip` from https://doi.org/10.7910/DVN/ABBWZG into the
repository root, then:

```bash
mkdir -p tmp
unzip -q ExposureGAT_TCGA_inputs.zip -d tmp/
mv tmp/ExposureGAT_TCGA_inputs tmp/tcga_data
```

Check the download (expected: `patients: 59`):

```bash
python3 -c "import pandas as pd; m = pd.read_csv('tmp/tcga_data/metadata_samples.csv'); print('patients:', m.patient_id.nunique())"
```

## 3. MSigDB gene sets (for the gene-signature step)

The gene sets are included in the repository: the full MSigDB v2026.1 file with **Entrez IDs**
(all collections), zipped because of GitHub's file-size limit:

```text
ExposureGAT/MSigDB/msigdb.v2026.1.Hs.entrez.gmt.zip
```

Nothing needs to be downloaded or unzipped: `python3 exposureGAT.py tcga` uses the `MSigDB/` folder
automatically, and the gene-signature script reads the zipped file directly. (An unzipped
`msigdb.v2026.1.Hs.entrez.gmt` in the same folder also works. Source: https://www.gsea-msigdb.org.)

The 49 signatures used in the paper (smoking-, nicotine- and lung cancer-related; 44 from C2 and 5 from
C5) are listed in `resources/signature_catalog.csv`. If the `MSigDB/` folder is missing, the
gene-signature step is skipped and everything else runs.

## 4. Run the analysis

```bash
python3 exposureGAT.py tcga
```

The run uses the published settings (600 epochs, learning rate 0.003, hidden dimension 8, 10 initializations,
20 patient subsets, 20 exposure permutations, seed 42, weighted node score with gamma 0.10 and top 10%) and
takes several hours on a CPU. To keep it running after closing the terminal:

```bash
nohup python3 exposureGAT.py tcga > tcga_run.log 2>&1 &
tail -f tcga_run.log
```

When it finishes, results are in `results/tcga/`. Where each manuscript result comes from:

| Manuscript item | File |
|---|---|
| Weighted node scores, top 10% genes (main output) | `csv/weighted_node_score/weighted_node_scores.csv` (CpG nodes: `weighted_node_scores_CpG.csv`) |
| Figure 3A (fitted vs recorded exposure rank) | `plots/figures/Fitted_vs_recorded_exposure_rank.*` |
| Figure 3B (top 30 genes, signed contributions) | `plots/weighted_node_score/Top_genes_signed_contributions_heatmap.*`, data `csv/weighted_node_score/top_genes_signed_contributions.csv` |
| Figure 3C (exposure alignment of top 5/10/20%) | `plots/weighted_node_score/Top_gene_exposure_alignment.*`, data `csv/weighted_node_score/top_gene_alignment_bins.csv` |
| Supplementary Figure 6 (training losses) | `plots/figures/TCGA_Supp6_training.*` |
| Supplementary Figures 7–8 (node and edge priorities) | `plots/figures/TCGA_Supp7_*`, `TCGA_Supp8_*` |
| Supplementary Figure 9 (weighted scores of signature genes) | `plots/gene_signatures/Weighted_score_density_*` |
| Supplementary Figure 10 (ranking agreement) | `plots/label_controls/Ranking_agreement_real_vs_shuffled.*` |
| Supplementary Figure 11 (aligned-gene counts; top 5% agreement) | `plots/label_controls/Exposure_aligned_gene_counts.*`, `Fixed_top5_agreement_with_shuffled_fits.*` |
| Supplementary Figures 12–13 (signature representation) | `plots/gene_signatures/Smoking_signatures.*`, `Lung_cancer_signatures.*` |
| Supplementary Data 1 (direction of change, top 10%) | `csv/genes/direction_of_change_top_genes.csv` |
| Gene-signature results (all 49 sets) | `csv/gene_signatures/gene_signature_table.csv` |
| Ranks of smoking-anchor genes and CpGs | `csv/genes/marker_ranks.csv` |
| Held-out agreement and ranking stability | `csv/reproducibility/Held_out_summary.csv` |
| Real- vs shuffled-label statistics | `csv/label_controls/statistical_comparisons.csv` |

To start again from scratch, remove the previous results first: `rm -rf results/tcga`.

## 5. Regenerate tables and figures only (no refitting)

After step 4, every post-fit step can be rerun from `results/tcga/model/` in a few minutes:

```bash
M=results/tcga/model; OUT=results/tcga/rerun; mkdir -p $OUT
python3 src/scoring/weighted_node_score.py --model $M --out $OUT/weighted_node_score
python3 src/analysis/label_controls.py --model $M --weighted $OUT/weighted_node_score --out $OUT/label_controls
Rscript src/plotting/plot_label_controls.R $OUT/label_controls
python3 src/analysis/direction_of_change.py --data tmp/tcga_data --weighted $OUT/weighted_node_score --out $OUT/genes
python3 src/analysis/marker_ranks.py --model $M --weighted $OUT/weighted_node_score --out $OUT/genes \
  --genes 57491,1543,1545,5348,64005 --cpgs cg05575921,cg03636183,cg21161138,cg19816560,cg05951221,cg06125433,cg25648203
python3 src/analysis/gene_signatures.py --weighted $OUT/weighted_node_score --gmt MSigDB --out $OUT/gene_signatures
Rscript src/plotting/plot_gene_signatures.R $OUT/gene_signatures $OUT/weighted_node_score
Rscript src/plotting/plot_results.R tmp/tcga_data $M $OUT/figures
Rscript src/analysis/held_out_summary.R $M $OUT/reproducibility
```

## 6. Optional: training-length check

Reproduces the choice of 600 epochs by comparing real-label and shuffled-label fits up to 800 epochs
(about 30 minutes):

```bash
python3 src/diagnose_exposure.py --data_dir tmp/tcga_data --out_dir results/tcga_training_length \
  --stage opt --lrs 0.003 --epochs 800 --checkpoints 400 500 600 800 --threads 4
```

---

Settings, all options and output descriptions: [README.md](README.md) and [FileDescription.md](FileDescription.md).
