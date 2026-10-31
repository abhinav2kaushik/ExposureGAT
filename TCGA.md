# Reproducing the TCGA analysis

Step-by-step commands to download the processed TCGA-LUAD/LUSC inputs, run ExposureGAT with the published
settings, and generate all tables and figures. Run every command from the repository root (`ExposureGAT/`).

Processed input data: **https://doi.org/DOI.XXXXXXXXX** (file `ExposureGAT_TCGA_inputs.zip`)

---

## 1. Get the code and install requirements

```bash
git clone https://github.com/abhinav2kaushik/ExposureGAT.git
cd ExposureGAT

python3 -m pip install -r requirements.txt
Rscript -e 'install.packages(c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite"), repos = "https://cloud.r-project.org")'
```

## 2. Download the input data

The inputs must be unpacked into `tmp/tcga_data/`, so that this file exists:

```text
ExposureGAT/tmp/tcga_data/metadata_samples.csv
```

**Option A: automatic.** `exposureGAT.sh` downloads and unpacks the archive itself (step 3) when given the
direct download link of `ExposureGAT_TCGA_inputs.zip` from the DOI record:

```bash
bash exposureGAT.sh tcga --url "<direct link to ExposureGAT_TCGA_inputs.zip>"
```

**Option B: manual.** Download `ExposureGAT_TCGA_inputs.zip` from https://doi.org/DOI.XXXXXXXXX into the
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

## 3. Run the analysis

```bash
bash exposureGAT.sh tcga
```

With the data already in `tmp/tcga_data/`, no download takes place. The run uses the published settings
(600 epochs, learning rate 0.003, hidden dimension 8, 10 initializations, 20 patient subsets,
20 exposure permutations, seed 42) and takes several hours on a CPU.

To keep it running after closing the terminal:

```bash
nohup bash exposureGAT.sh tcga > tcga_run.log 2>&1 &
tail -f tcga_run.log
```

When it finishes, results are in:

```text
results/tcga/
├── model/     fitted outputs
├── csv/       summary tables
├── plots/     figures (PDF and PNG)
└── logs/      settings and step logs
```

To start again from scratch, remove the previous results first: `rm -rf results/tcga`.

## 4. Regenerate tables and figures only (no refitting)

After step 3, the summaries and figures can be rebuilt from `results/tcga/model/` in a few minutes:

```bash
OUT=results/tcga/replot
Rscript src/plotting/plot_results.R tmp/tcga_data results/tcga/model "$OUT/figures"
Rscript src/analysis/exposure_associated_nodes.R tmp/tcga_data "$OUT/exposure_associated" results/tcga/model
Rscript src/analysis/held_out_summary.R results/tcga/model "$OUT/reproducibility"
```

Node heatmaps (genes and CpGs; signed contributions and node priorities):

```bash
VALUE=signed   Rscript src/analysis/node_heatmaps.R tmp/tcga_data "$OUT/heatmaps" results/tcga/model
VALUE=priority Rscript src/analysis/node_heatmaps.R tmp/tcga_data "$OUT/heatmaps" results/tcga/model
```

## 5. Optional: training-length check

Reproduces the choice of 600 epochs by comparing real-label and shuffled-label fits up to 800 epochs
(about 30 minutes):

```bash
python3 src/diagnose_exposure.py --data_dir tmp/tcga_data --out_dir results/tcga_training_length \
  --stage opt --lrs 0.003 --epochs 800 --checkpoints 400 500 600 800 --threads 4
```

---

Settings, all options and output descriptions: [README.md](README.md) and [FileDescription.md](FileDescription.md).
