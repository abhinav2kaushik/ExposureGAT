#!/usr/bin/env bash
# =============================================================================
# ExposureGAT: single entry point
# =============================================================================
# Usage:
#   bash exposureGAT.sh demo                         simulate a small dataset and run everything
#   bash exposureGAT.sh tcga [--url URL]             download processed TCGA inputs and reproduce the paper
#   bash exposureGAT.sh run --data DIR [options]     run on your own input folder
#
# Options (run/tcga/demo):
#   --name NAME            results/NAME (default: demo | tcga | name of the data folder)
#   --epochs N             training epochs              (default: 600; demo 100)
#   --lr X                 learning rate                (default: 0.003)
#   --hidden N             hidden dimension             (default: 8)
#   --initializations N    random-seed fits             (default: 10; demo 3)
#   --subsets N            80% patient-subset fits      (default: 20; demo 4)
#   --permutations N       shuffled-exposure fits       (default: 20; demo 4)
#   --top_k N              per-patient top nodes        (default: 100; demo 30)
#   --exposure_col NAME    metadata exposure column     (default: smoking_intensity)
#   --scaling MODE         per_gene | shared            (default: per_gene; demo shared)
#   --seed N               base random seed             (default: 42)
#   --no-plots             skip the R summaries and figures
#   --python PATH / --rscript PATH   interpreters (or set PYTHON / RSCRIPT)
#
# Outputs (results/NAME/):
#   model/   raw fitted outputs from exposureGAT.py (one folder per fit + summary tables)
#   csv/     summary tables, one subfolder per step
#   plots/   figures (PDF and PNG), one subfolder per step
#   logs/    log of every step
# Intermediate and downloaded files go to tmp/.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"
# Processed TCGA inputs (59-patient analysis cohort). Set to the archive's public URL.
TCGA_DATA_URL="${TCGA_DATA_URL:-https://REPLACE-WITH-ARCHIVE-URL/ExposureGAT_TCGA_inputs.zip}"

usage() { sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
die() { echo "ERROR: $*" >&2; exit 1; }
step() { echo; echo "== $* =="; }

MODE="${1:-}"; [ -z "$MODE" ] && usage 1; shift || true
case "$MODE" in demo|tcga|run) ;; -h|--help|help) usage 0 ;; *) die "Unknown command '$MODE' (use demo, tcga or run)";; esac

# ---- defaults per mode --------------------------------------------------------
DATA=""; NAME=""; URL="$TCGA_DATA_URL"; PLOTS=1
EPOCHS=600; LR=0.003; HIDDEN=8; INITS=10; SUBSETS=20; PERMS=20; TOPK=100; SEED=42
EXPOSURE=smoking_intensity; SCALING=per_gene
if [ "$MODE" = demo ]; then EPOCHS=100; INITS=3; SUBSETS=4; PERMS=4; TOPK=30; SCALING=shared; fi

while [ $# -gt 0 ]; do
  case "$1" in
    --data) DATA="$2"; shift 2;;            --name) NAME="$2"; shift 2;;
    --url) URL="$2"; shift 2;;              --epochs) EPOCHS="$2"; shift 2;;
    --lr) LR="$2"; shift 2;;                --hidden) HIDDEN="$2"; shift 2;;
    --initializations) INITS="$2"; shift 2;; --subsets) SUBSETS="$2"; shift 2;;
    --permutations) PERMS="$2"; shift 2;;   --top_k) TOPK="$2"; shift 2;;
    --exposure_col) EXPOSURE="$2"; shift 2;; --scaling) SCALING="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;            --no-plots) PLOTS=0; shift;;
    --python) PYTHON="$2"; shift 2;;        --rscript) RSCRIPT="$2"; shift 2;;
    -h|--help) usage 0;;                    *) die "Unknown option: $1";;
  esac
done

# ---- dependency checks --------------------------------------------------------
"$PYTHON" -c 'import numpy, pandas, scipy, torch, torch_geometric' 2>/dev/null \
  || die "Python packages missing. Run: $PYTHON -m pip install -r requirements.txt"
if [ "$PLOTS" = 1 ]; then
  command -v "$RSCRIPT" >/dev/null 2>&1 || die "Rscript not found (install R or use --no-plots)"
  "$RSCRIPT" -e 'p <- c("dplyr","tidyr","readr","ggplot2","patchwork","jsonlite"); m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]; if (length(m)) stop(paste("Install R packages:", paste(m, collapse = ", ")))' \
    || die "R packages missing (see README)"
fi

mkdir -p "$HERE/tmp" "$HERE/results"

# ---- input data ---------------------------------------------------------------
case "$MODE" in
  demo)
    NAME="${NAME:-demo}"
    DATA="$HERE/tmp/demo_data"
    if [ ! -f "$DATA/metadata_samples.csv" ]; then
      step "Simulating demo dataset"
      "$PYTHON" "$HERE/src/simulation/simulate_data.py" --out_dir "$DATA" --seed 1
    fi;;
  tcga)
    NAME="${NAME:-tcga}"
    DATA="$HERE/tmp/tcga_data"
    if [ ! -f "$DATA/metadata_samples.csv" ]; then
      step "Downloading processed TCGA inputs"
      case "$URL" in *REPLACE-WITH*) die "TCGA archive URL not set. Use --url URL or set TCGA_DATA_URL.";; esac
      curl -fL --retry 3 -o "$HERE/tmp/tcga_inputs.zip" "$URL" || die "Download failed: $URL"
      rm -rf "$HERE/tmp/tcga_unzip"; mkdir -p "$HERE/tmp/tcga_unzip"
      unzip -q "$HERE/tmp/tcga_inputs.zip" -d "$HERE/tmp/tcga_unzip"
      META="$(find "$HERE/tmp/tcga_unzip" -name metadata_samples.csv -not -path '*/__MACOSX/*' | head -n 1)"
      [ -n "$META" ] || die "metadata_samples.csv not found in the downloaded archive"
      rm -rf "$DATA"; mv "$(dirname "$META")" "$DATA"; rm -rf "$HERE/tmp/tcga_unzip"
    fi;;
  run)
    [ -n "$DATA" ] || die "run requires --data DIR"
    [ -f "$DATA/metadata_samples.csv" ] || die "No metadata_samples.csv in $DATA"
    DATA="$(cd "$DATA" && pwd)"
    NAME="${NAME:-$(basename "$DATA")}";;
esac

OUT="$HERE/results/$NAME"
[ -e "$OUT/model" ] && die "results/$NAME already exists. Use --name to choose a new name, or remove that folder."
mkdir -p "$OUT/model" "$OUT/csv" "$OUT/plots" "$OUT/logs"
TMP="$HERE/tmp/$NAME"; rm -rf "$TMP"; mkdir -p "$TMP"
cat > "$OUT/logs/settings.txt" <<EOF
command: $MODE
data: $DATA
epochs=$EPOCHS lr=$LR hidden=$HIDDEN initializations=$INITS subsets=$SUBSETS permutations=$PERMS
top_k=$TOPK exposure_col=$EXPOSURE scaling=$SCALING seed=$SEED
date: $(date)
EOF

# Moves tables to csv/STEP and figures to plots/STEP.
organize() {
  local src="$1" step="$2"
  mkdir -p "$OUT/csv/$step" "$OUT/plots/$step"
  find "$src" -maxdepth 1 -type f \( -name '*.csv' -o -name '*.txt' -o -name '*.json' \) -exec mv {} "$OUT/csv/$step/" \;
  find "$src" -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) -exec mv {} "$OUT/plots/$step/" \;
  rmdir "$OUT/plots/$step" 2>/dev/null || true
}

# ---- 1. input audit -----------------------------------------------------------
step "1/5 Auditing inputs"
"$PYTHON" "$HERE/src/audit_inputs.py" --data_dir "$DATA" --out_dir "$TMP/audit" > "$OUT/logs/1_audit.log" 2>&1 \
  || { cat "$OUT/logs/1_audit.log"; die "Input audit failed"; }
organize "$TMP/audit" input_audit

# ---- 2. model fitting ---------------------------------------------------------
step "2/5 Fitting ExposureGAT ($INITS initializations, $SUBSETS subsets, $PERMS permutations, $EPOCHS epochs)"
"$PYTHON" -u "$HERE/exposureGAT.py" --data_dir "$DATA" --out_dir "$OUT/model" \
  --scaling "$SCALING" --exposure_col "$EXPOSURE" --epochs "$EPOCHS" --lr "$LR" --hidden "$HIDDEN" \
  --initializations "$INITS" --subsets "$SUBSETS" --permutations "$PERMS" --top_k "$TOPK" --seed "$SEED" \
  2>&1 | tee "$OUT/logs/2_model.log"

if [ "$PLOTS" = 0 ]; then
  echo; echo "Done (no plots). Results: results/$NAME"; exit 0
fi

# ---- 3. figures ---------------------------------------------------------------
step "3/5 Figures and summary tables"
"$RSCRIPT" "$HERE/src/plotting/plot_results.R" "$DATA" "$OUT/model" "$TMP/figures" > "$OUT/logs/3_figures.log" 2>&1 \
  || { tail -n 20 "$OUT/logs/3_figures.log"; die "Figure generation failed (see logs/3_figures.log)"; }
organize "$TMP/figures" figures

# ---- 4. exposure-associated nodes ---------------------------------------------
step "4/5 Exposure-associated nodes"
"$RSCRIPT" "$HERE/src/analysis/exposure_associated_nodes.R" "$DATA" "$TMP/exposure_associated" "$OUT/model" \
  > "$OUT/logs/4_exposure_associated.log" 2>&1 \
  || { tail -n 20 "$OUT/logs/4_exposure_associated.log"; die "Exposure-associated analysis failed"; }
organize "$TMP/exposure_associated" exposure_associated

# ---- 5. reproducibility across seeds and splits -------------------------------
step "5/5 Reproducibility across seeds and patient subsets"
"$RSCRIPT" "$HERE/src/analysis/held_out_summary.R" "$OUT/model" "$TMP/reproducibility" \
  > "$OUT/logs/5_reproducibility.log" 2>&1 \
  || { tail -n 20 "$OUT/logs/5_reproducibility.log"; die "Reproducibility summary failed"; }
organize "$TMP/reproducibility" reproducibility

rm -rf "$TMP"
echo
echo "Done. Results in results/$NAME:"
echo "  model/  raw fitted outputs     csv/  summary tables     plots/  figures     logs/  run logs"
