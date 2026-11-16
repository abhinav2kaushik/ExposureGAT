# =============================================================================
# ExposureGAT workflow
# =============================================================================
#   snakemake --cores 8                    cohort analysis (TCGA) + simulation benchmark
#   snakemake --cores 8 cohort             cohort analysis only   -> results/<cohort name>/
#   snakemake --cores 8 simulation         simulation benchmark   -> results/<simulation name>/
#   snakemake --cores 4 --configfile config/demo.yaml     quick test on small simulated data
#   snakemake -n                           dry run: list the steps without running them
# Settings are in config/config.yaml.

import os
import sys

configfile: "config/config.yaml"

C = config["cohort"]
S = config["simulation"]
COHORT = f"results/{C['name']}"
SIM = f"results/{S['name']}"
DATA = C["data_dir"]
SEEDS = [str(s) for s in S["seeds"]]
RUN_SIGNATURES = bool(C["gmt_dir"]) and os.path.isdir(C["gmt_dir"])
RUN_MARKERS = bool(C["marker_genes"] or C["marker_cpgs"])
if C["gmt_dir"] and not RUN_SIGNATURES:
    print(f"Note: gene-set folder '{C['gmt_dir']}' not found; the gene-signature step "
          "(Supplementary Figures 9, 12, 13) is skipped.", file=sys.stderr)


def fit_options(c):
    """Model settings shared by every fit."""
    return (f"--epochs {c['epochs']} --lr {c['lr']} --hidden {c['hidden']} "
            f"--initializations {c['initializations']} --subsets {c['subsets']} "
            f"--permutations {c['permutations']} --top_k {c['top_k']} --seed {c['seed']}")


# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
COHORT_FIGURES = (
    [f"{COHORT}/figures/Figure_3{p}.pdf" for p in "ABC"]
    + [f"{COHORT}/figures/Supplementary_Figure_{n}.pdf" for n in ["06", "07", "08", "10", "11A", "11B"]]
    + ([f"{COHORT}/figures/Supplementary_Figure_{n}.pdf" for n in ["09A", "09B", "12", "13"]] if RUN_SIGNATURES else [])
)
COHORT_TABLES = [
    f"{COHORT}/weighted_node_score/weighted_node_scores.csv",
    f"{COHORT}/label_controls/statistical_comparisons.csv",
    f"{COHORT}/reproducibility/Held_out_summary.csv",
    f"{COHORT}/genes/direction_of_change_top_genes.csv",
] + ([f"{COHORT}/genes/marker_ranks.csv"] if RUN_MARKERS else []) \
  + ([f"{COHORT}/gene_signatures/gene_signature_table.csv"] if RUN_SIGNATURES else [])

SIM_FIGURES = (
    [f"{SIM}/figures/Figure_2{p}.pdf" for p in "BCD"]
    + [f"{SIM}/figures/Supplementary_Figure_{n}.pdf" for n in ["01", "02", "03", "04", "05A", "05B"]]
)
SIM_TABLES = [f"{SIM}/tables/AUROC_summary.csv"] + \
    [f"{SIM}/seed_{s}/weighted_node_score/weighted_node_scores.csv" for s in SEEDS]


rule all:
    input: COHORT_FIGURES + COHORT_TABLES + SIM_FIGURES + SIM_TABLES

rule cohort:
    input: COHORT_FIGURES + COHORT_TABLES

rule simulation:
    input: SIM_FIGURES + SIM_TABLES


# =============================================================================
# Cohort analysis
# =============================================================================
# Unpack the processed TCGA inputs (Harvard Dataverse, see TCGA.md) if they are not there yet.
if C["zip"] and not os.path.exists(f"{DATA}/metadata_samples.csv"):
    rule unpack_inputs:
        input: C["zip"]
        output: f"{DATA}/metadata_samples.csv"
        shell:
            """
            rm -rf {DATA}.unzip && mkdir -p {DATA}.unzip && unzip -q {input} -d {DATA}.unzip
            src=$(dirname "$(find {DATA}.unzip -name metadata_samples.csv -not -path '*__MACOSX*' | head -n 1)")
            mkdir -p {DATA} && mv "$src"/* {DATA}/ && rm -rf {DATA}.unzip
            """

rule audit_inputs:
    input: f"{DATA}/metadata_samples.csv"
    output: directory(f"{COHORT}/input_audit")
    log: f"{COHORT}/logs/audit_inputs.log"
    shell: "python3 src/audit_inputs.py --data_dir {DATA} --out_dir {output} > {log} 2>&1"

rule fit_model:
    # 10 initializations, 20 patient-subset fits and 20 exposure-permutation fits.
    input: f"{DATA}/metadata_samples.csv", f"{COHORT}/input_audit"
    output: directory(f"{COHORT}/model")
    log: f"{COHORT}/logs/fit_model.log"
    params: options=fit_options(C)
    shell:
        "python3 -u src/exposuregat_model.py --data_dir {DATA} --out_dir {output} "
        "--scaling {C[scaling]} --exposure_col {C[exposure_col]} {params.options} > {log} 2>&1"

rule weighted_node_score:
    # Main output, plus Figure 3B (top-gene heatmap) and Figure 3C (exposure alignment).
    input: f"{COHORT}/model"
    output:
        table=f"{COHORT}/weighted_node_score/weighted_node_scores.csv",
        fig_b=f"{COHORT}/figures/Figure_3B.pdf",
        fig_c=f"{COHORT}/figures/Figure_3C.pdf"
    log: f"{COHORT}/logs/weighted_node_score.log"
    params: symbols=f"--symbols {C['gene_symbols']}" if C["gene_symbols"] else ""
    shell:
        "python3 src/scoring/weighted_node_score.py --model {input} --out {COHORT}/weighted_node_score "
        "--gamma {C[gamma]} --top_fraction {C[top_fraction]} --figures {COHORT}/figures {params.symbols} > {log} 2>&1"

rule label_controls:
    # Real-label versus shuffled-label fits: Supplementary Figures 10-11.
    input: model=f"{COHORT}/model", weighted=f"{COHORT}/weighted_node_score/weighted_node_scores.csv"
    output:
        table=f"{COHORT}/label_controls/statistical_comparisons.csv",
        figures=[f"{COHORT}/figures/Supplementary_Figure_{n}.pdf" for n in ["10", "11A", "11B"]]
    log: f"{COHORT}/logs/label_controls.log"
    shell:
        """
        python3 src/analysis/label_controls.py --model {input.model} --weighted {COHORT}/weighted_node_score \
            --out {COHORT}/label_controls > {log} 2>&1
        Rscript src/figures/label_control_figures.R {COHORT}/label_controls {COHORT}/figures >> {log} 2>&1
        """

rule held_out_summary:
    # Held-out agreement and ranking stability across seeds and patient subsets (tables).
    input: f"{COHORT}/model"
    output: f"{COHORT}/reproducibility/Held_out_summary.csv"
    log: f"{COHORT}/logs/held_out_summary.log"
    shell: "Rscript src/analysis/held_out_summary.R {input} {COHORT}/reproducibility > {log} 2>&1"

rule direction_of_change:
    # Measured tumor-normal direction of change of the top genes (table).
    input: f"{COHORT}/weighted_node_score/weighted_node_scores.csv"
    output: f"{COHORT}/genes/direction_of_change_top_genes.csv"
    log: f"{COHORT}/logs/direction_of_change.log"
    shell:
        "python3 src/analysis/direction_of_change.py --data {DATA} --weighted {COHORT}/weighted_node_score "
        "--out {COHORT}/genes --top_fraction {C[top_fraction]} "
        "--reference_label {C[reference_label]} --transition_label {C[transition_label]} > {log} 2>&1"

rule marker_ranks:
    # Positions of prespecified marker genes and CpGs in the rankings (table).
    input: model=f"{COHORT}/model", weighted=f"{COHORT}/weighted_node_score/weighted_node_scores.csv"
    output: f"{COHORT}/genes/marker_ranks.csv"
    log: f"{COHORT}/logs/marker_ranks.log"
    shell:
        "python3 src/analysis/marker_ranks.py --model {input.model} --weighted {COHORT}/weighted_node_score "
        "--out {COHORT}/genes --genes '{C[marker_genes]}' --cpgs '{C[marker_cpgs]}' > {log} 2>&1"

rule gene_signatures:
    # MSigDB signature representation among the top genes: Supplementary Figures 9, 12 and 13.
    input: f"{COHORT}/weighted_node_score/weighted_node_scores.csv"
    output:
        table=f"{COHORT}/gene_signatures/gene_signature_table.csv",
        figures=[f"{COHORT}/figures/Supplementary_Figure_{n}.pdf" for n in ["09A", "09B", "12", "13"]]
    log: f"{COHORT}/logs/gene_signatures.log"
    shell:
        """
        python3 src/analysis/gene_signatures.py --weighted {COHORT}/weighted_node_score --gmt {C[gmt_dir]} \
            --catalog {C[catalog]} --out {COHORT}/gene_signatures --top_fraction {C[top_fraction]} > {log} 2>&1
        Rscript src/figures/gene_signature_figures.R {COHORT}/gene_signatures {COHORT}/weighted_node_score \
            {COHORT}/figures >> {log} 2>&1
        """

rule cohort_figures:
    # Figure 3A and Supplementary Figures 6-8.
    input: f"{COHORT}/model"
    output: [f"{COHORT}/figures/{n}.pdf" for n in ["Figure_3A", "Supplementary_Figure_06", "Supplementary_Figure_07", "Supplementary_Figure_08"]]
    log: f"{COHORT}/logs/cohort_figures.log"
    shell: "Rscript src/figures/cohort_figures.R {DATA} {input} {COHORT}/figures > {log} 2>&1"


# =============================================================================
# Simulation benchmark
# =============================================================================
rule simulate_dataset:
    # One independent simulated dataset per seed (graph, node sets, exposure and noise all resampled).
    output: f"{SIM}/data/seed_{{seed}}/metadata_samples.csv"
    log: f"{SIM}/logs/simulate_seed_{{seed}}.log"
    shell:
        "python3 src/simulation/simulate_three_set.py --out_dir {SIM}/data/seed_{wildcards.seed} --seed {wildcards.seed} "
        "--n_patients {S[patients]} --n_genes {S[genes]} --n_a {S[set_a]} --n_c {S[set_c]} --degree {S[degree]} > {log} 2>&1"

rule fit_simulated:
    input: f"{SIM}/data/seed_{{seed}}/metadata_samples.csv"
    output: directory(f"{SIM}/seed_{{seed}}/model")
    log: f"{SIM}/logs/fit_seed_{{seed}}.log"
    params: options=fit_options(S)
    shell:
        "python3 -u src/exposuregat_model.py --data_dir {SIM}/data/seed_{wildcards.seed} --out_dir {output} "
        "--scaling shared --exposure_col smoking_intensity {params.options} > {log} 2>&1"

rule weighted_simulated:
    input: f"{SIM}/seed_{{seed}}/model"
    output: f"{SIM}/seed_{{seed}}/weighted_node_score/weighted_node_scores.csv"
    log: f"{SIM}/logs/weighted_seed_{{seed}}.log"
    shell:
        "python3 src/scoring/weighted_node_score.py --model {input} --out {SIM}/seed_{wildcards.seed}/weighted_node_score "
        "--gamma {S[gamma]} --top_fraction {S[top_fraction]} > {log} 2>&1"

rule simulation_manifest:
    # One row per dataset; the first seed is the illustrated dataset.
    input: expand(f"{SIM}/seed_{{seed}}/model", seed=SEEDS)
    output: f"{SIM}/manifest.csv"
    run:
        with open(output[0], "w") as f:
            f.write("seed,data_dir,analysis_dir\n")
            for s in SEEDS:
                f.write(f"{s},{os.path.abspath(f'{SIM}/data/seed_{s}')},{os.path.abspath(f'{SIM}/seed_{s}/model')}\n")

rule simulation_tables:
    # AUROCs, top-k recovery, repetition controls and design checks.
    input: f"{SIM}/manifest.csv"
    output: f"{SIM}/tables/AUROC_summary.csv"
    log: f"{SIM}/logs/simulation_tables.log"
    shell: "python3 src/simulation/simulation_summary.py --manifest {input} --out {SIM}/tables --top-k 50 > {log} 2>&1"

rule simulation_figures:
    # Figure 2B-D and Supplementary Figures 1-5.
    input: f"{SIM}/manifest.csv"
    output: SIM_FIGURES
    log: f"{SIM}/logs/simulation_figures.log"
    shell: "Rscript src/figures/simulation_figures.R {input} {SIM}/figures > {log} 2>&1"
