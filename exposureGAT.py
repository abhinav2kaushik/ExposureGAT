#!/usr/bin/env python3
"""ExposureGAT: patient-specific, exposure-associated node and edge prioritization
from paired multi-omics data.

Commands:
  python3 exposureGAT.py demo                      simulate a small dataset and run the full analysis
  python3 exposureGAT.py simulation                reproduce the simulation benchmark (ten datasets)
  python3 exposureGAT.py tcga [--url URL]          reproduce the published TCGA analysis (see TCGA.md)
  python3 exposureGAT.py run --data DIR [options]  run on your own input folder

Each analysis fits the model with its sensitivity analyses (random initializations, 80% patient
subsets, shuffled-exposure fits) and writes results to results/<name>/:
  model/   fitted outputs        csv/   summary tables        plots/   figures        logs/   run logs
The main output is the weighted node score (csv/weighted_node_score/weighted_node_scores.csv), the
exposure-aligned node ranking (gene nodes; other node types, e.g. CpGs, in files with
the type as suffix). Downloaded, simulated and intermediate files go to tmp/.

Use "python3 exposureGAT.py <command> --help" for the options of each command.
"""
import argparse
import datetime
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'src'
# Direct download link of the processed TCGA inputs (ExposureGAT_TCGA_inputs.zip, see TCGA.md).
TCGA_DATA_URL = os.environ.get('TCGA_DATA_URL', '')
R_PACKAGES = ['dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite', 'ggpubr']
# Smoking anchors used in the TCGA analysis (Entrez: AHRR, CYP1A1, CYP1B1, 5348 = FXYD1, 64005 = MYO1G).
TCGA_MARKER_GENES = '57491,1543,1545,5348,64005'
TCGA_MARKER_CPGS = 'cg05575921,cg03636183,cg21161138,cg19816560,cg05951221,cg06125433,cg25648203'

PUBLISHED = dict(epochs=600, lr=0.003, hidden=8, initializations=10, subsets=20, permutations=20,
                 top_k=100, scaling='per_gene')
DEMO = dict(epochs=100, lr=0.003, hidden=8, initializations=3, subsets=4, permutations=4,
            top_k=30, scaling='shared')
SIMULATION = dict(epochs=60, lr=0.003, hidden=8, initializations=10, subsets=20, permutations=20,
                  top_k=100, scaling='shared')


def fail(message):
    sys.exit(f'ERROR: {message}')


def header(text):
    print(f'\n== {text} ==', flush=True)


def run_logged(command, log_file, echo=False):
    """Run a command, writing all output to log_file (and to the screen when echo=True)."""
    with open(log_file, 'w') as log:
        process = subprocess.Popen([str(c) for c in command], stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            log.write(line)
            if echo:
                print(line, end='', flush=True)
        if process.wait() != 0:
            if not echo:
                print(''.join(Path(log_file).read_text().splitlines(True)[-20:]))
            fail(f'step failed; see {log_file.relative_to(ROOT)}')


def check_dependencies(rscript, plots):
    missing = []
    for module in ['numpy', 'pandas', 'scipy', 'torch', 'torch_geometric']:
        try:
            __import__(module)
        except ImportError:
            missing.append(module)
    if missing:
        fail(f'missing Python packages: {", ".join(missing)}. Run: {sys.executable} -m pip install -r requirements.txt')
    if plots:
        if shutil.which(rscript) is None:
            fail(f'"{rscript}" not found. Install R, set --rscript, or use --no-plots.')
        check = ('p <- c(%s); m <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]; '
                 'if (length(m)) { cat(paste(m, collapse = ", ")); quit(status = 1) }') % \
                ', '.join(f'"{p}"' for p in R_PACKAGES)
        result = subprocess.run([rscript, '-e', check], capture_output=True, text=True)
        if result.returncode != 0:
            fail(f'missing R packages: {result.stdout.strip() or "see R output"}. '
                 f'In R: install.packages(c({", ".join(repr(p) for p in R_PACKAGES)}))')


def simulate(data, seed, **size):
    if not (data / 'metadata_samples.csv').exists():
        command = [sys.executable, SRC / 'simulation' / 'simulate_three_set.py', '--out_dir', data, '--seed', seed]
        for key, value in size.items():
            command += [f'--{key}', value]
        subprocess.run([str(c) for c in command], check=True, stdout=subprocess.DEVNULL)
    return data


def prepare_demo_data():
    header('Simulating demo dataset')
    return simulate(ROOT / 'tmp' / 'demo_data', 1, n_patients=60, n_genes=80, n_a=16, n_c=16, degree=6)


def prepare_tcga_data(url):
    data = ROOT / 'tmp' / 'tcga_data'
    if (data / 'metadata_samples.csv').exists():
        return data
    if not url:
        fail('TCGA inputs not found in tmp/tcga_data/. Download them as described in TCGA.md, '
             'or pass the direct download link with --url.')
    header('Downloading processed TCGA inputs')
    archive = ROOT / 'tmp' / 'tcga_inputs.zip'
    archive.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, archive)
    except Exception as error:
        fail(f'download failed ({error}): {url}')
    unpacked = ROOT / 'tmp' / 'tcga_unzip'
    shutil.rmtree(unpacked, ignore_errors=True)
    with zipfile.ZipFile(archive) as z:
        z.extractall(unpacked)
    found = [p for p in unpacked.rglob('metadata_samples.csv') if '__MACOSX' not in p.parts]
    if not found:
        fail('metadata_samples.csv not found in the downloaded archive')
    shutil.rmtree(data, ignore_errors=True)
    shutil.move(str(found[0].parent), data)
    shutil.rmtree(unpacked, ignore_errors=True)
    archive.unlink(missing_ok=True)
    return data


def organize(source, out, step):
    """Move tables to csv/<step>/ and figures to plots/<step>/."""
    for f in sorted(Path(source).glob('*')):
        if not f.is_file():
            continue
        if f.suffix in ('.csv', '.txt', '.json', '.gz', '.xlsx'):
            target = out / 'csv' / step
        elif f.suffix in ('.pdf', '.png'):
            target = out / 'plots' / step
        else:
            continue
        target.mkdir(parents=True, exist_ok=True)
        shutil.move(str(f), target / f.name)


def fit(a, data, out, logs, echo=True):
    run_logged([sys.executable, '-u', SRC / 'exposuregat_model.py', '--data_dir', data, '--out_dir', out / 'model',
                '--scaling', a.scaling, '--exposure_col', a.exposure_col, '--epochs', a.epochs, '--lr', a.lr,
                '--hidden', a.hidden, '--initializations', a.initializations, '--subsets', a.subsets,
                '--permutations', a.permutations, '--top_k', a.top_k, '--seed', a.seed],
               logs / '2_model.log', echo=echo)


def settings(a, data, logs):
    (logs / 'settings.txt').write_text(
        f'command: {a.command}\ndata: {data}\nepochs={a.epochs} lr={a.lr} hidden={a.hidden} '
        f'initializations={a.initializations} subsets={a.subsets} permutations={a.permutations}\n'
        f'top_k={a.top_k} exposure_col={a.exposure_col} scaling={a.scaling} seed={a.seed}\n'
        f'weighted node score: gamma={a.gamma} top_fraction={a.top_fraction}\n'
        f'date: {datetime.datetime.now():%Y-%m-%d %H:%M}\n')


def pipeline(a, data, name):
    out = ROOT / 'results' / name
    if (out / 'model').exists():
        fail(f'results/{name} already exists. Choose another --name, or remove that folder.')
    logs = out / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    work = ROOT / 'tmp' / 'work' / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    settings(a, data, logs)
    py = sys.executable

    header('1 Auditing inputs')
    run_logged([py, SRC / 'audit_inputs.py', '--data_dir', data, '--out_dir', work / 'audit'], logs / '1_audit.log')
    organize(work / 'audit', out, 'input_audit')

    header(f'2 Fitting ExposureGAT ({a.initializations} initializations, {a.subsets} subsets, '
           f'{a.permutations} permutations, {a.epochs} epochs)')
    fit(a, data, out, logs)

    header('3 Weighted node score (main output)')
    command = [py, SRC / 'scoring' / 'weighted_node_score.py', '--model', out / 'model', '--out', work / 'weighted',
               '--gamma', a.gamma, '--top_fraction', a.top_fraction]
    run_logged(command + (['--no-plots'] if a.no_plots else []), logs / '3_weighted_node_score.log', echo=True)
    weighted = out / 'csv' / 'weighted_node_score'
    organize(work / 'weighted', out, 'weighted_node_score')

    header('4 Exposure-label controls (real-label vs shuffled-label fits)')
    run_logged([py, SRC / 'analysis' / 'label_controls.py', '--model', out / 'model', '--weighted', weighted,
                '--out', work / 'controls'], logs / '4_label_controls.log')
    if not a.no_plots and (work / 'controls' / 'statistical_comparisons.csv').exists():
        run_logged([a.rscript, SRC / 'plotting' / 'plot_label_controls.R', work / 'controls'], logs / '4_label_controls_plots.log')
    organize(work / 'controls', out, 'label_controls')

    header('5 Direction of molecular change and marker genes')
    import csv
    with open(data / 'metadata_samples.csv') as handle:
        tissues = {row.get('tissue_type') for row in csv.DictReader(handle)}
    if {a.reference_label, a.transition_label} <= tissues and (data / 'X_meth_gene_features.csv').exists():
        run_logged([py, SRC / 'analysis' / 'direction_of_change.py', '--data', data, '--weighted', weighted, '--out', work / 'direction',
                    '--top_fraction', a.top_fraction, '--reference_label', a.reference_label, '--transition_label', a.transition_label],
                   logs / '5_direction_of_change.log')
    if a.markers or a.marker_cpgs:
        run_logged([py, SRC / 'analysis' / 'marker_ranks.py', '--model', out / 'model', '--weighted', weighted, '--out', work / 'direction',
                    '--genes', a.markers, '--cpgs', a.marker_cpgs], logs / '5_marker_ranks.log')
    organize(work / 'direction', out, 'genes')

    if a.gmt:
        header('6 Gene-signature representation (MSigDB)')
        run_logged([py, SRC / 'analysis' / 'gene_signatures.py', '--weighted', weighted, '--gmt', a.gmt, '--catalog', a.catalog,
                    '--out', work / 'signatures', '--top_fraction', a.top_fraction], logs / '6_gene_signatures.log')
        if not a.no_plots:
            run_logged([a.rscript, SRC / 'plotting' / 'plot_gene_signatures.R', work / 'signatures', weighted],
                       logs / '6_gene_signatures_plots.log')
        organize(work / 'signatures', out, 'gene_signatures')
    else:
        print('\n(Gene-signature step skipped; pass --gmt DIR with MSigDB Entrez GMT files to run it.)')

    if not a.no_plots:
        header('7 Model figures and reproducibility summaries')
        run_logged([a.rscript, SRC / 'plotting' / 'plot_results.R', data, out / 'model', work / 'figures'], logs / '7_figures.log')
        organize(work / 'figures', out, 'figures')
        run_logged([a.rscript, SRC / 'analysis' / 'held_out_summary.R', out / 'model', work / 'reproducibility'],
                   logs / '7_reproducibility.log')
        organize(work / 'reproducibility', out, 'reproducibility')

    shutil.rmtree(work, ignore_errors=True)
    print(f'\nDone. Results in results/{name}/')
    print('  Main output: csv/weighted_node_score/weighted_node_scores.csv')
    print('  model/  fitted outputs   csv/  summary tables   plots/  figures   logs/  run logs')


def simulation(a):
    out = ROOT / 'results' / (a.name or 'simulation')
    if out.exists():
        fail(f'{out.relative_to(ROOT)} already exists. Choose another --name, or remove that folder.')
    out.mkdir(parents=True)
    rows = ['seed,data_dir,analysis_dir,weighted_node_score_dir']
    for seed in a.seeds:
        header(f'Simulated dataset {seed} ({a.patients} patients, {a.genes} genes)')
        data = simulate(ROOT / 'tmp' / 'simulation' / f'seed_{seed}', seed, n_patients=a.patients, n_genes=a.genes,
                        n_a=a.set_a, n_c=a.set_c, degree=a.degree)
        target = out / f'seed_{seed}'
        logs = target / 'logs'; logs.mkdir(parents=True)
        settings(a, data, logs)
        fit(a, data, target, logs, echo=False)
        run_logged([sys.executable, SRC / 'scoring' / 'weighted_node_score.py', '--model', target / 'model',
                    '--out', target / 'weighted_node_score', '--gamma', a.gamma, '--top_fraction', a.top_fraction]
                   + (['--no-plots'] if a.no_plots else []), logs / '3_weighted_node_score.log')
        rows.append(f'{seed},"{data}","{target / "model"}","{target / "weighted_node_score"}"')
    (out / 'manifest.csv').write_text('\n'.join(rows) + '\n')
    header('Figure 2 and simulation supplementary data')
    run_logged([sys.executable, SRC / 'simulation' / 'plot_simulation.py', '--manifest', out / 'manifest.csv',
                '--out', out / 'figures', '--top-k', 50], out / 'figures.log')
    if not a.no_plots:
        first = rows[1].split(',')[0]
        run_logged([a.rscript, SRC / 'simulation' / 'plot_simulation_supplement.R', ROOT / 'tmp' / 'simulation' / f'seed_{first}',
                    out / f'seed_{first}' / 'model', out / 'supplementary_figures', out / 'manifest.csv'], out / 'supplementary_figures.log')
    print(f'\nDone. Results in {out.relative_to(ROOT)}/ (AUROC_summary.csv and topk_recovery.csv in figures/)')


def add_options(p, defaults):
    g = p.add_argument_group('analysis settings')
    g.add_argument('--name', help='results are written to results/NAME')
    g.add_argument('--epochs', type=int, default=defaults['epochs'], help='training epochs per fit (default: %(default)s)')
    g.add_argument('--lr', type=float, default=defaults['lr'], help='learning rate (default: %(default)s)')
    g.add_argument('--hidden', type=int, default=defaults['hidden'], help='hidden dimension (default: %(default)s)')
    g.add_argument('--initializations', type=int, default=defaults['initializations'],
                   help='fits on all patients with different seeds (default: %(default)s)')
    g.add_argument('--subsets', type=int, default=defaults['subsets'],
                   help='fits on random 80%% patient subsets (default: %(default)s)')
    g.add_argument('--permutations', type=int, default=defaults['permutations'],
                   help='fits with shuffled exposure (default: %(default)s)')
    g.add_argument('--top_k', type=int, default=defaults['top_k'], help='top nodes per patient for stability summaries (default: %(default)s)')
    g.add_argument('--exposure_col', default='smoking_intensity', help='exposure column in metadata_samples.csv (default: %(default)s)')
    g.add_argument('--scaling', choices=['per_gene', 'shared'], default=defaults['scaling'], help='feature scaling (default: %(default)s)')
    g.add_argument('--seed', type=int, default=42, help='base random seed (default: %(default)s)')
    g.add_argument('--gamma', type=float, default=0.10, help='weight of paired-transition percentile in the weighted node score (default: %(default)s)')
    g.add_argument('--top_fraction', type=float, default=0.10, help='fraction of nodes of each type reported as top nodes (default: %(default)s)')
    g.add_argument('--gmt', help='folder with MSigDB GMT files (Entrez IDs); msigdb.v2026.1.Hs.entrez.gmt alone is sufficient')
    g.add_argument('--catalog', default=str(ROOT / 'resources' / 'signature_catalog.csv'), help='signature catalog (default: resources/signature_catalog.csv)')
    g.add_argument('--markers', default=defaults.get('markers', ''), help='comma-separated Entrez IDs of marker genes to locate in the rankings')
    g.add_argument('--marker_cpgs', default=defaults.get('marker_cpgs', ''), help='comma-separated CpG probe IDs to locate')
    g.add_argument('--reference_label', default='NT', help='tissue_type of the reference state (default: %(default)s)')
    g.add_argument('--transition_label', default='TP', help='tissue_type of the transition state (default: %(default)s)')
    g.add_argument('--no-plots', dest='no_plots', action='store_true', help='skip figures (tables are still written)')
    g.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'), help='Rscript executable (default: %(default)s)')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', metavar='{demo,simulation,tcga,run}')
    demo = commands.add_parser('demo', help='simulate a small dataset and run the full analysis')
    add_options(demo, DEMO)
    sim = commands.add_parser('simulation', help='reproduce the simulation benchmark (Figure 2)')
    sim.add_argument('--seeds', type=int, nargs='+', default=list(range(1, 11)), help='dataset seeds (default: 1-10)')
    for key, default in [('patients', 250), ('genes', 400), ('set_a', 80), ('set_c', 80), ('degree', 8)]:
        sim.add_argument(f'--{key}', type=int, default=default, help=f'(default: {default})')
    add_options(sim, SIMULATION)
    tcga = commands.add_parser('tcga', help='reproduce the published TCGA analysis (see TCGA.md)')
    tcga.add_argument('--url', default=TCGA_DATA_URL, help='direct download link of the processed TCGA inputs')
    add_options(tcga, dict(PUBLISHED, markers=TCGA_MARKER_GENES, marker_cpgs=TCGA_MARKER_CPGS))
    run = commands.add_parser('run', help='run on your own input folder')
    run.add_argument('--data', required=True, help='input folder (see FileDescription.md)')
    add_options(run, PUBLISHED)
    a = parser.parse_args()
    if a.command is None:
        parser.print_help()
        sys.exit(1)

    check_dependencies(a.rscript, not a.no_plots)
    if a.command == 'simulation':
        simulation(a)
        return
    if a.command == 'tcga' and a.gmt is None and (ROOT / 'MSigDB').is_dir():
        a.gmt = str(ROOT / 'MSigDB')
    if a.command == 'demo':
        data, name = prepare_demo_data(), a.name or 'demo'
    elif a.command == 'tcga':
        data, name = prepare_tcga_data(a.url), a.name or 'tcga'
    else:
        data = Path(a.data).resolve()
        if not (data / 'metadata_samples.csv').exists():
            fail(f'no metadata_samples.csv in {data}')
        name = a.name or data.name
    pipeline(a, data, name)


if __name__ == '__main__':
    main()
