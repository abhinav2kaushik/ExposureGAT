#!/usr/bin/env python3
"""ExposureGAT: patient-specific, exposure-associated node and edge prioritization
from paired multi-omics data.

Commands:
  python3 exposureGAT.py demo                      simulate a small dataset and run the full analysis
  python3 exposureGAT.py tcga [--url URL]          reproduce the published TCGA analysis (see TCGA.md)
  python3 exposureGAT.py run --data DIR [options]  run on your own input folder

Each command fits the model with its sensitivity analyses (random initializations, 80% patient
subsets, shuffled-exposure fits) and then writes summary tables and figures to results/<name>/:
  model/   fitted outputs        csv/   summary tables        plots/   figures        logs/   run logs
Downloaded, simulated and intermediate files go to tmp/.

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
R_PACKAGES = ['dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'jsonlite']

PUBLISHED = dict(epochs=600, lr=0.003, hidden=8, initializations=10, subsets=20, permutations=20,
                 top_k=100, scaling='per_gene')
DEMO = dict(epochs=100, lr=0.003, hidden=8, initializations=3, subsets=4, permutations=4,
            top_k=30, scaling='shared')


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


def prepare_demo_data():
    data = ROOT / 'tmp' / 'demo_data'
    if not (data / 'metadata_samples.csv').exists():
        header('Simulating demo dataset')
        subprocess.run([sys.executable, SRC / 'simulation' / 'simulate_data.py', '--out_dir', data, '--seed', '1'], check=True)
    return data


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
        if f.suffix in ('.csv', '.txt', '.json'):
            target = out / 'csv' / step
        elif f.suffix in ('.pdf', '.png'):
            target = out / 'plots' / step
        else:
            continue
        target.mkdir(parents=True, exist_ok=True)
        shutil.move(str(f), target / f.name)


def pipeline(a, data, name):
    out = ROOT / 'results' / name
    if (out / 'model').exists():
        fail(f'results/{name} already exists. Choose another --name, or remove that folder.')
    logs = out / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    work = ROOT / 'tmp' / 'work' / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    (logs / 'settings.txt').write_text(
        f'command: {a.command}\ndata: {data}\nepochs={a.epochs} lr={a.lr} hidden={a.hidden} '
        f'initializations={a.initializations} subsets={a.subsets} permutations={a.permutations}\n'
        f'top_k={a.top_k} exposure_col={a.exposure_col} scaling={a.scaling} seed={a.seed}\n'
        f'date: {datetime.datetime.now():%Y-%m-%d %H:%M}\n')

    steps = 2 if a.no_plots else 5
    header(f'1/{steps} Auditing inputs')
    run_logged([sys.executable, SRC / 'audit_inputs.py', '--data_dir', data, '--out_dir', work / 'audit'], logs / '1_audit.log')
    organize(work / 'audit', out, 'input_audit')

    header(f'2/{steps} Fitting ExposureGAT ({a.initializations} initializations, {a.subsets} subsets, '
           f'{a.permutations} permutations, {a.epochs} epochs)')
    run_logged([sys.executable, '-u', SRC / 'exposuregat_model.py', '--data_dir', data, '--out_dir', out / 'model',
                '--scaling', a.scaling, '--exposure_col', a.exposure_col, '--epochs', a.epochs, '--lr', a.lr,
                '--hidden', a.hidden, '--initializations', a.initializations, '--subsets', a.subsets,
                '--permutations', a.permutations, '--top_k', a.top_k, '--seed', a.seed],
               logs / '2_model.log', echo=True)

    if not a.no_plots:
        header('3/5 Figures and summary tables')
        run_logged([a.rscript, SRC / 'plotting' / 'plot_results.R', data, out / 'model', work / 'figures'], logs / '3_figures.log')
        organize(work / 'figures', out, 'figures')
        header('4/5 Exposure-associated nodes')
        run_logged([a.rscript, SRC / 'analysis' / 'exposure_associated_nodes.R', data, work / 'exposure_associated', out / 'model'],
                   logs / '4_exposure_associated.log')
        organize(work / 'exposure_associated', out, 'exposure_associated')
        header('5/5 Reproducibility across seeds and patient subsets')
        run_logged([a.rscript, SRC / 'analysis' / 'held_out_summary.R', out / 'model', work / 'reproducibility'],
                   logs / '5_reproducibility.log')
        organize(work / 'reproducibility', out, 'reproducibility')

    shutil.rmtree(work, ignore_errors=True)
    print(f'\nDone. Results in results/{name}/')
    print('  model/  fitted outputs   csv/  summary tables   plots/  figures   logs/  run logs')


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
    g.add_argument('--no-plots', dest='no_plots', action='store_true', help='fit only; skip R summaries and figures')
    g.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'), help='Rscript executable (default: %(default)s)')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', metavar='{demo,tcga,run}')
    demo = commands.add_parser('demo', help='simulate a small dataset and run the full analysis')
    add_options(demo, DEMO)
    tcga = commands.add_parser('tcga', help='reproduce the published TCGA analysis (see TCGA.md)')
    tcga.add_argument('--url', default=TCGA_DATA_URL, help='direct download link of the processed TCGA inputs')
    add_options(tcga, PUBLISHED)
    run = commands.add_parser('run', help='run on your own input folder')
    run.add_argument('--data', required=True, help='input folder (see FileDescription.md)')
    add_options(run, PUBLISHED)
    a = parser.parse_args()
    if a.command is None:
        parser.print_help()
        sys.exit(1)

    check_dependencies(a.rscript, not a.no_plots)
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
