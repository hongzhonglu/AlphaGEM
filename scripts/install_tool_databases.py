"""Install omitted eggNOG databases through the bundled tool's own methods."""

import argparse
import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.verify_resources import sha256


TAXA = {'fungi': 'Fungi', 'bacteria': 'Bacteria', 'mammalia': 'Mammalia'}
# eggNOG 5.0 full protein snapshot. Old bundles contained a small partial file.
PROTEINS_SHA256 = 'efc3685ebcbfb1d73359231dbb823cd4cf8bfc80c1ca1b98c7117d1c4674b8f2'


def load_tool(tool, name, data):
    sys.path.insert(0, str(tool))
    spec = importlib.util.spec_from_file_location(f'alphagem_{name}', tool / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.args = SimpleNamespace(simulate=False, force=False)
    # Use the official functions and endpoints over HTTPS, with checked exits.
    module.BASE_URL = module.BASE_URL.replace('http://', 'https://', 1)
    module.run = lambda command: subprocess.run(command, shell=True, check=True, cwd=data)
    return module


def check_sqlite(path):
    with sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True) as database:
        if not database.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1").fetchone():
            raise ValueError(f'No tables in downloaded database: {path}')


def check_diamond(path):
    subprocess.run(['diamond', 'dbinfo', '--db', str(path)], check=True)


def install(root, taxa):
    tool = root / 'tools/eggnog_mapper'
    if not (tool / 'download_eggnog_data.py').is_file():
        raise FileNotFoundError('Extract the AlphaGEM resource bundle before installing tool databases')
    data = tool / 'data'
    data.mkdir(parents=True, exist_ok=True)
    download = load_tool(tool, 'download_eggnog_data', data)
    create = load_tool(tool, 'create_dbs', data)
    from eggnogmapper.version import __DB_VERSION__, __VERSION__
    if (__VERSION__, __DB_VERSION__) != ('2.1.12', '5.0.2'):
        raise ValueError('This installer requires bundled eggNOG-mapper 2.1.12 / database 5.0.2')
    for filename, method in [('eggnog.db', download.download_annotations),
                             ('eggnog.taxa.db', download.download_taxa)]:
        target = data / filename
        if not target.is_file():
            print(f'Downloading {filename} with eggNOG-mapper', flush=True)
            method('.')
        check_sqlite(target)
    missing = []
    for taxon in taxa:
        target = data / f'{taxon}.dmnd'
        if target.is_file():
            check_diamond(target)
        else:
            missing.append(taxon)
    if missing:
        proteins = data / 'e5.proteomes.faa'
        if not proteins.is_file():
            create.download_proteins('.')
        print('Verifying the complete eggNOG source FASTA', flush=True)
        if sha256(proteins) != PROTEINS_SHA256:
            raise ValueError('e5.proteomes.faa is not the pinned complete snapshot. Move the partial/other file aside and rerun; it will not be used to build databases.')
        info = data / 'e5.taxid_info.tsv'
        if not info.is_file():
            create.download_taxid_info('.')
        for taxon in missing:
            with tempfile.TemporaryDirectory(dir=data, prefix=f'{taxon}-build-') as temporary:
                temporary = Path(temporary)
                selected = create.parse_taxa_table(str(info), None, TAXA[taxon])
                if not selected:
                    raise ValueError(f'No taxa selected for {taxon}')
                fasta = temporary / 'proteins.faa'
                create.parse_proteins(str(fasta), str(proteins), selected)
                if not fasta.stat().st_size:
                    raise ValueError(f'No proteins selected for {taxon}')
                pending = temporary / f'{taxon}.dmnd'
                create.create_diamond_db(str(pending), str(fasta))
                check_diamond(pending)
                pending.replace(data / f'{taxon}.dmnd')
    report = {'mapper': __VERSION__, 'database': __DB_VERSION__, 'taxa': taxa,
              'annotation_and_taxonomy_sqlite_readable': True,
              'taxon_databases_checked_by_diamond': True}
    (data / 'alphagem-database-install.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--taxa', nargs='+', choices=sorted(TAXA), default=['fungi'])
    parser.add_argument('--dry-run', action='store_true', help='show the tool methods without downloads')
    args = parser.parse_args()
    if args.dry_run:
        print(json.dumps({'mapper': '2.1.12', 'database': '5.0.2', 'taxa': args.taxa,
                          'download_methods': ['download_annotations', 'download_taxa',
                                               'download_proteins', 'download_taxid_info'],
                          'build_methods': ['parse_taxa_table', 'parse_proteins', 'create_diamond_db'],
                          'complete_proteins_sha256': PROTEINS_SHA256}, indent=2))
        return
    try:
        report = install(args.root.resolve(), list(dict.fromkeys(args.taxa)))
    except (OSError, ValueError, sqlite3.DatabaseError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Tool database installation failed: {error}\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
