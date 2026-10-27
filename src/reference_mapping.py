"""Validate accession-to-locus evidence without inventing model gene identifiers."""

import json
from pathlib import Path

import pandas as pd


def load_reference_loci(refname):
    path = Path('data_available') / f'{refname}.xlsx'
    if not path.is_file():
        raise FileNotFoundError(f'Reference mapping table missing: {path}; install the matching reference resources')
    table = pd.read_excel(path)
    required = {'Entry', 'Gene Names (ordered locus)'}
    if not required <= set(table.columns):
        raise ValueError(f'Reference mapping table {path} requires columns: {sorted(required)}')
    if table['Entry'].isna().any() or table['Entry'].duplicated().any():
        raise ValueError(f'Reference mapping table {path} has blank or duplicate accessions')
    return {
        accession: (locus.strip() if isinstance(locus, str) else None)
        for accession, locus in zip(table['Entry'], table['Gene Names (ordered locus)'])
    }


def map_reference_pairs(frame, reference_column, target_column, refname, audit_path):
    loci = load_reference_loci(refname)
    missing = set(frame[reference_column]) - loci.keys()
    if missing:
        raise ValueError(f'Reference proteins missing ordered locus mapping: {sorted(missing)}')
    mapped = frame[reference_column].map(loci)
    unavailable = mapped.isna() | mapped.eq('')
    excluded = frame.loc[unavailable, [reference_column, target_column]].drop_duplicates()
    audit_path = Path(audit_path)
    audit_path.write_text(json.dumps({
        'reference': refname, 'reason': 'reference protein has no ordered locus mapping',
        'excluded_pair_count': int(unavailable.sum()),
        'excluded_pairs': [
            {'reference_accession': accession, 'target_gene': target}
            for accession, target in excluded.itertuples(index=False, name=None)
        ],
    }, indent=2) + '\n')
    if unavailable.any():
        print(
            f'Excluded {int(unavailable.sum())} homolog pairs without reference loci; see {audit_path}',
            flush=True,
        )
    result = frame.loc[~unavailable].copy()
    result[reference_column] = mapped.loc[~unavailable]
    return result
