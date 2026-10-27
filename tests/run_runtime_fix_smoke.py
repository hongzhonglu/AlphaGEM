"""Real E. coli subset smoke for embedding provenance, mapping and GPR transfer.

Uses installed ESM weights, PLM model, BLAST, DIAMOND and the COPT solver.
This is a small engineering smoke, not a complete species reconstruction.
"""

import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src'), str(ROOT / 'plmsearchtools/plmsearch/plmsearch_util')]

import cobra
import pandas as pd
import torch
from Bio import SeqIO
from Bio.Data import CodonTable

from src import fasta_handle, gapfilling, model_build1, model_build2, plmsearch_embedding
from src.plmsearch_pipeline import PlmsearchSettings, run_homolog_search
from src.solver_config import select_solver


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', default='runtimefix_ecoli_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    args = parser.parse_args()
    os.chdir(ROOT)
    torch.set_num_threads(4)
    name = args.name
    work = ROOT / 'working' / name
    work.mkdir(parents=True, exist_ok=False)
    select_solver('copt')
    # Resource bundles omit generated indexes; prepare this smoke's reference.
    subprocess.run(['makeblastdb', '-in', 'data_available/ecoli.fasta', '-dbtype', 'prot'], check=True)
    accessions = ['A0A385XJN2', 'A0A385XK32', 'P0A7Z0', 'P0DQD7']
    references = {r.id.split('|')[1]: r for r in SeqIO.parse('data_available/ecoli.fasta', 'fasta')}
    source = work / 'input.fasta'
    source.write_text(''.join(f'>smoke_{a}\n{references[a].seq}\n' for a in accessions))
    fasta_handle.handle(source, name)
    statuses = [plmsearch_embedding.embedding_generate(name, nogpu=True)['status']]
    statuses.append(plmsearch_embedding.embedding_generate(name, nogpu=True)['status'])
    output = work / f'{name}_embedding.pkl'
    previous_bytes = output.read_bytes()
    previous_vectors = pickle.loads(previous_bytes)
    # Change a residue while preserving the identifier, then run real inference.
    fasta = work / f'{name}.fasta'
    original = fasta.read_text()
    lines = original.splitlines()
    lines[1] = lines[1][:-1] + ('A' if lines[1][-1] != 'A' else 'G')
    fasta.write_text('\n'.join(lines) + '\n')
    statuses.append(plmsearch_embedding.embedding_generate(name, nogpu=True)['status'])
    assert output.read_bytes() != previous_bytes
    assert statuses == ['generated', 'reused', 'generated']
    with output.open('rb') as handle:
        embeddings = pickle.load(handle)
    changed_id = 'sp|smoke_A0A385XJN2|smoke_A0A385XJN2'
    assert not torch.equal(previous_vectors[changed_id], embeddings[changed_id])
    assert set(embeddings) == {f'sp|smoke_{a}|smoke_{a}' for a in accessions}
    assert all(t.shape == (1280,) and t.device.type == 'cpu' and torch.isfinite(t).all()
               for t in embeddings.values())
    run_homolog_search(name, 'ecoli', PlmsearchSettings())
    matrix = pd.read_excel(work / f'matrix_homolog{name}.xlsx').iloc[:, -2:]
    pairs = {tuple(row) for row in matrix.itertuples(index=False, name=None)}
    required = {('b3036', 'smoke_A0A385XJN2'), ('b4223', 'smoke_A0A385XK32'),
                ('b2914', 'smoke_P0A7Z0')}
    assert required <= pairs, (required - pairs)
    assert 'b303' not in set(matrix.iloc[:, 0]) and 'b422' not in set(matrix.iloc[:, 0])
    assert 'smoke_P0DQD7' not in set(matrix.iloc[:, 1])
    exclusions = {
        stage: json.loads(path.read_text())
        for stage, path in {
            'plm': work / f'{name}_ss_homolog0.8_unmapped_reference.json',
            'bbh': work / f'{name}_bbh_unmapped_reference.json',
        }.items()
    }
    assert all(any(row['reference_accession'] == 'P0DQD7'
                   for row in audit['excluded_pairs']) for audit in exclusions.values())
    model_build1.modelbuild('ecoli-GEM.xml', name)
    model_build2.modelbuild('ecoli-GEM.xml', name)
    draft1 = cobra.io.load_yaml_model(work / f'model1{name}.yml')
    draft2 = cobra.io.load_yaml_model(work / f'model2{name}.yml')
    allowed_genes = {f'smoke_{a}' for a in accessions}
    assert {g.id for g in draft1.genes} <= allowed_genes
    assert {g.id for g in draft2.genes} <= allowed_genes
    assert 'smoke_P0A7Z0' in draft2.genes
    reference = cobra.io.read_sbml_model('models/ecoli-GEM.xml')
    assert not {'b3036', 'b4223'} & {g.id for g in reference.genes}
    for draft in [draft1, draft2]:
        for reaction in draft.reactions:
            ref = reference.reactions.get_by_id(reaction.id)
            assert {m.id: c for m, c in reaction.metabolites.items()} == {
                m.id: c for m, c in ref.metabolites.items()}
            assert reaction.bounds == ref.bounds
    # Run both external alignment routes and attach evidence to a reference reaction.
    accession = 'P0A7Z0'
    SeqIO.write([references[accession]], work / 'ecoli_gap.fasta', 'fasta')
    reaction = sorted(reference.genes.get_by_id('b2914').reactions, key=lambda r: r.id)[0]
    expected_genes = (reaction.gpr.genes - {'b2914'}) | {'smoke_P0A7Z0'}
    alignment_results = {}
    for mode in ['diamond', 'tblastn']:
        if mode == 'tblastn':
            # A controlled coding fragment, back-translated from the actual protein.
            # This verifies parser/tool compatibility; it is not a sampled genome.
            codons = {amino: codon for codon, amino in CodonTable.unambiguous_dna_by_id[11].forward_table.items()}
            coding = ''.join(codons[amino] for amino in str(references[accession].seq))
            (work / f'{name}_non_anno.fasta').write_text(
                f'>sp|smoke_{accession}|smoke_{accession}\n{coding}\n')
        gapfilling.run_gapgene_alignment(name, 'ecoli')
        tiny = cobra.Model('gapgene-smoke')
        tiny.add_reactions([reaction.copy()])
        tiny, unmatched = gapfilling.read_gapgenes(
            tiny, [references[accession].id], [reaction.id], 'ecoli', name)
        assert unmatched == []
        assert tiny.reactions.get_by_id(reaction.id).gpr.genes == expected_genes
        assert {g.id for g in tiny.genes} == expected_genes
        hits = (work / f'outfile{name}.csv').read_text()
        (work / f'gapgene-{mode}.tsv').write_text(hits)
        alignment_results[mode] = {'reaction': reaction.id,
                                   'gpr': tiny.reactions.get_by_id(reaction.id).gene_reaction_rule,
                                   'has_blast_comments': hits.startswith('#')}
    report = {
        'name': name, 'scope': 'E. coli four-protein engineering smoke; not a full reconstruction',
        'solver': tiny.solver.interface.__name__, 'embedding_statuses': statuses,
        'embedding_shapes': {k: list(v.shape) for k, v in embeddings.items()},
        'mapping_pairs': sorted(pairs),
        'unmapped_reference_audits': exclusions,
        'draft1': {'reactions': len(draft1.reactions), 'genes': sorted(g.id for g in draft1.genes)},
        'draft2': {'reactions': len(draft2.reactions), 'genes': sorted(g.id for g in draft2.genes)},
        'stoichiometry_and_bounds_match_reference': True,
        'gapgene_alignments': alignment_results,
        'nucleotide_test_input': 'synthetic coding fragment back-translated from actual E. coli protein',
        'fasta_sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(),
        'artifacts': str(work),
    }
    (work / 'runtime_fix_smoke.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
