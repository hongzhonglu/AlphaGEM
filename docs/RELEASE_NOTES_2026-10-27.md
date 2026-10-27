# Runtime maintenance update — 2026-10-27

This update repairs existing reconstruction workflows and packages their runtime
resources. It preserves the Git-tracked reaction pool and the formal legacy
general-to-specific mapping snapshot. The experimental reaction-pool builder,
SQLite annotation catalogue and specificity backend are not included.

## Included changes

- Expose PLM/BLAST homolog thresholds and record input, parameter and resource
  provenance. Validate target embedding caches against FASTA, checkpoint and
  output hashes, and audit reference proteins that have no ordered locus.
- Support explicit Gurobi/COPT selection with optlang-copt 0.2.1; configure
  gapfilling time, memory and MIP-gap limits and validate incumbent solutions.
- Replace complete GPR gene tokens, remove unused reference gene objects, and
  propagate the supplied structure directory through structural alignment.
- Handle Foldseek capabilities, CLEAN interpreter/checkpoint compatibility and
  empty Rhea results. Run the existing Rhea annotation method through an
  independent legacy module without enabling the experimental backend.
- Verify packaged resources and install omitted eggNOG databases through the
  bundled tool's download, taxon-selection and database-building methods.

## Resource bundle

- Filename: `data_source.tar.xz`
- Compressed bytes: **2,753,960,504** (about 2.75 GB / 2.56 GiB).
- Resource files: **21,949**; uncompressed bytes: **9,276,816,548**.
- SHA-256: `c97fff632a12c76522e5cf109bef582e18f6bc76efeb3d417fbafa4083816283`.

The matching checksum is committed at the repository root. The download URL
will be added after the maintainer uploads the file. See
[RESOURCE_INSTALLATION.md](RESOURCE_INSTALLATION.md) for contents and commands.
Large downloadable eggNOG databases, generated indexes, target-species results
and experimental annotation resources are excluded. Required trained models,
reference snapshots and the dedicated DeepProZyme DIAMOND reference are retained.

## Validation

- 109 targeted tests passed in the AlphaGEM environment.
- Syntax checks, shell syntax and CLI help passed.
- A four-protein E. coli smoke run using extracted bundle resources checked
  embedding generation/cache reuse/invalidation, homolog mapping, draft GPRs,
  reaction stoichiometry/bounds and real DIAMOND/tblastn gapfill gene matching.
  The tblastn fixture is a back-translated coding fragment, not a full genome.
- Actual bundled eggNOG APIs passed offline taxon-selection, FASTA-filtering
  and native DIAMOND database checks; no large database downloads were run.
- DeepEC tokenizer loading passed offline in the deepectransformer environment.
- All archive member hashes and the complete XZ stream were verified.

These checks do not constitute a new full-species reconstruction benchmark or
full CLEAN/DeepEC inference run.
