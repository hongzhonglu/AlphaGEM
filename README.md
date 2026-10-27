AlphaGEM
======
_Han, W., Xiao, L., Sun, H. et al. AlphaGEM enables precise genome-scale metabolic modelling by integrating protein structure alignment with deep-learning-based dark metabolism mining. Nat Commun (2026). https://doi.org/10.1038/s41467-026-75549-w_

 >## Platform: Ubuntu only (tested on Ubuntu 20.04/22.04). #

Installation:
======
## Install codes and datas: #
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEM.git AlphaGEM

Use the **2026-10-27 runtime resource bundle** matching
`data_source.tar.xz.sha256` in this checkout. Its new download link will be added
after upload. The [previous Figshare bundle](https://dx.doi.org/10.6084/m9.figshare.30630524)
is an older snapshot and does not match this update.

See [resource extraction, verification and tool database installation](docs/RESOURCE_INSTALLATION.md).
Downloadable eggNOG databases are excluded from the archive and prepared using
the bundled tool's own methods. The existing Rhea annotation rules and reaction
pool remain in use.

## Prepare envrionment: #
It need to create two enviroment as DeepECtransformer uses the old version of pytorch.

* __1.create environment for AlphaGEM:__

      conda env create -f AlphaGEM.yml

* __2.create enviroment for deepectransformer:__
  
      conda env create -f deepecenvironment.yml
  
* __3.activate AlphaGEM environment:__
  
      conda activate AlphaGEM
  
* __4.install PyTorch/CUDA/cuDNN:__
 
   Please install a compatible combination of `PyTorch`, `CUDA`, and `cuDNN` for your system yourself. (Recommended `PyTorch` version: >= 2.4.0)
  
   Use the official selector to obtain the correct command for your OS and CUDA version: [PyTorch Local Installation](https://pytorch.org/get-started/locally/)
  
   >Tested platform for this project: `PyTorch 2.9.0` + `CUDA 13.0` on `Nvidia GeForce RTX 5080` and `RTX 5070ti`.

* __5.install ESM (provides `esm` package used by embeddings):__

      pip install fair-esm==2.0.0

* __6.prepare the source data:__
  
      python setup.py install

* __7.Select and license an optimizer (Gurobi or COPT):__

  Choose either Gurobi or COPT and configure a valid license for your chosen
  optimizer. Both are supported through COBRApy/optlang.

  ### Option A: Gurobi

  **Installation** (in the AlphaGEM environment):

  ```bash
  conda install -n AlphaGEM -c gurobi gurobi=11.0.3
  ```

  **License:** Obtain a valid license following the
  [Gurobi licensing guide](https://www.gurobi.com/documentation/).
  For a license key retrieved with `grbgetkey`, run:

  ```bash
  conda run -n AlphaGEM grbgetkey <YOUR-LICENSE-KEY>
  ```

  **Selection:** Use `--solver gurobi` to require Gurobi.

  ### Option B: COPT

  **Installation** (in the AlphaGEM environment):

  ```bash
  conda run -n AlphaGEM python -m pip install coptpy==8.0.6 optlang-copt==0.2.1
  ```

  These packages are already included in `AlphaGEM.yml`. The COBRApy/optlang
  adapter is [optlang-copt](https://github.com/HanWeishang/optlang-copt) `0.2.1`.
  Version `0.2.1` is available on
  [PyPI](https://pypi.org/project/optlang-copt/0.2.1/).

  **License:** Obtain and configure a valid license following the
  [COPT licensing documentation](https://guide.coap.online/). Installing the
  Python packages alone does not configure a license for AlphaGEM-sized models.

  **Selection:** Use `--solver copt` to require COPT.

  ### Automatic selection and shared solver settings

  `--solver auto` (default) tests Gurobi first, then COPT. Detection solves LP and
  binary MIP problems with more than 2,000 variables to reject missing/invalid or
  size-limited licenses before reconstruction. Use `--solver gurobi` or
  `--solver copt` to require that backend; an unavailable explicit choice stops
  immediately. `ALPHAGEM_SOLVER=gurobi|copt` can also select the backend when the
  CLI uses `auto`; an explicit CLI choice takes precedence. GLPK/CPLEX are not
  automatic fallbacks.

  The main entry point sets COBRApy's global solver before creating any models
  and propagates it through `ALPHAGEM_SOLVER`. Drafts, reference/pool models,
  gapfilling, pFBA and knockout checks all use the selected solver, including
  standalone reaction-addition/gapfilling calls. Each run records the solver and
  package versions in `working/<name>/solver.json`. `src/main.py` forwards to the
  same pipeline and accepts the historical mode/argument spellings.

  `optlang-copt 0.2.1` implements variable-removal, quadratic-objective and
  context-rollback fixes inside the adapter. LP/MIP models delete variables
  directly; current/former QP models use a safe rebuild. Copies and rebuilds
  preserve solver configuration and nondefault native parameters, including
  memory limits and relative MIP gaps. AlphaGEM retains its process-local
  fallback only for `0.2.0`; it does not patch `0.2.1`. Infeasible
  drafts retain their solver status and expose missing primal/objective values
  as `NaN`, so COBRA can represent them and proceed to gapfilling without
  interpreting COPT's no-solution sentinel as growth. If a gapfilling solve hits
  a resource limit, an incumbent is accepted only after finite-value checks and
  growth validation. Missing/invalid solutions stop the workflow; final exports
  require finite positive growth.

  CLEAN and its ESM extraction subprocess use the AlphaGEM Python interpreter.
  For legacy ESM1b checkpoints under PyTorch 2.6+, the extraction child temporarily
  allowlists `argparse.Namespace`; weights-only loading remains enabled. Tool
  failures stop the workflow instead of reusing stale annotation output.

  Use the Foldseek `9.427df8a` and DIAMOND `2.1.10` versions declared in
  `AlphaGEM.yml`. Older Foldseek releases lack the required alignment-normalized
  TM-score/probability fields; older DIAMOND releases lack eggNOG's `--iterate`.
  AlphaGEM checks Foldseek output capabilities before starting the search.

Test:
======
   
   ## Create GEM for _C. albicans_ #
   
   Structure alignment mode
   
      python AlphaGEM.py --name candida --mode structure_alignment --fasta ./working/candida/candida.fasta --refname yeast
      
   PLMsearch mode
   
      python AlphaGEM.py --name candida --mode plmsearch --fasta ./working/candida/candida.fasta --refname yeast

   To rerun both Candida workflows with COPT in independent output directories:

      conda run --no-capture-output -n AlphaGEM python AlphaGEM.py --name candida_copt_plmsearch --mode plmsearch --fasta data_available/candida.fasta --refname yeast --solver copt
      conda run --no-capture-output -n AlphaGEM python AlphaGEM.py --name candida_copt_structure --mode structure_alignment --fasta data_available/candida.fasta --structure struct_data/taryeast/candida --refname yeast --solver copt

   Replace `--solver copt` with `--solver gurobi` to require Gurobi, or omit it
   for automatic detection. These commands rerun homolog search, draft transfer,
   annotation, reaction addition and gapfilling. Target ESM embeddings are reused
   only when their manifest matches the current FASTA, checkpoint and embedding
   file hashes. Missing, legacy or mismatched caches are regenerated. The workflow uses the original Rhea-based annotation rules and the existing reaction pool.
   Structure mode requires the downloaded Candida/reference PDBs.

   Backend regression tests (unavailable licensed backends are reported as skips):

      conda run -n AlphaGEM python -m pytest -q tests/test_solver_config.py tests/test_clean_compat.py tests/test_foldseek_compat.py tests/test_plmsearch_pipeline.py
 
Usage:
======

 ## Parameters: #
  - **--solver**: `auto` (default), `gurobi`, or `copt`.
  - **--mode**: Workflow mode. Default: `structure_alignment`. Choices: `structure_alignment`, `plmsearch`
  - **--plm-threshold**, **--plm-score**: PLMSearch homolog score, fraction from 0 to 1. Default: `0.8`.
  - **--plm-coverage**: Minimum pairwise BLAST alignment coverage of both proteins in a PLM candidate. Default: `0.5`.
  - **--plm-identity**: Minimum pairwise BLAST identity of a PLM candidate, in percent (0–100). Default: `20`.
  - **--bbh-coverage**: Minimum query coverage in each direction of the independent BLAST search. Default: `0.25`.
  - **--bbh-identity**: Minimum identity for the independent BLAST best-hit branch, in percent. Default: `60`.
  - **--bbh-direction**: `direct` (default, includes one-way best hits) or `bidirect` (reciprocal best hits only).
  - **--refname**: Reference species. Choices: `ecoli`, `yeast`, `strco`, `human`, `synechocystis`. Required
  - **--name**: Job/species name used to create `working/<name>`. Required
  - **--fasta**: Target species FASTA file path. Required
  - **--list**: Path to a list file of structures and gene names. Default: empty
  - **--structure**: Directory containing structure files. Default: empty
  - **--cleanuse**: Whether CLEAN has been used. To enable, pass `--cleanuse True`; to disable, omit the flag
  - **--USalignfilter**: USalign filter threshold. Default: 0.5
  - **--TMscore**: Structure alignment filter threshold(foldseek). Default: 0.7
  - **--upTMscore**: Safe TMscore(foldseek). Default: 0.9
  - **--TMscoretrans**: Transporter TMscore filter (foldseek). Default: 0.7
  - **--coverage**: Coverage filter threshold(foldseek). Default: 0.8
  - **--upcoverage**: Safe coverage(foldseek). Default: 0.95
  - **--coveragetrans**: Transporter coverage filter(foldseek). Default: 0.8
  - **--pLDDT**: pLDDT filter threshold. Default: 70
  - **--esp**: Clustering parameter. Default: 1
  - **--grothmedium**: Growth medium. Default: `min`; choices: `min`, `full`
  - **--gapfill-time-limit**: Maximum seconds for each gapfilling MIP solve. Default: `3600`.
  - **--gapfill-soft-mem-limit**: Memory limit per gapfilling MIP solve, in GB (Gurobi soft / COPT hard). Default: `8`.
  - **--gapfill-mip-gap**: Accepted relative MIP gap for gapfilling. Default: `0.01`.

PLMSearch uses the union of two evidence branches: PLM candidates supported by
pairwise BLAST, and an independent proteome BLAST best-hit search. PLM score,
PLM coverage, and both identity cutoffs use strict `>` comparisons; BBH coverage
uses `>=`. BBH retains the historical highest-identity hit selection; `direct`
includes nonreciprocal hits with a covered reverse hit. These options apply only
to homolog transfer in `--mode plmsearch`; they do not change Rhea search or
structure-alignment thresholds.

The defaults above restore the historical `src/main.py` homolog settings. Earlier
`AlphaGEM.py` used conflicting hardcoded cutoffs: PLM `0.95`, coverage `0.9`,
identity `70`, and BLAST coverage `0.8`, identity `70`. Those stringent identity
cutoffs discard many distant fungal homologs. To reproduce that earlier setup,
append:

```bash
--plm-threshold 0.95 --plm-coverage 0.9 --plm-identity 70 \
--bbh-coverage 0.8 --bbh-identity 70 --bbh-direction direct
```

The restored settings are a starting point for homolog discovery, not a
species-independent validation of enzyme specificity. Inspect paralogs, GPR
assignments and pathway support before interpreting added reactions. Each new
PLM run prints its parameters and records input hashes, screening settings,
completion status and homolog counts in `working/<name>/plmsearch_parameters.json`.
BLAST best hits are recomputed against the current FASTAs, preventing reuse of a
coverage-only cache from a different reference or target. Lowering the PLM score
reruns similarity inference; refiltering a previously truncated `0.95` result
cannot recover candidates below that score. Target embeddings have a companion
`<name>_embedding.manifest.json`; changes to the normalized FASTA or checkpoint,
a missing manifest, or a changed embedding file trigger regeneration. Packaged
reference embeddings remain versioned input assets and are hashed in the run
manifest. Use a new `--name` for comparisons to preserve earlier results.
The first run after this update recomputes target embeddings that have no manifest.
Reference tables are checked before ESM loading. An accession missing from the
table is an error; proteins listed with no ordered locus are excluded from PLM
and BLAST homolog transfer and recorded in `*_unmapped_reference.json` files.
Their target genes can still enter the later non-homolog annotation stage.
GPR transfer matches complete gene tokens, and exported drafts remove reference
gene objects that no longer participate in any reaction.
The compatibility entry point `src/main.py` shares these CLI defaults.

`--grothmedium min` runs only the minimal-medium gapfilling pass. The optional
full-medium pass runs only with `--grothmedium full`. Gapfilling uses solver
tolerances of `1e-6`; both backends apply the configurable time and relative
MIP-gap limits above. Memory uses Gurobi `SoftMemLimit` in GB or COPT `MemLimit`
in MB (the CLI GB value is multiplied by 1024). The COPT memory limit is hard;
its behavior is not identical to Gurobi's soft limit. Validation failures at an
optimal solution retry with `1e-9` tolerances. Resource-limit configuration errors
stop the run instead of silently dropping the limits. Presolve is explicitly
enabled so optlang does not automatically repeat an entire solve after a resource
limit. Each call uses one native solve; a separate validation retry receives its
own time budget. Resource-limited incumbents still require growth validation.
Gapfill gene alignments stop on external command failure and publish their
output only after success. BLAST comment lines and valid empty results are
handled explicitly; malformed results stop the run. The bit-score cutoff remains 200.

 ## Examples: #
  - **Structure alignment mode (default)**
  
    > **Pre-requisites for Structure Alignment:**
    > If you perform structure alignment, please provide:
    > 1. A directory containing structure data via the `--structure` parameter.
    > 2. A genome-wide table file via the `--list` parameter. 
    >
    > **List File Format Requirements:**
    > *   The table must contain the following columns: **Entry**, **Entry Name**, **Structure Data**.
    > *   The content of the **Entry** column must strictly match the gene IDs provided in your `--fasta` file.

        python AlphaGEM.py \
          --mode structure_alignment \
          --name my_species \
          --refname yeast \
          --fasta ./my_species.fasta \
          --structure ./structure_dir \
          --list ./my_species.xlsx \
          --cleanuse True \
          --TMscore 0.7 --coverage 0.8 --pLDDT 70

  - **PLMSearch mode**

        python AlphaGEM.py \
          --mode plmsearch \
          --name my_species \
          --refname ecoli \
          --fasta ./my_species.fasta
