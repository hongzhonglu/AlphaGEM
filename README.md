AlphaGEM
======
_AlphaGEM Enables Precise Genome-Scale Metabolic Modelling by Integrating Protein Structure Alignment with deep-learning-based Dark Metabolism Mining_
_doi: https://doi.org/10.1101/2025.07.21.665674_
 >## Platform: Ubuntu only (tested on Ubuntu 20.04/22.04). #

 
Installation:
======
## Install codes and datas: #
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEMs.git AlphaGEM

You can download supplyment resources from figshare:

[Data Source Figshare](https://dx.doi.org/10.6084/m9.figshare.30630524)

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

* __7.Gurobi license (if you don’t have one yet)__

  We recommend obtaining a valid Gurobi license. Otherwise, adjust the code to use a different optimizer. To retrieve your license:

      grbgetkey <YOUR-LICENSE-KEY>
  
  For details, see the Gurobi licensing guide: https://www.gurobi.com/documentation/

Test:
======
   
   ## Create GEM for _C. albicans_ #
   
   Structure alignment mode
   
      python AlphaGEM.py --name candida --mode structure_alignment --fasta ./working/candida/candida.fasta --refname yeast
      
   PLMsearch mode
   
      python AlphaGEM.py --name candida --mode plmsearch --fasta ./working/candida/candida.fasta --refname yeast
 
Usage:
======

 ## Parameters: #
  - **--mode**: Workflow mode. Default: `structure_alignment`. Choices: `structure_alignment`, `plmsearch`
  - **--refname**: Reference species. Choices: `ecoli`, `yeast`, `strco`, `human`, `synechocystis`. Required
  - **--name**: Job/species name used to create `working/<name>`. Required
  - **--fasta**: Target species FASTA file path. Required
  - **--list**: Path to a list file of structures and gene names. Default: empty
  - **--structure**: Directory containing structure files. Default: empty
  - **--cleanuse**: Whether CLEAN has been used. To enable, pass `--cleanuse True`; to disable, omit the flag
  - **--TMscore**: Structure alignment filter threshold. Default: 0.7
  - **--upTMscore**: Safe TMscore. Default: 0.9
  - **--TMscoretrans**: Transporter TMscore filter. Default: 0.7
  - **--coverage**: Coverage filter threshold. Default: 0.8
  - **--upcoverage**: Safe coverage. Default: 0.9
  - **--coveragetrans**: Transporter coverage filter. Default: 0.8
  - **--pLDDT**: pLDDT filter threshold. Default: 70
  - **--esp**: Clustering parameter. Default: 1
  - **--grothmedium**: Growth medium. Default: `min`; choices: `min`, `full`

 ## Tutorials: #

 ## Examples: #
  - **Structure alignment mode (default)**
  
    > **Pre-requisites for Structure Alignment:**
    > If you perform structure alignment, please provide:
    > 1. A directory containing structure data via the `--structure` parameter.
    > 2. A genome-wide table file via the `--list` parameter. 
    >
    > **List File Format Requirements:**
    > *   The table must contain the following columns: **Entry**, **Entry Name**, **Structure**.
    > *   The content of the **Entry** column must strictly match the gene IDs provided in your `--fasta` file.

        python AlphaGEM.py \
          --mode "structure alignment" \
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