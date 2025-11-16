AlphaGEM
======
 ## Platform: Ubuntu only (tested on Ubuntu 20.04/22.04). #

 
Installation:
======
## Install codes and datas: #
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEMs.git AlphaGEM

You can download supplyment resources from PanBaidu:

## Prepare envrionment: #
It need to create two enviroment as DeepECtransformer uses the old version of pytorch.

* 1.create environment for AlphaGEM:

      conda env create -f AlphaGEM.yml

* 2.create enviroment for deepectransformer:
  
      conda env create -f deepecenvironment.yml
  
* 3.activate AlphaGEM environment:
  
      conda activate AlphaGEM
  
* 4.prepare the source data: 
  
      python setup.py install
     

Usage:
======

 ## Parameters: #
  - **--mode**: Workflow mode. Default: `structure alignment`. Choices: `structure alignment`, `plmsearch`
  - **--refname**: Reference species. Choices: `ecoli`, `yeast`, `strco`, `human`. Required
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

 ## Examples: #
  - **Structure alignment mode (default)**

        python AlphaGEM.py \
          --mode "structure alignment" \
          --name my_species \
          --refname yeast \
          --fasta ./input.fa \
          --cleanuse True \
          --TMscore 0.7 --coverage 0.8 --pLDDT 70

  - **PLMSearch mode**

        python AlphaGEM.py \
          --mode plmsearch \
          --name my_species \
          --refname ecoli \
          --fasta ./input.fa
