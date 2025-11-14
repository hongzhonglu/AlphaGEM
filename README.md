Installation:
======
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEMs.git AlphaGEM

You can download supplyment resources from PanBaidu:

  ## envrionment #
you should create two enviroment as Deepectransformer use the old version of transformer.

* 1.create environment for AlphaGEM:

      conda create -n AlphaGEM 

* 2.create enviroment for deepectransformer:
  
      conda create -n deepectransformer
  
* 3.prepare the source data: 
  
      bash ./setup.sh

* 4.prepare dataset for eggnog-mapper:

      python ./tools/eggnog_mapper/download_eggnog_database.py     

Usage:
======
## Basic Use #
     python AlphaGEM.py --name species_name --refname reference_species_name --fasta fasta_file --cleanuse True

  Parameters:
  ========
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

  Examples:
  ========
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
