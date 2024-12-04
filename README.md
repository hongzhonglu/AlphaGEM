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

     

Usage:
======
## Basic Use #
     main.py --name species_name --refname reference_species_name --fasta fasta_file --cleanuse True
