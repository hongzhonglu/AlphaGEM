Installation:
======
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEMs.git AlphaGEM

You can download supplyment resources from PanBaidu:

  ## envrionment #

* 1.create environment for AlphaGEM

      conda create -n AlphaGEM python==3.10

* 2.You can install the following packages through conda for AlphaGEM:
  
      conda install diamond -c bioconda
  
* 3.the following packages through pip
  
      pip install biopython
      pip install cobra
      pip install numpy
      pip install pandas
  
## use CLEAN #

create enviroment for clean:
     

Usage:
======
## Basic Use #
     main.py --name species_name --refname reference_species_name --fasta fasta_file --cleanuse True
