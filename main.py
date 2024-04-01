from config import refname
from config import name
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import time
t1=time.time()
import orthofinder_datahandle
import US_align_knock1
import US_align_choose
import foldseek_knock1
import foldseek_choose
import foldseek_choose2
import homolog_concat2
import model_build1
import model_build2
import use_eggnog
import kegg_find
import model_reaction
import gapfilling
path_taryeast_structure=''
path_reference_structure=''
if refname=='ecoli':
    refmodel='iML1515.json'
if refname=='yeast':
    refmodel='yeast-GEM.xml'
if refname=='strco':
    refmodel='Sco-GEM.xml'
if refname=='human':
    refmodel='Human-GEM.xml'
def main():
    orthofinder_datahandle.datahandel(name,refname)
    US_align_knock1.US_align_find(name, path_taryeast_structure,refname,path_reference_structure)
    US_align_choose.US_align_choose()
    foldseek_knock1.foldseekfind(path_taryeast_structure, name, refname)
    foldseek_choose.foldseek_choose(name, refname, path_taryeast_structure)
    foldseek_choose2.foldseek_choose2(refname)
    homolog_concat2.homo()
    model_build1.modelbuild(refmodel,name)
    model_build2.modelbuild(refmodel, name)
    use_eggnog.eggnog(name,refname, 1)
    kegg_find.reaction_get()
    kegg_find.genes_get(name)
    model_reaction.modelreaction(name)
    gapfilling.gapfill(name, refname)

if __name__=='__main__':
    main()
t2=time.time()
print(t2-t1)
print('s')