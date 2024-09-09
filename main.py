import fasta_handle
import argparse
import warnings
import time
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
import nonhomogene_anno_score_from_eggandclean
import add_nonhomo_reaction_score
import gapfilling
from config import refmodel

warnings.filterwarnings("ignore", category=FutureWarning)
t1=time.time()
path_taryeast_structure=''
path_reference_structure=''


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--refname', type=str, help='reference name')
    parser.add_argument('--name', type=str, help='target GEMs name')
    parser.add_argument('--fasta', type=str, help='your target species genome')
    parser.add_argument('--cleanuse', type=bool, help='whether you have used CLEAN')
    args = parser.parse_args()
    fasta=args.fasta
    refname=args.refname
    refmodel=''
    name=args.name
    cleanuse=args.cleanuse
    if refname == 'ecoli':
        refmodel = 'iML1515.json'
    if refname == 'yeast':
        refmodel = 'yeast-GEM.xml'
    if refname == 'strco':
        refmodel = 'Sco-GEM.xml'
    if refname == 'human':
        refmodel = 'Human-GEM.xml'
    fasta_handle.handle(fasta,name)
    orthofinder_datahandle.datahandel(name,refname)
    US_align_knock1.US_align_find(name, path_taryeast_structure,refname,path_reference_structure)
    US_align_choose.US_align_choose(name)
    foldseek_knock1.foldseekfind(path_taryeast_structure, name, refname)
    foldseek_choose.foldseek_choose(name, refname, path_taryeast_structure)
    foldseek_choose2.foldseek_choose2(refname,refmodel,name)
    homolog_concat2.homo(name)
    model_build1.modelbuild(refmodel,name)
    model_build2.modelbuild(refmodel, name)
    use_eggnog.eggnog(name,refname, 1)
    #kegg_find.reaction_get(name)
    #kegg_find.genes_get(name)
    #model_reaction.modelreaction(name)
    nonhomogene_anno_score_from_eggandclean.nonhome(name,cleanuse)
    add_nonhomo_reaction_score.model_reaction(name,refname)

    gapfilling.gapfill(name, refname)

if __name__=='__main__':
    main()
t2=time.time()
print(t2-t1)
print('s')
