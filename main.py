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
import use_clean
import use_eggnog
import multiple_annotition_for_nonhomogene
import add_nonhomo_reaction_score
import gapfilling
import use_deepectransformer

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
    parser.add_argument('--TMscore', type=float, help='filter TMscore')
    parser.add_argument('--TMscoretrans', type=float, help='filter TMscoretrans')
    parser.add_argument('--coverage', type=float, help='filter coverage')
    parser.add_argument('--coveragetrans',type=float, help='filter coveragetrans')
    parser.add_argument('--pLDDT', type=float, help='filter pLDDT')
    args = parser.parse_args()

    fasta=args.fasta
    refname=args.refname
    refmodel=''
    name=args.name
    TMscore=0.8
    TMscoretrans=0.8
    coverage=0.9
    coveragetrans=0.8
    pLDDT=70
    TMscore = args.TMscore
    TMscoretrans=args.TMscoretrans
    coverage = args.coverage
    cleanuse=args.cleanuse
    coveragetrans=args.coveragetrans
    pLDDT=args.pLDDT
    if refname == 'ecoli':
        refmodel = 'iML1515.json'
    if refname == 'yeast':
        refmodel = 'yeast-GEM.xml'
    if refname == 'strco':
        refmodel = 'Sco-GEM.xml'
    if refname == 'human':
        refmodel = 'Human-GEM.xml'
    #fasta_handle.handle(fasta,name)
    #orthofinder_datahandle.datahandel(name,refname)
    #US_align_knock1.US_align_find(name, path_taryeast_structure,refname,path_reference_structure)
    US_align_choose.US_align_choose(name)
    #foldseek_knock1.foldseekfind(path_taryeast_structure, name, refname)
    #foldseek_choose.foldseek_choose(name, refname, path_taryeast_structure)
    foldseek_choose2.foldseek_choose2(refname,refmodel,name,TMscoretrans,coveragetrans,TMscore,coverage,pLDDT)
    homolog_concat2.homo(name)
    model_build1.modelbuild(refmodel,name)
    model_build2.modelbuild(refmodel, name)
    use_eggnog.eggnog(name,refname, 1)
    use_clean.clean_result(name)
    use_deepectransformer.use_deepectransformer(name)
    #kegg_find.reaction_get(name)
    #kegg_find.genes_get(name)
    #model_reaction.modelreaction(name)
    multiple_annotition_for_nonhomogene.nonhome(name,cleanuse,True)
    add_nonhomo_reaction_score.model_reaction(name,refname)

    gapfilling.gapfill(name, refname)

if __name__=='__main__':
    main()
t2=time.time()
print(t2-t1)
print('s')
