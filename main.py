import fasta_handle
import argparse
import warnings
import time
import plmsearch_embedding
import plmsearch
import ss_predict_choose
import mkdir
import cov_filter
import orthofinder_datahandle
import US_align_knock1
import US_align_choose
import foldseek_knock1
import foldseek_choose
import foldseek_choose2
import foldseekcluster
import homolog_concat2
import model_build1
import model_build2
import use_clean
import use_eggnog
import use_rhea
import multiple_annotition_for_nonhomogene
import add_nonhomo_reaction_score
import gapfilling
import use_deepectransformer
import add_reactions_based_pool
import os
import generate_list
from C_sourcegapseqmodel import medium

warnings.filterwarnings("ignore", category=FutureWarning)
t1=time.time()
path_taryeast_structure=''
path_reference_structure=''


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--mode',type=str,default='structure alignment',choices=['structure alignment','plmsearch'])
    parser.add_argument('--refname', type=str, help='reference name')
    parser.add_argument('--name', type=str, help='target GEMs name')
    parser.add_argument('--fasta', type=str, help='your target species genome')
    parser.add_argument('--list', type=str,default='', help='lists of structures and genenames')
    parser.add_argument('--structure', type=str,default='', help='files you store structures')
    parser.add_argument('--cleanuse', type=bool, help='whether you have used CLEAN')
    parser.add_argument('--TMscore', type=float,default=0.7,help='filter TMscore')
    parser.add_argument('--upTMscore', type=float, default=0.9, help='safe TMscore')
    parser.add_argument('--TMscoretrans', type=float,default=0.7, help='filter TMscoretrans')
    parser.add_argument('--coverage', type=float,default=0.8, help='filter coverage')
    parser.add_argument('--upcoverage', type=float, default=0.9, help='safe coverage')
    parser.add_argument('--coveragetrans',type=float,default=0.8, help='filter coveragetrans')
    parser.add_argument('--pLDDT', type=float,default=70, help='filter pLDDT')
    parser.add_argument('--spe', type=float, default=1, help='cluster spe')
    parser.add_argument('--grothmedium',type=str,default='min', help='grothmedium',choices=['min','full'])
    args = parser.parse_args()
    mode=args.mode
    fasta=args.fasta
    refname=args.refname
    refmodel=''
    name=args.name
    listfile=args.list
    structurefile=args.structure
    TMscore = args.TMscore
    TMscoretrans=args.TMscoretrans
    upTMscore=args.upTMscore
    coverage = args.coverage
    cleanuse=args.cleanuse
    coveragetrans=args.coveragetrans
    upcoverage=args.upcoverage
    pLDDT=args.pLDDT
    spe=args.spe
    grothmedium=args.grothmedium
    if refname == 'ecoli':
        refmodel = 'iML1515.json'
    if refname == 'yeast':
        refmodel = 'yeast-GEM.xml'
    if refname == 'strco':
        refmodel = 'Sco-GEM.xml'
    if refname == 'human':
        refmodel = 'Human-GEM.xml'
    try:
        os.mkdir(f'./working/{name}')
    except FileExistsError:
        pass
    fasta_handle.handle(fasta,name)
    plmsearch_embedding.embedding_generate(name)
    if mode=='structure alignment':
        generate_list.generare_list_with_structure(name,listfile,structurefile)
        orthofinder_datahandle.datahandel(name,refname)
        US_align_knock1.US_align_find(name, path_taryeast_structure,refname,path_reference_structure)
        US_align_choose.US_align_choose(name)
        foldseek_knock1.foldseekfind(path_taryeast_structure, name, refname)
        foldseek_choose.foldseek_choose(name, refname, path_taryeast_structure)
        foldseek_choose2.foldseek_choose2(refname,refmodel,name,TMscoretrans,coveragetrans,TMscore,coverage,pLDDT)
        foldseekcluster.cluster(name,spe,uptm=upTMscore,upcov=upcoverage)
        homolog_concat2.homo(name)
    if mode=='plmsearch':
        threshold = 0.8
        cov_thre = 0.25
        cov_threshold = 0.5
        pid_thre = 60
        id_threshold = 20
        direction = "direct"
        generate_list.generare_list_without_structure(name)
        plmsearch.ss_predictor(refname, name)
        ss_predict_choose.ss_predict_choose(name, refname, threshold)
        cov_filter.result_blast(name, refname, threshold)
        cov_filter.cov_filter(name, refname, threshold, cov_threshold, id_threshold)
        cov_filter.bbh(name, refname, threshold, cov_thre, pid_thre, direction)
        cov_filter.merge(name, refname, threshold, cov_threshold, id_threshold, cov_thre, pid_thre, direction)
        mkdir.mvfile(name, refname, threshold, cov_threshold, id_threshold, cov_thre, pid_thre, direction)
    model_build1.modelbuild(refmodel,name)
    model_build2.modelbuild(refmodel, name)
    use_eggnog.eggnog(name,refname, 1)
    use_clean.clean_result(name)
    use_deepectransformer.use_deepectransformer(name)
    use_rhea.rhea(name,refname)
    #kegg_find.reaction_get(name)
    #kegg_find.genes_get(name)
    #model_reaction.modelreaction(name)
    multiple_annotition_for_nonhomogene.nonhome(name,cleanuse,True)
    add_reactions_based_pool.model_reaction(name,refname)
    gapfilling.gapfill(name, refname,grothmedium)


if __name__=='__main__':
    main()
t2=time.time()
print(t2-t1)
print('s')
