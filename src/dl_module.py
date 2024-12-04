# from config import refname
# from config import name
# import warnings
# warnings.filterwarnings("ignore", category=FutureWarning)
import time
t1=time.time()
# import orthofinder_datahandle
# import US_align_knock1
# import US_align_choose
# import foldseek_knock1
# import foldseek_choose
# import foldseek_choose2
# import homolog_concat2
# import model_build1
# import model_build2
# import use_eggnog
# import kegg_find
# import model_reaction
# import gapfilling
import plmsearch
import plmsearch_embedding
import mkdir
# import ss_predict_choose
# import cov_filter
import os
# path_taryeast_structure=''
# path_reference_structure=''



def dl_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction):
    filename=f"./tmpe/{name}.fasta"
    filelist=os.listdir(f"./working/")
    if name not in filelist:
        mkdir.makeworkdir(name,filename)
    file_detail=os.listdir(f"./working/{name}/")
    if f"{name}_embedding.pkl" not in file_detail:
        plmsearch_embedding.embedding_generate(name)
    if f"{name}_ss_predictor_filter{str(threshold)}" not in file_detail:
        plmsearch.ss_predictor(refname, name, threshold)


    # ss_predict_choose.ss_predict_choose(name,refname,threshold)
    # cov_filter.result_blast(name,refname,threshold)
    # cov_filter.cov_filter(name,refname,threshold,cov_threshold,id_threshold)
    # cov_filter.bbh(name,refname,threshold,cov_thre,pid_thre,direction)
    # cov_filter.merge(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
    # mkdir.mvfile(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
    #cov_filter.evaluate(name,refname,threshold,cov_threshold,cov_thre,pid_thre,direction)
    # orthofinder_datahandle.datahandel(name,refname)
    # US_align_knock1.US_align_find(name, path_taryeast_structure,refname,path_reference_structure)
    # US_align_choose.US_align_choose()
    # foldseek_knock1.foldseekfind(path_taryeast_structure, name, refname)
    # foldseek_choose.foldseek_choose(name, refname, path_taryeast_structure)
    # foldseek_choose2.foldseek_choose2(refname)
    # homolog_concat2.homo(name)
    # model_build1.modelbuild(refmodel,name)
    # model_build2.modelbuild(refmodel, name)
    # use_eggnog.eggnog(name,refname, 1)
    # kegg_find.reaction_get()
    # kegg_find.genes_get(name)
    # model_reaction.modelreaction(name)
    # gapfilling.gapfill(name, refname)

