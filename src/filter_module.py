# from config import refname
# from config import name
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
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
# import plmsearch
# import plmsearch_embedding
import mkdir
import ss_predict_choose
import cov_filter
import os
# path_taryeast_structure=''
# path_reference_structure=''

# threshold=0.8
# cov_thre=0.25
# cov_threshold=0.5
# pid_thre=60
# id_threshold=20
# direction="direct"
# if refname=='ecoli':
#     refmodel='iML1515.json'
# if refname=='yeast':
#     refmodel='yeast-GEM.xml'
# if refname=='strco':
#     refmodel='Sco-GEM.xml'
# if refname=='human':
#     refmodel='Human-GEM.xml'
#
# refname="yeast"



def filter_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction):
    # filelist=os.lisdir(f"./working/"):
    # if name not in filelist:
    #     mkdir.makeworkdir(name)
    # plmsearch_embedding.embedding_generate(name)
    # plmsearch.ss_predictor(name,refname,threshold)

    # filelist=os.lisdir(f"./working/"):
    # if name not in filelist:
    #     mkdir.makeworkdir(name)
    ss_predict_choose.ss_predict_choose(name,refname,threshold)
    cov_filter.result_blast(name,refname,threshold)
    cov_filter.cov_filter(name,refname,threshold,cov_threshold,id_threshold)
    cov_filter.bbh(name,refname,threshold,cov_thre,pid_thre,direction)
    cov_filter.merge(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
    mkdir.mvfile(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
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

if __name__=='__main__':
    main()
t2=time.time()
print(t2-t1)
print('s')