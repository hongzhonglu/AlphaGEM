# from config import refname
# from config import name
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

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
import dl_module
# import ss_predict_choose
# import cov_filter
import os
import pandas as pd
import tqdm
# path_taryeast_structure=''
# path_reference_structure=''
import use_rhea
import mkdir


threshold=0.8
cov_thre=0.25
cov_threshold=0.5
pid_thre=60
id_threshold=20
direction="direct"
refname="yeast"
if refname=='ecoli':
    refmodel='iML1515.json'
if refname=='yeast':
    refmodel='yeast-GEM.xml'
if refname=='strco':
    refmodel='Sco-GEM.xml'
if refname=='human':
    refmodel='Human-GEM.xml'


# species_df=pd.read_excel("./1000_construction/mapping_df4.xlsx")
import pickle
with open("./tmpe/fungi_list3.pkl","rb") as f:
    strain_list=pickle.load(f)
# # for index,row in tqdm.tqdm(species_df.iterrows()):
# #     if index>300:
# #         name=row["index"]
# #         #if not os.path.exists(f"./working/{name}/")
# #         dl_module.dl_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
# c=0
# for index in tqdm.tqdm(["mouse"]):
#     name=index
#     # if not os.path.exists(f"./working/{name}/"):
#     dl_module.dl_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)

for index in tqdm.tqdm(strain_list[168:]):
    name=index
    # if not os.path.exists(f"./working/{name}/"):
    dl_module.dl_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)

for i in ["mouse"]:
    name=i
    mkdir.mvfile(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
    use_rhea.rhea(name,refname,threshold=0.9)
