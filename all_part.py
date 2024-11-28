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
#import addreaction_module
# import ss_predict_choose
# import cov_filter
import os
import pandas as pd
import tqdm
import filter_module
import premodelbuild_module
import addreaction_module
# path_taryeast_structure=''
# path_reference_structure=''

threshold=0.8
cov_thre=0.25
cov_threshold=0.5
pid_thre=60
id_threshold=20
direction="direct"
refname="human"
if refname=='ecoli':
    refmodel='iML1515.json'
if refname=='yeast':
    refmodel='yeast-GEM.xml'
if refname=='strco':
    refmodel='Sco-GEM.xml'
if refname=='human':
    refmodel='Human-GEM.xml'


species_df=pd.read_excel("./1000_construction/mapping_df4.xlsx")
already=list(species_df["index"])
import pickle
with open("./1000_construction/y1000_indexname.pkl","rb") as f:
    strain_list=pickle.load(f)
for index in tqdm.tqdm(["mouse"]):
    name=index
    if name not in already:
        filelist=os.listdir(f"./working/{name}/")
        if any(item.endswith(f"filter{threshold}") for item in filelist):
            if not any(item.endswith("direct.xlsx") for item in filelist):
                filter_module.filter_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
            # if not any(item.endswith(f"tarmodel{name}.yml") for item in filelist):
            #     premodelbuild_module.prebuild_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
            # if not any(item.endswith(f"text.xml") for item in filelist):
            #     addreaction_module.prebuild_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)

# for i in ["candida"]:
#     name=i
#     filelist = os.listdir(f"./working/{name}/")
#     if any(item.endswith(f"filter{threshold}") for item in filelist):
#         if not any(item.endswith("direct.xlsx") for item in filelist):
#             filter_module.filter_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)
#         if not any(item.endswith(f"tarmodel{name}.yml") for item in filelist):
#             premodelbuild_module.prebuild_module(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction)