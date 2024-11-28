import os
import pandas as pd
import shutil
import numpy as np
from Bio import SeqIO
from plmsearch import rhea_search
from ss_predict_choose import ss_predict_choose_rhea

def process_rhea(name,refname,threshold,left):
    rhea_homo=pd.read_excel(f"./working/{name}/{name}_rhea_homolog{str(threshold)}.xlsx",index_col=0)

    header=["ID","Query"]
    rhea_homo.columns=header
    rhea_homo = rhea_homo[~rhea_homo["Query"].isin(left)]
    gene_rhea=pd.read_csv(f"./rhea/rhea2uniprot_sprot.tsv",sep="\t")
    gene2rheaid=gene_rhea.set_index("ID")["RHEA_ID"].to_dict()
    rhea_homo["Rhea_id"]=rhea_homo["ID"].map(gene2rheaid)
    rhea_homo.to_excel(f"./working/{name}/{name}_rheaid_plmsearch.xlsx")

    return  rhea_homo
def rhea(name,refname,threshold=0.9):
    path=os.getcwd()
    tary = []
    homo = pd.read_excel(f'./working/{name}/matrix_homolog{name}.xlsx')
    homo2 = []
    left = []
    
    for record in SeqIO.parse(f"./working/{name}/{name}.fasta","fasta"):
        tary.append(record.id)

    for i in range(len(homo.index)):
        homo2.append(homo.iat[i, 2])

    for i in tary:
        if i not in homo2:
            left.append(i)
    left = pd.DataFrame(left)

    rhea_search(name,refname,threshold)
    ss_predict_choose_rhea(name,threshold)
    process_rhea(name,refname,threshold,left)

    # with open(f'./working/{name}/test{name}.emapper.annotations') as f:
    #     content = f.read()
    #     contents = content.split('\n')
    # genes = pd.DataFrame()
    # for i in contents:
    #     j = i.split('\t')
    #     j = np.transpose(pd.DataFrame((np.array(j))))
    #     genes = pd.concat([genes, j])
    # genes.index = range(len(genes.index))
    # for i in range(len(genes.index)):
    #     try:
    #         genes.iat[i, 0] = genes.iat[i, 0].split('|')[2]
    #     except:
    #         continue
    # genes1 = []
    # for i in range(len(genes.index)):
    #     genes1.append(genes.iat[i, 0])
    # left.insert(loc=1, column='EC', value='')
    # left.insert(loc=2, column='Rec', value='')
    # for i in range(len(left.index)):
    #     try:
    #         left.iat[i, 1] = genes.iat[genes1.index(left.iat[i, 0]), 10]
    #         left.iat[i, 2] = genes.iat[genes1.index(left.iat[i, 0]), 14]
    #     except:
    #         left.iat[i, 1], left.iat[i, 2] = '-', '-'
    # eggnogannop = pd.DataFrame()
    # for i in range(len(left.index)):
    #     if left.iat[i, 2] != '-':
    #         eggnogannop = pd.concat([eggnogannop, left.iloc[[i], [0, 1, 2]]])
    # eggnogannop.index = range(len(eggnogannop.index))
    # eggnogannop.to_excel(f'./working/{name}/eggec2{name}.xlsx')





# def rhea_search(name,refname,threshold=0.8,save_model_path="./plmsearch/model/plmsearch.sav"):
#     input_query_embedding =f"./working/{name}/{name}_embedding.pkl"
#     if torch.cuda.is_available() == False:
#         print("GPU selected but none of them is available.")
#         device = "cpu"
#     else:
#         print("We have", torch.cuda.device_count(), "GPUs in total!, we will use as you selected")
#         device = f'cuda:{0}'
#     with open(input_query_embedding, 'rb') as handle:
#         query_embedding_dic = pickle.load(handle)
#
#     model = plmsearch(embed_dim=1280)
#     model.load_pretrained(save_model_path)
#     model.eval()
#     model_methods = model
#     if (device != "cpu"):
#         model = nn.DataParallel(model, device_ids=[0])
#         model_methods = model.module
#     model.to(device)
#
#     output_search_result = f"./working/{name}/{name}_ss_predictor_filter{str(threshold)}"
#
#     with open(output_search_result, 'w') as f:
#         pass
#
#     batch_size = 50
#     query_keys_all = list(query_embedding_dic.keys())
#
#     search_dict = None
#
#     for i in range(1,7):
#         input_target_embedding=(f"./rhea/reah_prot{i}.pkl")
#         with open(input_target_embedding, 'rb') as handle:
#             target_embedding_dic = pickle.load(handle)
#
#
#         for i in range(0, len(query_keys_all), batch_size, desc="Search query proteins batch by batch"):
#             batch_keys = query_keys_all[i:i + batch_size]
#             batch_query_embedding_dic = {k: query_embedding_dic[k] for k in batch_keys}
#             # get sorted pairlist
#             if (save_model_path != None):
#                 search_result = plmsearch_search(batch_query_embedding_dic, target_embedding_dic, device, model_methods,
#                                                  search_dict)
#             else:
#                 search_result = esm_similarity_search(batch_query_embedding_dic, target_embedding_dic, device)
#             with open(output_search_result, 'a') as f:
#                 for protein in search_result:
#                     for pair in search_result[protein]:
#                         if float(pair[1]) > threshold:
#                             f.write(f"{protein}\t{pair[0]}\t{pair[1]}\n")
#     return None
#



# def ss_predictor(refname,name,threshold=0.8,save_model_path="./plmsearch/model/plmsearch.sav"):
#     input_query_embedding = f"./working/{name}/{name}_embedding.pkl"
#     input_target_embedding=f"./working/{refname}_embedding.pkl"
#     with open(input_query_embedding, 'rb') as handle:
#         query_embedding_dic = pickle.load(handle)
#     with open(input_target_embedding, 'rb') as handle:
#         target_embedding_dic = pickle.load(handle)
#     if torch.cuda.is_available() == False:
#         print("GPU selected but none of them is available.")
#         device = "cpu"
#     else:
#         print("We have", torch.cuda.device_count(), "GPUs in total!, we will use as you selected")
#         device = f'cuda:{0}'
#
#     model = plmsearch(embed_dim=1280)
#     model.load_pretrained(save_model_path)
#     model.eval()
#     model_methods = model
#     if (device != "cpu"):
#         model = nn.DataParallel(model, device_ids=[0])
#         model_methods = model.module
#     model.to(device)
#
#     output_search_result = f"./working/{name}/{name}_ss_predictor_filter{str(threshold)}"
#
#     with open(output_search_result, 'w') as f:
#         pass
#
#     batch_size = 50
#     query_keys_all = list(query_embedding_dic.keys())
#
#     search_dict = None
#
#     for i in trange(0, len(query_keys_all), batch_size, desc="Search query proteins batch by batch"):
#         batch_keys = query_keys_all[i:i + batch_size]
#         batch_query_embedding_dic = {k: query_embedding_dic[k] for k in batch_keys}
#         # get sorted pairlist
#         if (save_model_path != None):
#             search_result = plmsearch_search(batch_query_embedding_dic, target_embedding_dic, device, model_methods,
#                                              search_dict)
#         else:
#             search_result = esm_similarity_search(batch_query_embedding_dic, target_embedding_dic, device)
#         with open(output_search_result, 'a') as f:
#             for protein in search_result:
#                 for pair in search_result[protein]:
#                     if float(pair[1]) > threshold:
#                         f.write(f"{protein}\t{pair[0]}\t{pair[1]}\n")
