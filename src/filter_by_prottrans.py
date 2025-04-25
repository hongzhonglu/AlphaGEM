import gc
import tqdm
import cobra
from transformers import T5Tokenizer, T5EncoderModel
import torch
import pickle
import re
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import pandas as pd
from Bio import SeqIO
import findrefnames

def embedding(sequences):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Load the tokenizer
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False, legacy=False)

    # Load the model
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(device)

    # only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
    if device == torch.device("cpu"):
        model.to(torch.float32)

    # prepare your protein sequences as a list
    sequence_examples = sequences

    # replace all rare/ambiguous amino acids by X and introduce white-space between all amino acids
    sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]

    # tokenize sequences and pad up to the longest sequence in the batch
    ids = tokenizer(sequence_examples, add_special_tokens=True, padding="longest")

    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)

    # generate embeddings
    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
    embs=[]
    for number,sequence in enumerate(sequences):
        embs.append(embedding_repr.last_hidden_state[number,:len(sequence)].mean(dim=0).numpy())
    return embs

def findtargethomos(data,dataindex):
    pca = PCA(n_components=2)  # 将数据降至 2 维，便于可视化
    data_reduced = pca.fit_transform(data)
    target_data=data[0]
    distances = pairwise_distances([target_data], data, metric='euclidean')[0]
    distance_threshold = 1.1*min(distances[1:])
    close_points = distances <= distance_threshold
    return dataindex[close_points]

def filter_by_prottrans(name,refname):
    betch_size=1
    homologs = pd.read_excel(f'working/{name}/matrix_homolog{name}_preforprottrans.xlsx')
    groupedhomologs={homos[0]:list(homos[1][1]) for homos in homologs.groupby(0)}
    if refname == 'yeast':
        refmodel = cobra.io.read_sbml_model('models/yeast-GEM.xml')
    if refname == 'ecoli':
        refmodel = cobra.io.load_json_model('models/iML1515.json')
    if refname == 'strco':
        refmodel = cobra.io.read_sbml_model('models/Sco-GEM.xml')
    if refname == 'human':
        refmodel = cobra.io.read_sbml_model('models/Human-GEM.xml')
    genes=[gene.id for gene in refmodel.genes]
    final_homologs=[]
    sortedgenes={}
    for i in groupedhomologs.keys():
        if i in genes:
            sortedgenes[i]=groupedhomologs[i]
        else:
            for gene in groupedhomologs[i]:
                final_homologs.append([i,gene])
    indexes=[]
    genes2embedding=[]
    for i in sortedgenes.keys():
        if len(list(set(sortedgenes[i])))>1:
            genes2embedding.append((i,list(set(sortedgenes[i]))))
        else:
            final_homologs.append([i,sortedgenes[i][0]])
    indexseq=[]
    records1={record.id.split('|')[1]:str(record.seq) for record in SeqIO.parse(f'./data_available/{refname}.fasta','fasta')}
    records2={record.id.split('|')[1]:str(record.seq) for record in SeqIO.parse(f'./working/{name}/{name}.fasta','fasta')}
    findrefnames.predata(refname)
    for genes in genes2embedding:
        indexes.append(findrefnames.findmodelname(genes[0]))
        indexes=indexes+genes[1]
    indexes=list(set(indexes))
    for i in indexes:
        try:
          indexseq.append(records1[i])
        except:
            try:
               indexseq.append(records2[i])
            except:
                indexseq.append('PROT')
    indexemb=[]
    for seqs in tqdm.tqdm(range(0,len(indexseq),betch_size)):
        a=embedding(indexseq[seqs:seqs+betch_size])
        indexemb+=a
        del a
        gc.collect()
    embeddings={}
    for i in range(len(indexes)):
        embeddings[indexes[i]]=indexemb[i]
    del indexemb,indexseq,records1,records2,indexes
    with open('working/{name}/{name}.emb'.format(name=name), 'wb') as handle:
        # noinspection PyTypeChecker
        pickle.dump(embeddings, handle)
    # with open('working/{name}/{name}.emb'.format(name=name), 'rb') as handle:
    #     embeddings=pickle.load(handle)
    for genes in genes2embedding:
        try:
           data=[embeddings[findrefnames.findmodelname(genes[0])]]
        except:
            data=embedding([records1[findrefnames.findmodelname(genes[0])]])
        dataindex=[genes[0]]
        for gene in genes[1]:
            if gene=='A2ASS6':
                continue
            try:
                data.append(embeddings[gene])
            except:
                try:
                    data.append(embedding([records2[gene]])[0])
                except:
                    continue
            dataindex.append(gene)
        dataindex=np.array(dataindex)
        homoss=findtargethomos(data,dataindex)[1:]
        print(homoss)
        for homomomo in homoss:
            final_homologs.append([genes[0],homomomo])
    pd.DataFrame(final_homologs).to_excel(f'working/{name}/matrix_homolog{name}.xlsx')


