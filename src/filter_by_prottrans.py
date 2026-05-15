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

def load_model():
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(f"Loading model on {device}...")
    
    # Load the tokenizer
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False, legacy=False)

    # Load the model
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(device)

    # Use half precision on GPU to save memory (significant reduction)
    if device.type == 'cuda':
        model = model.half()

    return model, tokenizer, device

def embedding(sequences, model, tokenizer, device):
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
    # Move to CPU once to avoid multiple device transfers if not needed immediately on GPU for other things
    # .float() might be needed if half precision was used
    last_hidden_states = embedding_repr.last_hidden_state.detach().cpu()
    
    for number, sequence in enumerate(sequences):
        # average over the sequence length (excluding padding if any, though attention mask handles computation, we need to correct length here)
        # Note: The original code used :len(sequence) which is correct for unpadded original length.
        # But wait, sequence_examples has spaces added. The tokenizer handles that.
        # The original code logic: embedding_repr.last_hidden_state[number,:len(sequence)]
        # However, len(sequence) is the length of the original amino acid string (e.g. "MVL...").
        # The tokenized input has spaces, so it's approx double length? No, tokenizer splits by space.
        # Let's double check T5Tokenizer. It treats "M V L" as 3 tokens?
        # ProtT5 tokenizer is character based mostly if space separated?
        # Actually better to stick to original logic: :len(sequence)
        # But wait, original code: 
        # sequence_examples = [" ".join(list...
        # ids = tokenizer(sequence_examples...
        # If I have "M V" (len 3 chars), tokenizer might produce tokens corresponding to M and V.
        # If strict preservation of original logic is desired:
        embs.append(last_hidden_states[number, :len(sequence)].mean(dim=0).numpy())
            
    return embs

def findtargethomos(data,dataindex,distance=1.2):
    pca = PCA(n_components=2)  # 将数据降至 2 维，便于可视化
    data_reduced = pca.fit_transform(data)
    target_data=data[0]
    distances = pairwise_distances([target_data], data, metric='euclidean')[0]
    distance_threshold = distance*min(distances[1:])
    close_points = distances <= distance_threshold
    return dataindex[close_points]

def findtargethomos_by_level(data, dataindex, level_map, low_distance=1.2, high_distance=None):
    if high_distance is None:
        high_distance = low_distance
    if len(data) <= 1:
        return dataindex
    target_data = data[0]
    distances = pairwise_distances([target_data], data, metric='euclidean')[0]
    min_distance = min(distances[1:])
    selected = [dataindex[0]]
    for idx in range(1, len(dataindex)):
        gene = dataindex[idx]
        level = level_map.get(gene, 'low')
        threshold_factor = high_distance if level == 'high' else low_distance
        if distances[idx] <= threshold_factor * min_distance:
            selected.append(gene)
    return np.array(selected)

def filter_by_prottrans(name,refname,distance,high_distance=1.8):
    batch_size=8
    model, tokenizer, device = load_model()
    
    # 追踪是否已切换到 CPU 模式（发生 OOM 后永久使用 CPU）
    use_cpu_mode = False
    current_device = device
    current_batch_size = batch_size
    
    homologs = pd.read_excel(f'working/{name}/matrix_homolog{name}_preforprottrans.xlsx')
    if 'confidence_level' not in homologs.columns:
        homologs['confidence_level'] = 'high'
    groupedhomologs = {}
    for refgene, group in homologs.groupby('refmodelgene'):
        groupedhomologs[refgene] = list(zip(group['tarmodelgene'], group['confidence_level']))
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
            for gene, _level in groupedhomologs[i]:
                final_homologs.append([i,gene])
    indexes=[]
    genes2embedding=[]
    for i in sortedgenes.keys():
        level_map = {}
        unique_genes = []
        for gene, level in sortedgenes[i]:
            if gene not in level_map:
                unique_genes.append(gene)
                level_map[gene] = level
            elif level == 'high':
                level_map[gene] = 'high'
        if len(unique_genes)>1:
            genes2embedding.append((i, unique_genes, level_map))
        else:
            final_homologs.append([i, unique_genes[0]])
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
    # ========== 两阶段处理策略 ==========
    # 第一阶段：GPU 批量处理，OOM 时逐条尝试，实在不行的记录下来
    # 第二阶段：CPU 逐条处理失败的序列
    
    indexemb = [None] * len(indexseq)  # 预分配，保持顺序
    failed_indices = []  # 记录 GPU 处理失败的序列索引
    
    # ===== 第一阶段：GPU 处理 =====
    if device.type == 'cuda':
        print(f"Phase 1: GPU batch processing (batch_size={batch_size})...")
        for batch_start in tqdm.tqdm(range(0, len(indexseq), batch_size), desc="GPU Processing"):
            batch_end = min(batch_start + batch_size, len(indexseq))
            batch_sequences = indexseq[batch_start:batch_end]
            batch_indices = list(range(batch_start, batch_end))
            
            try:
                # 尝试批量处理
                embeddings_batch = embedding(batch_sequences, model, tokenizer, device)
                for idx, emb in zip(batch_indices, embeddings_batch):
                    indexemb[idx] = emb
                    
            except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
                is_oom = isinstance(e, torch.cuda.OutOfMemoryError) or "out of memory" in str(e).lower()
                if is_oom:
                    print(f"\n GPU OOM at batch {batch_start}-{batch_end}, trying one by one...")
                    torch.cuda.empty_cache()
                    
                    # 逐条尝试 GPU 处理
                    for idx in batch_indices:
                        seq = indexseq[idx]
                        try:
                            emb = embedding([seq], model, tokenizer, device)[0]
                            indexemb[idx] = emb
                        except (torch.cuda.OutOfMemoryError, RuntimeError) as e2:
                            is_oom2 = isinstance(e2, torch.cuda.OutOfMemoryError) or "out of memory" in str(e2).lower()
                            if is_oom2:
                                # 单条也 OOM，留给 CPU
                                torch.cuda.empty_cache()
                                failed_indices.append(idx)
                            else:
                                raise e2
                else:
                    raise e
    else:
        # 没有 GPU，所有序列都放到 CPU 处理队列
        failed_indices = list(range(len(indexseq)))
    
    # ===== 第二阶段：CPU 逐条处理失败的序列 =====
    if failed_indices:
        print(f"\n Phase 2: CPU processing {len(failed_indices)} failed sequences (batch_size=1)...")
        
        # 切换模型到 CPU
        torch.cuda.empty_cache()
        gc.collect()
        cpu_device = torch.device('cpu')
        model.to(cpu_device).float()
        
        for idx in tqdm.tqdm(failed_indices, desc="CPU Processing"):
            seq = indexseq[idx]
            try:
                emb = embedding([seq], model, tokenizer, cpu_device)[0]
                indexemb[idx] = emb
            except Exception as e:
                print(f"\n Failed to embed sequence at index {idx}: {e}")
                indexemb[idx] = np.zeros(1024)  # ProtT5-XL 的 embedding 维度
        
        # 处理完后移回 GPU（用于后续操作）
        if device.type == 'cuda':
            model.to(device)
            model.half()
    
    # 当前使用的 device（用于后续 embedding 调用）
    current_device = device
    
    embeddings={}
    for i in range(len(indexes)):
        embeddings[indexes[i]]=indexemb[i]
    del indexemb,indexseq,indexes
    with open('working/{name}/{name}.emb'.format(name=name), 'wb') as handle:
        # noinspection PyTypeChecker
        pickle.dump(embeddings, handle)
    # with open('working/{name}/{name}.emb'.format(name=name), 'rb') as handle:
    #     embeddings=pickle.load(handle)
    for genes in genes2embedding:
        try:
           data=[embeddings[findrefnames.findmodelname(genes[0])]]
        except:
             # Pass model args to embedding (使用 current_device 以防已切换到 CPU)
             data=embedding([records1[findrefnames.findmodelname(genes[0])]], model, tokenizer, current_device)
        dataindex=[genes[0]]
        for gene in genes[1]:
            if gene=='A2ASS6':
                continue
            try:
                data.append(embeddings[gene])
            except:
                try:
                    data.append(embedding([records2[gene]], model, tokenizer, current_device)[0])
                except:
                    continue
            dataindex.append(gene)
        dataindex=np.array(dataindex)
        homoss=findtargethomos_by_level(
            data,
            dataindex,
            genes[2],
            low_distance=distance,
            high_distance=high_distance
        )[1:]
        for homomomo in homoss:
            final_homologs.append([genes[0],homomomo])
    del records1, records2
    pd.DataFrame(final_homologs).to_excel(f'working/{name}/matrix_homolog{name}.xlsx')
    
    # Clean up memory
    if 'model' in locals():
        del model
    if 'tokenizer' in locals():
        del tokenizer
    gc.collect()
    torch.cuda.empty_cache()

