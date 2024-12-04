import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import gc
def cluster(spe,clusteresp=0.6,uptm=0.90,upcov=0.90):
   data = pd.read_excel(f'working/{spe}/matrix_foldseek_filtered2_{spe}.xlsx')[['gene1','gene2','tmscore', 'coverage', 'sumscore']]
   data.reset_index(drop=True, inplace=True)
   # 对于每个name1进行聚类
   result = []
   for name, group in data.groupby('gene1'):
       if len(group.index)==1:
           group['cluster']=-1
           for index,row in group.iterrows():
              result.append(row)
           continue
       # 提取数值列并标准化，包含数值3
       values = group[['tmscore', 'coverage', 'sumscore']].values
       scaler = StandardScaler()
       values_scaled = scaler.fit_transform(values)
       # 使用DBSCAN进行聚类
       dbscan = DBSCAN(eps=clusteresp, min_samples=1)
       group['cluster'] = dbscan.fit_predict(values_scaled)
       # 找到非噪声类（-1表示噪声）
       non_noise_clusters = group[group['cluster'] != -1]
       noise_clusters = group[group['cluster'] == -1]
       if not non_noise_clusters.empty:
           non_noise_clusters['cluster'] = non_noise_clusters['cluster'].astype(int)
           highest_cluster = non_noise_clusters.groupby('cluster')['sumscore'].mean().idxmax()
           highest_elements = non_noise_clusters[non_noise_clusters['cluster'] == highest_cluster]
           for index,row in noise_clusters.iterrows():
               if row['sumscore'] > highest_elements['sumscore'].mean():
                   result.append(row)
                   continue
           for index,row in highest_elements.iterrows():
               result.append(row)
       for index,row in group.iterrows():
           if row['tmscore']>=uptm and row['coverage']>=upcov:
              result.append(row)
       del values, values_scaled, group, non_noise_clusters, highest_elements
       gc.collect()
       # #使用KMeans聚类
       # if len(group.index) >=3:
       #     kmeans = KMeans(n_clusters=2)  # 可以根据需要调整聚类数量
       # if len(group.index) >=7:
       #     kmeans = KMeans(n_clusters=3)
       # group['cluster'] = kmeans.fit_predict(values)
       # # 找到数值最高的聚类
       # highest_cluster = group.groupby('cluster').sum().idxmax()
       # # 提取最高类的元素
       # highest_elements = group[group['cluster'] == highest_cluster]
       # result.append(highest_elements)
   final_result = pd.concat(result,axis=1).T.reset_index(drop=True)[['gene1','gene2']]
   final_result.to_excel(f'working/{spe}/matrix_homo_part2{spe}.xlsx')
