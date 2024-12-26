import pandas as pd
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import tqdm
import matplotlib.pyplot as plt
import os
import glob
import tqdm
#from matplotlib_venn import venn2
def result_blast(name,refname,threshold):
    # 读取DataFrame
    yea = pd.read_excel(f"./data_available/{refname}.xlsx")
    gene_dict = yea.set_index('Gene Names (ordered locus)')['Entry'].to_dict()
    gene_dict2 = yea.set_index('Entry')['Gene Names (ordered locus)'].to_dict()

    df = pd.read_excel(f'./working/{name}/{name}_ss_homolog{threshold}.xlsx')  # 假设您的CSV文件名为gene_pairs.csv

    # 读取FASTA文件并创建基因名称到序列的映射
    def load_fasta(file):
        seq_dict = {}
        for record in SeqIO.parse(file, "fasta"):
            a=record.id
            if "|" in a:
                a=a.split("|")[1]
            seq_dict[a] = record.seq
        return seq_dict

    fasta1 = f"./data_available/{refname}.fasta"
    fasta2 = f"./working/{name}/{name}.fasta"

    seq_dict1 = load_fasta(fasta1)
    seq_dict2 = load_fasta(fasta2)

    # 进行BLAST比对并保存结果
    blast_results = []

    for index, row in df.iterrows():
        gene1 = gene_dict[row[0]]
        gene2 = row[1]

        if gene1 in seq_dict1 and gene2 in seq_dict2:
            seq1 = seq_dict1[gene1]
            seq2 = seq_dict2[gene2]

            # 创建临时FASTA文件
            with open(f"./working/{name}/temp1.fasta", "w") as f:
                f.write(f">{gene1}\n{seq1}\n")
            with open(f"./working/{name}/temp2.fasta", "w") as f:
                f.write(f">{gene2}\n{seq2}\n")

            # 运行BLAST比对
            blastp_cline = NcbiblastpCommandline(query=f"./working/{name}/temp1.fasta", subject=f"./working/{name}/temp2.fasta", outfmt=6, out=f"./working/{name}/temp_blast_results.txt")
            stdout, stderr = blastp_cline()

            # 解析BLAST结果，找出最小E-value的比对
            min_e_value = float('inf')
            best_result = None

            with open(f"./working/{name}/temp_blast_results.txt") as result_handle:
                for line in tqdm.tqdm(result_handle):
                    fields = line.strip().split('\t')
                    e_value = float(fields[10])
                    score = float(fields[11])
                    identities = float(fields[2])
                    align_length = int(fields[3])

                    if e_value < min_e_value:
                        min_e_value = e_value
                        best_result = {
                            "Gene1": row[0],
                            "Gene2": gene2,
                            "E-value": e_value,
                            "Score": score,
                            "Identity": identities,
                            "Alignment Length": align_length
                        }

            # 将最小E-value的比对结果添加到blast_results中
            if best_result:
                blast_results.append(best_result)
    # 保存结果到DataFrame
    result_df = pd.DataFrame(blast_results)
    query_lengths = get_gene_lens(f"{name}.fasta", in_folder=f'./working/{name}')
    subject_lengths = get_gene_lens(f"{refname}.fasta", in_folder="./data_available")

    for index, row in result_df.iterrows():
        benchlen = subject_lengths[gene_dict[row["Gene1"]]]
        tarlen = query_lengths[row["Gene2"]]
        cov1 = row["Alignment Length"] / benchlen
        cov2 = row["Alignment Length"] / tarlen
        result_df.loc[index, "COV1"] = cov1
        result_df.loc[index, "COV2"] = cov2

    result_df.to_csv(f"./working/{name}/{name}_ss_{str(threshold)}_blast_results.csv", index=False)
    # 绘制频数分布直方图
    plt.figure(figsize=(10, 6))
    result_df['COV1'].plot(kind='hist', bins=50, edgecolor='black', alpha=0.7)
    result_df['COV2'].plot(kind='hist', bins=50, edgecolor='red', alpha=0.7)
    plt.title(f'Frequency Distribution of Coverage(ss_{threshold})')
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    # plt.xticks([i/10 for i in range(11)])  # 设置x轴刻度
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    # 保存图像到指定位置
    output_path = f'./working/{name}/coverage.png'
    plt.savefig(output_path, format='png', dpi=300)
    return result_df


def get_gene_lens(query, in_folder='prots'):
    file = '%s/%s' % (in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = dict()

    for record in records:
        out[record.id.split('|')[1]] = len(record.seq)

    return out

def make_blast_db(id, folder='prots', db_type='prot'):

    out_file = '%s/%s.pin' % (folder, id)

    cmd_line = 'makeblastdb -in %s/%s -dbtype %s ' % (folder, id, db_type)

    print('making blast db with following command line...')
    print(cmd_line)
    os.system(cmd_line)


def run_blastp(seq, db, in_folder1='prots',in_folder2=None, out="result", outfmt=6, evalue=0.001, threads=1):

    print('blasting %s vs %s' % (seq, db))
    db = '%s/%s' % (in_folder2, db)
    seq = '%s/%s' % (in_folder1, seq)
    cmd_line = 'blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
               % (db, seq, out, evalue, outfmt, threads)

    print('running blastp with following command line...')
    print(cmd_line)
    exit_code = os.system(cmd_line)

    if exit_code == 0:
        print("Command executed successfully!")
    else:
        print(f"Command failed with exit code {exit_code}.")
    return out


def get_gene_lens_frame(query, in_folder='prots'):
    file = '%s/%s' % (in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []

    for record in records:
        out.append({'gene': record.id.split('|')[1], 'gene_length': len(record.seq)})

    out = pd.DataFrame(out)
    return out

def get_bbh(query, subject, in_fold1,in_fold2,cov_thre):
    #     #Utilize the defined protein BLAST function
    run_blastp(query, subject,in_folder1=in_fold1,in_folder2=in_fold2,out="%s/target_bench.txt"%in_fold1)
    run_blastp(subject, query,in_folder1=in_fold2,in_folder2=in_fold1,out="%s/bench_target.txt"%in_fold1)

    query_lengths = get_gene_lens_frame(query, in_fold1)
    subject_lengths = get_gene_lens_frame(subject, in_fold2)

    # Define the output file of this BLAST
    out_file = f'{in_fold1}/bbh_parsed_{str(cov_thre)}.csv'

    # Combine the results of the protein BLAST into a dataframe
    print('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd',
            'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh = pd.read_csv('%s/target_bench.txt' % (in_fold1), sep='\t', names=cols)
    bbh['gene']=bbh['gene'].apply(lambda x: x.split('|')[1])
    bbh['subject'] = bbh['subject'].apply(lambda x: x.split('|')[1])
    bbh = pd.merge(bbh, query_lengths)
    bbh['COV'] = bbh['alnLength'] / bbh['gene_length']

    bbh2 = pd.read_csv('%s/bench_target.txt' % (in_fold1), sep='\t', names=cols)
    bbh2['gene']=bbh2['gene'].apply(lambda x: x.split('|')[1])
    bbh2['subject']=bbh2['subject'].apply(lambda x: x.split('|')[1])
    bbh2 = pd.merge(bbh2, subject_lengths)
    bbh2['COV'] = bbh2['alnLength'] / bbh2['gene_length']

    bbh.to_csv(f"{in_fold1}/bbh_cov.csv")
    bbh2.to_csv(f"{in_fold1}/bbh2_cov.csv")
# define a function to get Bi-Directional BLASTp Best Hits

def bbh_analysis(query, subject, in_fold1,in_fold2,cov_thre):
    # #     #Utilize the defined protein BLAST function
    # run_blastp(query, subject,in_folder1=in_fold1,in_folder2=in_fold2,out="%s/target_bench.txt"%in_fold1)
    # run_blastp(subject, query,in_folder1=in_fold2,in_folder2=in_fold1,out="%s/bench_target.txt"%in_fold1)
    #
    # query_lengths = get_gene_lens_frame(query, in_fold1)
    # subject_lengths = get_gene_lens_frame(subject, in_fold2)
    #
    # # Define the output file of this BLAST
    out_file = f'{in_fold1}/bbh_parsed_{str(cov_thre)}.csv'
    #
    # # Combine the results of the protein BLAST into a dataframe
    # print('parsing BBHs for', query, subject)
    # cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd',
    #         'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    # bbh = pd.read_csv('%s/target_bench.txt' % (in_fold1), sep='\t', names=cols)
    # bbh = pd.merge(bbh, query_lengths)
    # bbh['COV'] = bbh['alnLength'] / bbh['gene_length']
    #
    # bbh2 = pd.read_csv('%s/bench_target.txt' % (in_fold1), sep='\t', names=cols)
    # bbh2 = pd.merge(bbh2, subject_lengths)
    # bbh2['COV'] = bbh2['alnLength'] / bbh2['gene_length']
    out = pd.DataFrame()
    bbh=pd.read_csv(f"{in_fold1}/bbh_cov.csv")
    bbh2=pd.read_csv(f"{in_fold1}/bbh2_cov.csv")
    # Filter the genes based on coverage
    bbh = bbh[bbh.COV >= cov_thre]
    bbh2 = bbh2[bbh2.COV >= cov_thre]
    # bbh.to_csv(f"{in_fold1}/bbh_cov.csv")
    # bbh2.to_csv(f"{in_fold1}/bbh2_cov.csv")
    # Delineate the best hits from the BLAST
    for g in bbh.gene.unique():
        res = bbh[bbh.gene == g]
        if len(res) == 0:
            continue
        best_hit = res.loc[res.PID.idxmax()].copy()
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene == best_gene]
        if len(res2) == 0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        if g == best_gene2:
            best_hit['BBH'] = '<=>'
        else:
            best_hit['BBH'] = '->'
        out = pd.concat([out, pd.DataFrame(best_hit).transpose()])

    # Save the final file to a designated CSV file
    out.to_csv(out_file)
    return out



def bbh(name,refname,threshold,cov_thre=0.25,pid_thre=50,direction="bidirect"):
    yea = pd.read_excel(f"./data_available/{refname}.xlsx")
    gene_dict = yea.set_index('Gene Names (ordered locus)')['Entry'].to_dict()
    gene_dict2 = yea.set_index('Entry')['Gene Names (ordered locus)'].to_dict()
    make_blast_db(f"{name}.fasta", f"./working/{name}")
    #make_blast_db(f"{refname}.fasta", "./working")
    query=f"{name}.fasta"
    subject=f"{refname}.fasta"
    infold1=f"./working/{name}/"
    infold2="./data_available/"
    if os.path.exists(f'./working/{name}/bbh_parsed_{str(cov_thre)}.csv'):
        bbh_df=pd.read_csv(f'./working/{name}/bbh_parsed_{str(cov_thre)}.csv')
    else:
        get_bbh(query, subject, infold1, infold2, cov_thre)
        bbh_df = bbh_analysis(query, subject, infold1, infold2, cov_thre)
    if direction=="bidirect":
        bbh_df = bbh_df[bbh_df["BBH"] == '<=>']
    else:
        pass
    bbh_df = bbh_df[bbh_df["PID"] > pid_thre]
    bbh_df['subject'] = bbh_df['subject'].replace(gene_dict2)
    bbh_df.rename(columns={"subject": 'Gene1', 'gene': 'Gene2'}, inplace=True)
    bbh_df.to_excel(f"./working/{name}/{name}_ss{threshold}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.xlsx")


def cov_filter(name,refname,threshold,cov_threshold,id_threshold):
    df=pd.read_csv(f"./working/{name}/{name}_ss_{str(threshold)}_blast_results.csv")
    df = df[df["COV1"] > cov_threshold]
    df = df[df["COV2"] > cov_threshold]
    df = df[df["Identity"] > id_threshold]
    df.to_csv(f"./working/{name}/{name}_ss_{str(threshold)}_blast_cov{str(cov_threshold)}_id{str(id_threshold)}.csv")

    return df

def merge(name,refname,threshold,cov_threshold,id_threshold,cov_thre,pid_thre,direction):
    bbh_df=pd.read_excel(f"./working/{name}/{name}_ss{threshold}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.xlsx")
    df=pd.read_csv(f"./working/{name}/{name}_ss_{str(threshold)}_blast_cov{str(cov_threshold)}_id{str(id_threshold)}.csv")
    merged_df = pd.merge(bbh_df, df, on=['Gene1', 'Gene2'], how="outer")
    merged_df.to_excel(f"./working/{name}/{name}_ss{threshold}_cov{str(cov_threshold)}_id{str(id_threshold)}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.xlsx", header=True)

# def evaluate(name,refname,threshold,cov_threshold,cov_thre,pid_thre,direction):
#     gxss = pd.read_excel(f"./working/{name}/{name}_ss{threshold}_cov{str(cov_threshold)}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.xlsx")
#     gs = pd.read_excel(f"./ju/juzhen/juzhen_homolog_{name}.xlsx")
#     ssc = set(gxss["Gene1"])
#     benc = set(gs[0])
#     # 定义两个集合
#     set1 = ssc
#     set2 = benc
#
#     # 计算独立部分和交集的数量
#     only_set1 = len(set1 - set2)
#     only_set2 = len(set2 - set1)
#     intersection = len(set1 & set2)
#
#     # 创建Venn图
#     venn = venn2(subsets=(only_set1, only_set2, intersection), set_labels=('ss_predictor', 'benchmark'))
#     plt.title(f"bbh_pid04_cov025({name})")
#     output_path = f'./working/{name}/{name}_ss{threshold}_cov{str(cov_threshold)}_bbh_pid{str(pid_thre)}_cov{str(cov_thre)}_{direction}.png'
#     plt.savefig(output_path, format='png', dpi=300)
#     # 显示图形
#     plt.show()