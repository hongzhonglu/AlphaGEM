import pandas as pd
from src.reference_mapping import map_reference_pairs


def ss_predict_choose(name, refname, threshold):
    file_path = f"./working/{name}/{name}_ss_predictor_filter{threshold}"
    pairs = []
    with open(file_path) as handle:
        for line in handle:
            query, target, score = line.rstrip().split("\t")
            if float(score) <= threshold:
                continue
            query = query.split()[0].split("|")
            target = target.split()[0].split("|")
            pairs.append((target[1] if len(target) > 1 else target[0],
                          query[1] if len(query) > 1 else query[0]))
    frame = pd.DataFrame(pairs, columns=[0, 1]).drop_duplicates()
    if refname != "pan":
        frame = map_reference_pairs(
            frame, 0, 1, refname,
            f'working/{name}/{name}_ss_homolog{threshold}_unmapped_reference.json',
        )
    frame.to_excel(f"./working/{name}/{name}_ss_homolog{threshold}.xlsx")


def ss_predict_choose_rhea(name,threshold,threshold1=0.9):
    juzhen_ss = pd.DataFrame()
    file_path=f"./working/{name}/{name}_rhea_search_filter{str(threshold1)}"
    with open(file_path,"r") as f:
        while True:
            line=f.readline()
            if line=="":
                break
            # try:
            #     query=line.strip().split("\t")[0].split("|")[1]
            # except:
            #     print("a",line,"b",count)
            # count+=1
            query = line.strip().split("\t")[0].split(" ")[0]
            tar=line.strip().split("\t")[1].split("|")[1]
            score = float(line.strip().split("\t")[2])
            if score>threshold:
                jp=pd.concat([pd.DataFrame([tar]),pd.DataFrame([query]),pd.DataFrame([score])],axis=1)
                juzhen_ss=pd.concat([jp,juzhen_ss])

    gx2 = pd.DataFrame()
    for i in range(len(juzhen_ss.index)):
        gx2 = pd.concat([gx2, pd.DataFrame({
            0: [juzhen_ss.iat[i, 0]],
            1: [juzhen_ss.iat[i, 1]]
        })])
    # if refname!="pan":
    #     yea = pd.read_excel(f"./ziyuan/{refname}.xlsx")
    #     yea1 = []
    #     yea2 = []
    #     for i in range(len(yea.index)):
    #         yea1.append(yea.iat[i, 0])
    #         yea2.append(yea.iat[i, 2])
    #     for i in range(len(gx2.index)):
    #         index = yea1.index(gx2.iat[i, 0])
    #         gx2.iat[i, 0] = yea2[index]
    rhea_geng=pd.read_csv("./rhea/rhea2uniprot_sprot.tsv")

    gx2.to_excel(f"./working/{name}/{name}_rhea_homolog{str(threshold)}.xlsx")
