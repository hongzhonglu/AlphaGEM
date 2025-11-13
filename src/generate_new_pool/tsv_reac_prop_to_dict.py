import csv
import os
import re
import pickle

def parse_side(side_str):
    # 解析如 "2 MNXM1@MNXD1 + 3 MNXM2@MNXD2"
    chems = []
    if not side_str:
        return chems
    for item in side_str.strip().split('+'):
        part = item.strip()
        if not part:
            continue
        m = re.match(r"([\-\d.]+)?\s*([A-Za-z0-9_]+@[A-Za-z0-9_]+)", part)
        if m:
            stoih = float(m.group(1)) if m.group(1) is not None else 1.0
            idcomp = m.group(2)
            if '@' in idcomp:
                chem_id, comp = idcomp.split('@', 1)
            else:
                chem_id, comp = idcomp, ''
            chems.append({'id': chem_id, 'comp': comp, 'stoih': stoih})
    return chems

def parse_ecnumber(ref_str):
    # 从reference字段最后一段提取ec号，如:  ... 6.3.1.2
    if not ref_str:
        return ''
    matches = re.findall(r"(\d+\.\d+\.\d+\.\d+)", ref_str)
    return matches[-1] if matches else ''

def main():
    tsv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data_available/reac_prop.tsv'))
    data = {}
    with open(tsv_path, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rxn_id = row.get('#ID') or row.get('id') or row.get('reaction')
            eq = row.get('mnx_equation') or row.get('equation')
            ref = row.get('reference')
            is_transport = row.get('is_transport') or row.get('is_transporter') or row.get('classifs')
            # 有的is_transport可能在最后一列，没有内容就是False，有T就True
            is_trans = str(is_transport).strip().upper() == 'T'

            if not rxn_id or not eq or eq.strip() == '=':
                continue

            # 反应式分割
            if '=' not in eq:
                reac_side, prod_side = eq.strip(), ''
            else:
                reac_side, prod_side = eq.split('=', 1)
            reactants = parse_side(reac_side)
            products = parse_side(prod_side)
            ecnumber = parse_ecnumber(ref)

            data[rxn_id] = {
                'reactants': reactants,
                'products': products,
                'is_transporter': is_trans,
                'ecnumber': ecnumber
            }
    print(f"共提取 {len(data)} 条反应。")
    for i, (rid, rdata) in enumerate(data.items()):
        if i < 5:
            print(rid, rdata)
    outpkl = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data_available/reac_prop_dict.pkl'))
    with open(outpkl, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Pickle输出: {outpkl}")

if __name__ == "__main__":
    main()
