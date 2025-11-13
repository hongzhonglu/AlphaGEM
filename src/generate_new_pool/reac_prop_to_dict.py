import csv
import os
import re
import pickle

def parse_side(side_str):
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
            chem_id, comp = idcomp.split('@', 1)
            chems.append({'id': chem_id, 'comp': comp, 'stoih': stoih})
    return chems

def parse_ecnumber(reference):
    if not reference:
        return ''
    match = re.search(r"(\d+\.\d+\.\d+\.\d+)", reference)
    return match.group(1) if match else ''

def main():
    tsv_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data_available/reac_prop.tsv'))
    output_dict = {}
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rxn_id = row.get('#ID', '').strip()
            equation = row.get('mnx_equation', '').strip()
            reference = row.get('reference', '').strip()
            is_transport = row.get('is_transport', '').strip().upper()

            if not rxn_id or not equation or equation == '=':
                continue

            if '=' not in equation:
                left, right = equation.strip(), ''
            else:
                left, right = equation.split('=', 1)

            reactants = parse_side(left)
            products = parse_side(right)
            ecnumber = parse_ecnumber(reference)
            transporter = (is_transport == 'T')

            output_dict[rxn_id] = {
                'reactants': reactants,
                'products': products,
                'is_transporter': transporter,
                'ecnumber': ecnumber
            }
    print(f"共提取 {len(output_dict)} 条反应样本。示例：")
    for i, (k, v) in enumerate(output_dict.items()):
        if i < 5:
            print(k, v)
    outpkl = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../data_available/reac_prop_dict.pkl'))
    with open(outpkl, 'wb') as f:
        pickle.dump(output_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Pickle已保存: {outpkl}")

if __name__ == "__main__":
    main()
