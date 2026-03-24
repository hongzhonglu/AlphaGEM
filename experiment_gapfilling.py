import sys
import os
import shutil
import cobra

sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
import gapfilling
import gapfilling_universal

def main():
    # 设置测试的物种和参考模型参数
    # 如果想测试其他物种，可以在这里修改
    name = "spombe"
    refname = "yeast"
    refmodelname = "yeast-GEM.xml"
    grothmedium = "min"
    
    print(f"Starting experiment for {name} using reference {refname}...")
    
    # 确保 working 目录下有必要的文件
    if not os.path.exists(f"working/{name}/{name}-GEM_withgaps.yml"):
        print(f"Error: working/{name}/{name}-GEM_withgaps.yml not found. Please run the pipeline until gapfilling step first.")
        return

    # 1. 运行 标准的 gapfilling
    print("\n--- Running gapfilling.py ---")
    gapfilling.gapfill(name, refname, refmodelname, grothmedium)
    # 保存结果，以免被覆盖
    std_xml = f"working/{name}/{name}-GEM_standard.xml"
    shutil.copy(f"working/{name}/{name}-GEM.xml", std_xml)
    print(f"Standard gapfilling model saved to {std_xml}")
    
    # 2. 运行 universal 的 gapfilling
    print("\n--- Running gapfilling_universal.py ---")
    gapfilling_universal.gapfill(name, refname, refmodelname, grothmedium)
    # 保存结果
    uni_xml = f"working/{name}/{name}-GEM_universal.xml"
    shutil.copy(f"working/{name}/{name}-GEM.xml", uni_xml)
    print(f"Universal gapfilling model saved to {uni_xml}")
    
    # 3. 比较两个模型
    print("\n--- Comparing the two models ---")
    base_model_path = f"working/{name}/{name}-GEM_withgaps.yml"
    model_base = cobra.io.load_yaml_model(base_model_path)
    model_std = cobra.io.read_sbml_model(std_xml)
    model_uni = cobra.io.read_sbml_model(uni_xml)
    
    rxns_base = set([r.id for r in model_base.reactions])
    rxns_std = set([r.id for r in model_std.reactions])
    rxns_uni = set([r.id for r in model_uni.reactions])
    
    genes_std = set([g.id for g in model_std.genes])
    genes_uni = set([g.id for g in model_uni.genes])
    
    mets_std = set([m.id for m in model_std.metabolites])
    mets_uni = set([m.id for m in model_uni.metabolites])
    
    gapfilled_std = rxns_std - rxns_base
    gapfilled_uni = rxns_uni - rxns_base
    
    print(f"Base model (with gaps): {len(rxns_base)} reactions")
    print(f"Standard model: {len(rxns_std)} reactions ({len(gapfilled_std)} reactions were gapfilled)")
    print(f"Universal model: {len(rxns_uni)} reactions ({len(gapfilled_uni)} reactions were gapfilled)")
    if len(gapfilled_uni) > len(gapfilled_std):
        print(f"Universal gapfilling added {len(gapfilled_uni) - len(gapfilled_std)} MORE reactions than Standard gapfilling.")
    elif len(gapfilled_std) > len(gapfilled_uni):
        print(f"Standard gapfilling added {len(gapfilled_std) - len(gapfilled_uni)} MORE reactions than Universal gapfilling.")
    else:
        print(f"Both methods added the same number of reactions ({len(gapfilled_std)}).")
        
    print(f"\nStandard model: {len(genes_std)} genes, {len(mets_std)} metabolites")
    print(f"Universal model: {len(genes_uni)} genes, {len(mets_uni)} metabolites")
    
    # 分析差异的反应
    only_uni = rxns_uni - rxns_std
    only_std = rxns_std - rxns_uni
    
    print(f"\nReactions ONLY in Universal ({len(only_uni)}):")
    for r in sorted(list(only_uni)):
        rxn = model_uni.reactions.get_by_id(r)
        print(f"  {r}: {rxn.name} [{rxn.reaction}]")
        
    print(f"\nReactions ONLY in Standard ({len(only_std)}):")
    for r in sorted(list(only_std)):
        rxn = model_std.reactions.get_by_id(r)
        print(f"  {r}: {rxn.name} [{rxn.reaction}]")
        
    # 分析差异的基因
    only_uni_genes = genes_uni - genes_std
    only_std_genes = genes_std - genes_uni
    
    print(f"\nGenes ONLY in Universal ({len(only_uni_genes)}):")
    for g in sorted(list(only_uni_genes)):
        print(f"  {g}")
        
    print(f"\nGenes ONLY in Standard ({len(only_std_genes)}):")
    for g in sorted(list(only_std_genes)):
        print(f"  {g}")

    # 获取目标值
    print("\n--- Optimization Results ---")
    obj_std = model_std.slim_optimize()
    obj_uni = model_uni.slim_optimize()
    print(f"Objective value (Standard): {obj_std}")
    print(f"Objective value (Universal): {obj_uni}")

if __name__ == '__main__':
    main()
