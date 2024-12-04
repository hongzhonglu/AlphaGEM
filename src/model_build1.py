import cobra
import pandas as pd
import numpy as np
import time


def gx_get(name):
    global gx3
    gx3 = pd.read_excel(f'working/{name}/matrix_homolog{name}.xlsx')
    return gx3
m = 0
allgenes=[]
gx4=[]


def gpr_tissue_solve():
    for i in range(len(ymod.reactions)):
        a=ymod.reactions[i].gene_reaction_rule.find(' or ')
        if ymod.reactions[i].gene_reaction_rule[:a].count('(')!=ymod.reactions[i].gene_reaction_rule[:a].count(')') and a!=-1:
           ymod.reactions[i].gene_reaction_rule=ymod.reactions[i].gene_reaction_rule[1:-1]
    return ymod


def premodel(ymod,canmod,cangenes):
    allgenes=[]
    delgene=[]
    for i in ymod.genes:
        allgenes.append(i.id)
    homogenes=list(gx3['refmodelgene'])
    for i in allgenes:
        if i not in homogenes:
            delgene.append(ymod.genes.get_by_id(i))
    cobra.manipulation.remove_genes(ymod, delgene)
    for b in range(len(ymod.reactions)):
        if (ymod.reactions[b].upper_bound != 0) or (ymod.reactions[b].lower_bound != 0):
            canmod.add_reactions([ymod.reactions[b]])
        elif (ymod.reactions[b].gene_reaction_rule == ''):
            canmod.add_reactions([ymod.reactions[b]])
    print(len(canmod.reactions))
    return canmod


def pregenes(cangenes,canmod):
    global gx4
    for a in range(len(canmod.genes)):
        cangenes.append(canmod.genes[a].id)
    for i in range(len(gx3.index)):
        gx4.append(gx3.iat[i, 1])
    return cangenes


def change_gpr_forsco(canmod,cangenes):
    print('start change gpr')
    for a in range(len(cangenes)):
        if cangenes[a] == 0:
            break
        pan = 0
        n2 = 0
        i = -1
        while True:
            rule_pre=''
            try:
                i = gx4[i + 1:].index(cangenes[a]) + i + 1
                pan += 1
            except:
                break
            if pan == 1:
                n2 = i
                rule_pre+=gx3.iat[i, 2]
            if pan > 1:
                rule_pre=rule_pre+' or '+gx3.iat[i, 2]
            for n1 in range(len(canmod.reactions)):
                if cangenes[a] in canmod.reactions[n1].gene_reaction_rule:
                    canmod.reactions[n1].gene_reaction_rule = canmod.reactions[n1].gene_reaction_rule.replace(cangenes[a],'('+rule_pre+')')
    return canmod


def str_chanege_for_gpr(old_gpr,oldgene,newgenes):
    new_gpr=''
    oldgprcell=old_gpr.split(' or ')
    gprcell=[]
    for a in oldgprcell:
        if oldgene in a:
            for genes in newgenes:
                if 'and' in a:
                    gprcell.append('('+a.replace(oldgene,genes).replace('(','').replace(')','')+')')
                else:
                    gprcell.append(a.replace(oldgene, genes).replace('(', '').replace(')', ''))
        else:
            if 'and' in a:
                gprcell.append('('+a.replace('(','').replace(')','')+')')
            else:
                gprcell.append(a.replace('(', '').replace(')', ''))
    new_gpr=' or '.join(list(set(gprcell)))
    return new_gpr


def change_gpr(canmod,cangenes,dict):
    print('start change gpr')
    for refgene in dict.keys():
        try:
          for reaction in canmod.genes.get_by_id(refgene).reactions:
            reaction.gene_reaction_rule=str_chanege_for_gpr(reaction.gene_reaction_rule,refgene,dict[refgene])
        except:
            continue
    return canmod



def modelbuild(refmodel,name=''):
    cangenes = []
    canmod = cobra.Model(f'{name}_model')
    global ymod
    if refmodel=='yeast-GEM.xml' or refmodel=='Sco-GEM.xml' or refmodel=='Human-GEM.xml':
        ymod = cobra.io.read_sbml_model(f'models/{refmodel}')
    if refmodel=='iML1515.json':
        ymod = cobra.io.load_json_model(f'models/{refmodel}')
    gx3=gx_get(name)
    gx3.columns=['number','refmodelgene','tarmodelgene']
    dict={}
    for groups in gx3.groupby('refmodelgene'):
        dict[groups[0]]=list(groups[1]['tarmodelgene'])
    start_time=time.time()
    if refmodel=='iML1515.json':ymod=gpr_tissue_solve()
    canmod=premodel(ymod,canmod,cangenes)
    cangenes = pregenes(cangenes,canmod)
    #delete_or()
    if refmodel!='Sco-GEM.xml':
        canmod=change_gpr(canmod,cangenes,dict)
    else:
        canmod=change_gpr_forsco(canmod,cangenes)
    cobra.io.save_yaml_model(canmod, f'working/{name}/model1{name}.yml')
    end_time = time.time()
    print('the time it cost is' + str(end_time - start_time) + 's')