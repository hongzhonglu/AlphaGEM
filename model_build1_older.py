import cobra
import pandas as pd
import numpy as np
import time


def gx_get(name):
    global gx3
    gx3 = pd.read_excel(f'juzhen/juzhen_homolog{name}.xlsx')
    return gx3
m = 0
allgenes=[]
gx4=[]
def tianjia(rule,gene1,gene2):
    if genefind(rule,gene1)=='':
        return rule
    while True:
        if rule.find(genefind(rule,gene1)+' or ')!=-1:
            rule=rule.replace(genefind(rule,gene1)+' or ','')+' or '+genefind(rule,gene1).replace(gene1,gene2)+' or '+genefind(rule,gene1)
        else:
            if rule.find(' or '+genefind(rule,gene1))!=-1:
                rule = rule.replace(genefind(rule, gene1), '') + genefind(rule, gene1).replace(gene1,gene2) + ' or ' + genefind(rule, gene1)
            else:
                rule ='('+genefind(rule, gene1).replace(gene1,gene2)+ ') or (' + genefind(rule, gene1)+')'
        if not (rule.find(genefind(rule, gene1)) < rule.find(genefind(rule, gene2)) or genefind(rule,
                                                                                         gene2) == '') and genefind(
            rule, gene1) != '':break
    try:
        rule=jianyan(rule,gene2,rule)
    except:
        return rule
    return rule


def jianyan(rule,gene2,rule2):
    if gene2=='':
        return rule
    if type(rule)!=str or type(rule2)!=str:
        return str(rule)
    if rule2.replace(f'{genefind(rule2, gene2)}', '', 1).find(f'{genefind(rule2, gene2)}')==-1:
        return rule
    elif rule2.replace(f'{genefind(rule2, gene2)}', '', 1).find(f'{genefind(rule, gene2)}')!=-1:
        if rule.find(f'{genefind(rule2, gene2)} or ') !=-1:
            rule=rule.replace(f'{genefind(rule2, gene2)} or ', '', 1)
            rule2 = rule2.replace(f'{genefind(rule2, gene2)} or ', '', 1)
        else:
            rule=rule.replace(f' or {genefind(rule2, gene2)}', '', 1)
            rule = rule.replace(f' and {genefind(rule2, gene2)}', '', 1)
            rule = rule.replace(f'{genefind(rule2, gene2)} and ', '', 1)
            rule2 = rule2.replace(f' or {genefind(rule2, gene2)}', '', 1)
            rule2 = rule2.replace(f' and {genefind(rule2,gene2)}','',1)
            rule2 = rule2.replace(f'{genefind(rule2, gene2)} and ', '', 1)
        rule=jianyan(rule,gene2,rule2)
        return rule


def gswitch(rule, gene1,gene2):
    if rule.find(f'{gene1}')==-1 or (rule.find(f'{gene1}')+len(gene1)!=len(rule) and rule[rule.find(f'{gene1}')+len(gene1)]=='-'):
        print(rule)
        return rule
    else:
        rule=rule.replace(f'{gene1}',f'{gene2}')
        print(rule)
        return rule



def genefind(rule, gene):
    a1 = rule.find(f'{gene}')
    if a1==-1:
        return ''
    gf=''
    b1 = rule.rfind('or', 0, a1)
    c1 = rule.find('or', a1)
    if b1 != -1 and c1 == -1:
            gf = rule[b1 + 3:]
    elif b1 == -1 and c1 != -1:
            gf = rule[:c1 - 1]
    elif b1 != -1 and c1 != -1:
        gf = rule[b1+3:c1 - 1]
    elif b1 == -1 and c1 == -1:
        gf = rule
    return gf


def knock(rule, gene):
    if genefind(rule,gene)!='':
        rule = rule.replace(f'{genefind(rule, gene)} or ', '')
    if genefind(rule, gene) != '':
        rule = rule.replace(f' or {genefind(rule, gene)}', '')
    if genefind(rule, gene) != '':
        rule = rule.replace(f'{genefind(rule, gene)}', '')
    return rule


def qiaochu(rule, gene,n):
    if rule=='None':
        print(gene+' - '+rule+' wrong')
    if rule.find(f'{gene}') != -1:
        rule = knock(rule, gene)
        rule = qiaochu(rule, gene,n)
        print(rule + '-' + gene+''+str(n))
        return rule
    elif rule.find(f'{gene}') == -1:
        return rule

def pregenes(cangenes,canmod):
    global gx4
    for a in range(len(canmod.genes)):
        cangenes.append(canmod.genes[a].id)
    for i in range(len(gx3.index)):
        gx4.append(gx3.iat[i, 1])
    return cangenes


def premodel(ymod,canmod,cangenes):
    for i in ymod.genes:
        allgenes.append(i.id)
    for a in range(len(allgenes)):
        pan = 0
        for i in range(len(gx3.index)):
            if allgenes[a] == gx3.iat[i, 1]:
                pan = 1
        if pan == 0:
            cobra.manipulation.remove_genes(ymod, [ymod.genes.get_by_id(allgenes[a])])
    for b in range(len(ymod.reactions)):
        if (ymod.reactions[b].upper_bound != 0) or (ymod.reactions[b].lower_bound != 0):
            canmod.add_reactions([ymod.reactions[b]])
        elif (ymod.reactions[b].gene_reaction_rule == ''):
            canmod.add_reactions([ymod.reactions[b]])
    print(len(canmod.reactions))

def delete_or():
    global canmod
    global cangenes
    global gx4
    for a in range(len(canmod.genes)):
        cangenes.iat[a, 0] = canmod.genes[a].id
    for i in range(len(gx3.index)):
        gx4.append(gx3.iat[i, 1])
    print("start delete the 'or' genes")
    for a in range(len(cangenes.index)):
        if cangenes.iat[a, 0] == 0:
            break
        try:
            gx4.index(cangenes.iat[a, 0])
        except:
            for n1 in range(len(canmod.reactions)):
                canmod.reactions[n1].gene_reaction_rule = str(
                    qiaochu(canmod.reactions[n1].gene_reaction_rule, cangenes.iat[a, 0], n1))

def gpr_tissue_solve():
    for i in range(len(ymod.reactions)):
        a=ymod.reactions[i].gene_reaction_rule.find(' or ')
        if ymod.reactions[i].gene_reaction_rule[:a].count('(')!=ymod.reactions[i].gene_reaction_rule[:a].count(')') and a!=-1:
           ymod.reactions[i].gene_reaction_rule=ymod.reactions[i].gene_reaction_rule[1:-1]

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


def change_gpr(canmod,cangenes):
    print('start change gpr')
    for a in range(len(cangenes)):
        if cangenes[a] == 0:
            break
        pan = 0
        n2 = 0
        i = -1
        while True:
            try:
                i = gx4[i + 1:].index(cangenes[a]) + i + 1
                pan += 1
            except:
                break
            if pan == 1:
                n2 = i
                for n1 in range(len(canmod.reactions)):
                    if cangenes[a] in canmod.reactions[n1].gene_reaction_rule:
                       canmod.reactions[n1].gene_reaction_rule = gswitch(canmod.reactions[n1].gene_reaction_rule,
                                                                      cangenes[a], gx3.iat[i, 2])
            if pan > 1:
                for n1 in range(len(canmod.reactions)):
                    if gx3.iat[n2, 2] in canmod.reactions[n1].gene_reaction_rule:
                       canmod.reactions[n1].gene_reaction_rule = tianjia(canmod.reactions[n1].gene_reaction_rule,
                                                                      gx3.iat[n2, 2], gx3.iat[i, 2])
    return canmod

def modelbuild(refmodel,name=''):
    cangenes = []
    canmod = cobra.Model(f'{name}_model')
    global ymod
    if refmodel=='yeast-GEM.xml' or refmodel=='Sco-GEM.xml' or refmodel=='Human-GEM.xml':
        ymod = cobra.io.read_sbml_model(f'models/{refmodel}')
    if refmodel=='iML1515.json':
        ymod = cobra.io.load_json_model(f'models/{refmodel}')
    gx_get(name)
    start_time=time.time()
    if refmodel=='iML1515.json':gpr_tissue_solve()
    premodel(ymod,canmod,cangenes)
    cangenes = pregenes(cangenes,canmod)
    #delete_or()
    if refmodel!='Sco-GEM.xml':
        canmod=change_gpr(canmod,cangenes)
    else:
        canmod=change_gpr_forsco(canmod,cangenes)
    cobra.io.write_sbml_model(canmod, f'models/tarmod{name}.xml')
    end_time = time.time()
    print('the time it cost is' + str(end_time - start_time) + 's')