import cobra
import pandas as pd
def cur(name):
    model = cobra.io.read_sbml_model(f'models/tarmodel__{name}text.xml')
    reactionlist=[]
    for reaction in model.reactions:
        if reaction.id.startswith('rhea'):
            reactionlist.append(reaction)
    fullreaction = []
    lessreaction=[]
    metabolitesless=[]
    for reaction in reactionlist:
        for metabolite in reaction.metabolites:
            if 'R' in metabolite.formula:
                lessreaction.append(reaction)
                metabolitesless.append(metabolite)
            else:
                fullreaction.append(reaction)
    lessreaction=list(set(lessreaction))
    fullreaction=list(set(fullreaction))