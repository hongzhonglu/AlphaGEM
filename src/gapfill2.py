from cobra.flux_analysis import gapfilling
import cobra
def gapfill(name,refmod):
    if refmod=='yeast-GEM.xml':
        refmodel = cobra.io.read_sbml_model(f'models/{refmod}')
    if refmod=='iML1515.json':
        refmodel = cobra.io.load_json_model(f'models/{refmod}')
    model = cobra.io.load_yaml_model(f'models/tarmodel_{name}.yml')
    medium = model.medium
    if refmod == 'yeast-GEM.xml':
        medium['r_1714'] = 1000
    if refmod=='iML1515.json':
        medium['EX_glc__D_e'] = 1000
    model.medium = medium
    gap = gapfilling.GapFiller(model, universal=refmodel, integer_threshold=1e-10, demand_reactions=False,
                               lower_bound=0.05)
    gap.model.tolerance = 1e-10
    gap.model.solver.configuration.tolerances.feasibility = 1e-10
    gap.model.solver.configuration.tolerances.integrality = 1e-10
    gap.model.solver.configuration.tolerances.optimality = 1e-10
    gap.model.slim_optimize()
    solution=gap.fill()
    for solution1 in solution:
        print(solution1)
        model.add_reactions(solution1)
    fba = model.optimize()
    medium = model.medium
    if refmod == 'yeast-GEM.xml':
        medium['r_1714'] = 1
    if refmod=='iML1515.json':
        medium['EX_glc__D_e'] = 10
    model.medium = medium
    cobra.io.write_sbml_model(model, f'models/tarmodel__{name}text.xml')
