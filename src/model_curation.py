import cobra
def curate_model(name,refname):
    model=cobra.io.read_sbml_model(f"./working/{name}/{name}-GEM.xml")
    for r in model.boundary:
        print(r.compartments)
    comps=[]
    for r in model.reactions:
        if r.reactants == [] or r.products == []:
            list(r.compartments)
            cc=''
            for i in list(r.compartments):
                cc+=i
            comps.append(cc)
            if cc=='c':
                print(r.name)
    set(comps)