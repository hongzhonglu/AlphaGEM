class compound():
    def __init__(self,keggid,pubchem,name='',formula='',charge=0):
        self.keggid=keggid
        self.name=name
        self.pubchem=pubchem
        self.formula=formula
        self.charge=charge

def reac_compound(file):
    1
import pandas as pd
import pubchempy
cpd=pd.read_csv('data_available/pubchem.csv', sep='\t', header=None)
cpd[0]=cpd[0].apply(lambda x:x[:][8:])
cpd[1]=cpd[1].apply(lambda x:x[:][4:])
cpd.to_excel('togivemetanetx.xlsx')
cpds=[compound(row[1],row[0]) for index,row in cpd.iterrows()]
a=pubchempy.get_compounds(pubchempy.Substance.from_sid(3972).cids[0])
pubchempy