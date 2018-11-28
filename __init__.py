import os
import re
import pandas as pd
from math import inf

path=os.path.abspath(__file__)
dirname,filename=os.path.split(path)
Init='\\'.join([dirname,'Init.xlsx'])

# periodic_table
PERIODIC_TABLE=pd.read_excel(Init,sheet_name='periodic_table')
PERIODIC_TABLE.set_index(keys='element', inplace=True)

# chemical bond: ionic bond,covalent bond and metallic bond
CHEMICAL_BOND=pd.read_excel(Init,sheet_name='chemical_bond')
def __reverse(x):
    if x in ['N-H..O','N..H-O']:
        return None
    else:
        x=re.findall('[A-Z][a-z]?|[-=Ξ]',x)
        x.reverse()
        return ''.join(x)
__CB=CHEMICAL_BOND.copy()
__CB['bond']=__CB['bond'].map(__reverse)
CHEMICAL_BOND=CHEMICAL_BOND.append(__CB.dropna()).drop_duplicates()
CHEMICAL_BOND=pd.Series(CHEMICAL_BOND['energy(KJ/mol)'].values,index=CHEMICAL_BOND['bond'])

MOLECULE_TABLE=pd.read_excel(Init,sheet_name='molecule')
MOLECULE_TABLE.loc[:,['bond','ionization']]=MOLECULE_TABLE.loc[:,['bond','ionization']].applymap(eval)
MOLECULE_TABLE.replace('Inf', inf, inplace=True)
MOLECULE_TABLE.set_index(keys=['formula','isomeride'], inplace=True)


#list of reactions
CHEMICAL_REACTION=pd.read_excel(Init,sheet_name='chemical_equation')
CHEMICAL_REACTION.set_index(keys='equation', inplace=True)
