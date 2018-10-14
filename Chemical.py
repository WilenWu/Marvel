import pandas as pd
import numpy as np
import math
from Marvel.Physics import UNITS
from scipy import constants
from sympy import symbols

UNITS.load_prefix_unit(['kJ','km','nm','μm'])
UNITS.load_unit('ATP',30.54,kJ=1,mol=-1)
UNITS.load_unit('calorie',constants.calorie,J=1)
UNITS.loaded

#元素周期表
PERIODIC_TABLE=pd.read_excel('Marvel/chemical.xlsx',sheet_name='periodic_table')
PERIODIC_TABLE.set_index(keys='element',inplace=True)

# chemical bond: ionic bond,covalent bond and metallic bond
MOLECULE_TABLE=pd.read_excel('Marvel/chemical.xlsx',sheet_name='molecule')
MOLECULE_TABLE.loc[:,'bond']=MOLECULE_TABLE.loc[:,'bond'].applymap(eval)
MOLECULE_TABLE.set_index(keys='formula', inplace=True)
CHEMICAL_BOND=pd.read_excel('Marvel/chemical.xlsx',sheet_name='bond')
CHEMICAL_BOND.set_index(keys='bond', inplace=True)
#--------------------molecule
class Molecule:
    def __init__(self,formula,ionic_electron=0):
        '''
        :param atoms: int
        :param ionic_electron: int
        '''
        self.formula=forluma
        self.atoms=get_atoms(formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge(self.atoms,ionic_electron,anti)
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy(formula)

    @staticmethod
    def get_atoms(formula):
        import re
        def split_formula(formula):
            split_form = re.findall('[A-Z][a-z]?|\(.+\)|\d+', formula)
            atoms = {}
            for i, key in enumerate(split_form):
                if key.isnumeric():
                    next
                elif i + 1 == len(split_form):
                    atoms[key] = atoms.get(key, 0) + 1
                elif split_form[i + 1].isnumeric():
                    value = int(split_form[i + 1])
                    atoms[key] = atoms.get(key, 0) + value
                else:
                    atoms[key] = atoms.get(key, 0) + 1
            return atoms

        atoms = split_formula(formula)

        for key in list(atoms):
            if not key.isalpha():
                value = atoms.pop(key)
                key = key.replace('(', '').replace(')', '')
                sub_atoms = split_formula(key)
                for key in list(sub_atoms):
                    atoms[key] = atoms.get(key, 0) + sub_atoms[key]*value

        return atoms

    def get_mass_charge(self,atoms,ionic_electron=0,anti=False):
        electron=Electron(anti)
        ionic_charge=electron.charge*ionic_electron
        ionic_mass=electron.mass*ionic_electron

        if not set(atoms.keys()).issubset(set(PERIODIC_TABLE.index)):
            return None,None

        mass=charge=0
        for i in atoms.keys():
            p,n=PERIODIC_TABLE.loc[i,['proton','neutron']]
            atom=Atom(p,n,anti=anti)
            mass+=atom.mass
            charge+=atom.charge

        return mass+ionic_mass,charge+ionic_charge

    def get_bonds_energy(self,formula):
        bonds=MOLECULE_TABLE.loc[formula,'bond']
        bonds_energy=0
        for i in bonds:
            bonds_energy+=CHEMICAL_BOND.loc[i,'energy(KJ/mol)']*bonds[i]

        return bonds,bonds_energy

    #电离
    def ionization(self):
        pass

    #离子平衡
    def ion_balance(self):
        pass



# 焓变
def enthalpy(reactant):
    energy=0
    for i in reactant:
        energy+=MOLECULE.loc[i,'bond_energy']*reactant[i]
    return energy

#----------------- Chemical_reaction
'''
Chemical reaction:
combination,A+B=AB(化合反应)
decomposition,AB=A+B(分解反应)
single displacement,AB+C=A+BC(置换反应)
double decomposition,AB+CD=AC+BD(复分解反应)
'''
REACTION = pd.DataFrame(data=[
    [{'H2': 2, 'O2': 1}, {'H2O': 2}, '2H2+O2=2H2O', {'temperature': 574}, 1],
    [{'N2': 1, 'O2': 1}, {'NO': 2}, 'N2+O2=2NO', {'thunder': True}, 1],
    [{'H2O': 2},{'H2': 2, 'O2': 1},'2H2O=2H2+O2', {'electricity': True}, 0.1],
    [{'CH4': 1, 'O2': 1}, {'NO': 2}, 'CH4+O2=CO2+2H2O', {'temperature': 538}, 1]
],columns=['reactant', 'product', 'equation', 'condition', 'reaction_rate']
)
REACTION['enthalpy_change'] = REACTION['product'].map(enthalpy) - REACTION['reactant'].map(enthalpy)

REACTION2= pd.DataFrame(data=[
    [{'C': 1, 'O2': 1}, {'CO2': 1}, 'C+O2=CO2', {'temperature': 574}, 0.2,328],
    [{'P4': 1, 'O2': 5}, {'P2O5': 2}, '2H2+O2=2H2O', {'temperature': 30}, 3093.2],
    [{'S': 1, 'O2': 1}, {'SO2': 2}, 'S+O2=SO2', {'temperature': 574}, 1,93],
    [{'Mg': 2, 'O2': 1}, {'MgO': 2}, '2Mg+O2=2MgO', {'temperature': 574}, 1,504]
],columns=['reactant', 'product', 'equation', 'condition', 'reaction_rate','enthalpy_change']
)
# REACTION=pd.concat([REACTION,REACTION2,biochemical], axis=1)
REACTION.set_index(keys='equation', inplace=True)
REACTION['ATPs'] = REACTION['enthalpy_change'] / CONSTANT.ATP

class Chemical_reaction:
    def __init__(self,**environment):
        self.delta_t=None

    # 反应速率，即时属性
    def reaction_rate(**environment):
        '''
        molecules,atoms or ions : float(mol/m3)
        reactant cover catalyst(催化剂),enzyme(酶),etc
        environment: temperature,thunder,etc
        '''
        m={}
        energy={}
        for k in REACTION.index:
            reactant_rate = [environment.get(i, 0)  for i, j in REACTION.loc[k]['reactant'].items()]
            condition_rate = [environment.get(i, 0) / j for i, j in REACTION.loc[k]['condition'].items()]
            condition_rate = [(j if j >= 1 else 0) for j in condition_rate]
            rate=np.prod(reactant_rate) * np.prod(condition_rate) * REACTION.loc[k]['reaction_rate']

            reactant_power={i:-1*rate*j for i, j in REACTION.loc[k]['reactant'].items()}
            product_power={i:rate*j for i, j in REACTION.loc[k]['product'].items()}
            energy_power={'enthalpy_power':rate*REACTION.loc[k]['enthalpy_change'],'ATPs_power':rate*REACTION.loc[k]['ATPs']}

            power=reactant_power.copy()
            power.update(product_power)
            m[k]=power
            energy[k]=energy_power

        return pd.DataFrame(m).fillna(0),pd.DataFrame(energy).fillna(0)
