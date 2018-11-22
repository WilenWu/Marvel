#-*- coding:utf8 -*-

import pandas as pd
import numpy as np
import re
from scipy import constants
from math import inf

#-----------------subatomic particles

class Electron:
    def __init__(self,anti=False):
        self.charge=1 if anti else -1
        self.mass=constants.m_e/constants.atomic_mass
    def __repr__(self):
        return 'Electron'

class Proton:
    def __init__(self, anti=False):
        self.charge = -1 if anti else 1
        self.mass=constants.m_p/constants.atomic_mass
    def __repr__(self):
        return 'Proton'

class Neutron:
    def __init__(self, anti=False):
        self.charge=0
        self.mass=constants.m_u/constants.atomic_mass
    def __repr__(self):
        return 'Neutron'

class Photon:
    def __init__(self,frequency):
        self.frequency=frequency
        self.c=constants.c
        self.wavelength=self.c/self.frequency

        self.energy=constants.h*self.frequency
        self.momentum=self.energy/self.c

        self.ν=self.frequency
        self.λ=self.wavelength
        self.E=self.energy
        self.h=self.momentum

    def __repr__(self):
        return 'Photon(ν={})'.format(self.frequency)

class Field:
    pass

class ElectricField(Field):
    pass

class PhotonField(Field):
    pass


#----------------atom
class Atom:
    def __init__(self,p,n,e=None,anti=False):
        e = p if e is None else e
        self.extra_nuclear_electron=e #核外电子数

        proton = Proton(anti)
        neutron = Neutron(anti)
        electron = Electron(anti)

        self.protons = p
        self.neutrons = n
        self.electrons = e

        self.nucleus_mass=proton.mass*p + neutron.mass*n
        self.mass = self.nucleus_mass + electron.mass*e

        self.nuclear_charge=proton.charge*p #核电荷数
        self.charge = proton.charge*p + electron.charge*e

        self.anti=anti

    def nuclear_reaction(self):
        'fission and fusion(裂变和聚变)'
        pass

    def __repr__(self):
        return 'Atom({p},{n})'.format(p=self.protons,n=self.neutrons)


#----------------periodic_table
PERIODIC_TABLE=pd.read_excel('Init.xlsx',sheet_name='periodic_table')
PERIODIC_TABLE.set_index(keys='element', inplace=True)
ATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron']),axis=1)
antiATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron'],anti=True),axis=1)

#--------------------molecule
class Molecule:
    def __init__(self,formula,anti=False):
        '''
        :param formula: str
        '''
        self.formula=formula
        self.name=MOLECULE_TABLE['name'].get(self.formula,None)
        self.standard_density=MOLECULE_TABLE['standard density'].get(self.formula,None)
        self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge()
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy()  #kJ/mol

    @staticmethod
    def get_atoms(formula):
        def split_formula(formula):
            split_form = re.findall('[A-Z][a-z]?|\(.+\)|\d+', formula)
            atoms = {}
            for i, key in enumerate(split_form):
                if key.isnumeric() or key in list('^+-'):
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

    def get_mass_charge(self,ionic_electron=0):
        anti=self.anti
        atoms=self.atoms

        electron=Electron(anti)
        ionic_charge=electron.charge*ionic_electron
        ionic_mass=electron.mass*ionic_electron

        if not set(atoms.keys()).issubset(set(ATOMS.index)):
            return None,None

        mass=charge=0
        for i,j in atoms.items():
            atom=ATOMS[i] if not anti else antiATOMS[i]
            mass+=atom.mass*j
            charge+=atom.charge*j

        return mass+ionic_mass,charge+ionic_charge

    def get_bonds_energy(self):
        bonds=MOLECULE_TABLE.loc[self.formula,'bond']
        if bonds is None:
            return None,None
        bonds_energy=0
        for i,j in bonds.items():
            bonds_energy+=CHEMICAL_BOND[i]*j
        return bonds,bonds_energy

    def __repr__(self):
        return "Molecule({})".format(self.formula)

class Ion(Molecule):
    def __init__(self,formula,name=None,anti=False):
        '''
        :param formula: str
        '''
        self.formula=formula
        self.name=name
        self.ionic_electron=self.__ionic_electron()
        self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge(self.ionic_electron)
    def __ionic_electron(self):
        charge=re.findall('\^\d+[+-]|[+-]',self.formula)[0].replace('^','')
        sign=-1 if charge[-1]=='+' else 1
        n=1 if len(charge)==1 else int(charge[:-1])
        return n*sign

    def __repr__(self):
        return "Ion({})".format(self.formula)

# chemical bond: ionic bond,covalent bond and metallic bond
CHEMICAL_BOND=pd.read_excel('Init.xlsx',sheet_name='chemical_bond')
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

MOLECULE_TABLE=pd.read_excel('Init.xlsx',sheet_name='molecule')
MOLECULE_TABLE.loc[:,['bond','ionization']]=MOLECULE_TABLE.loc[:,['bond','ionization']].applymap(eval)
MOLECULE_TABLE.replace('Inf', inf, inplace=True)
MOLECULE_TABLE.set_index(keys='formula', inplace=True)

MOLECULES=MOLECULE_TABLE.index.to_series().map(Molecule)
antiMOLECULES=MOLECULE_TABLE.index.to_series().map(lambda x:Molecule(x,anti=True))
#----------------Chemical_reaction
class ChemicalReaction:
    def __init__(self,equation):
        self.equation=equation
        self.reactant,self.product=self.get_stoichiometric_number()

        self.A=self.chemical_property('A')
        self.Ea=self.chemical_property('Ea(kJ/mol)')#kJ/mol

        self.reaction_enthalpies=self.reaction_enthalpies() #kJ/mol
        self.ΔrHmΘ=self.reaction_enthalpies #kJ/mol
        self.reaction_heat=self.reaction_enthalpies #kJ/mol
        self.Qp=self.reaction_heat #kJ/mol

    def chemical_property(self,property):
        'get chemical properties'
        return CHEMICAL_REACTION[property].get(self.equation,None)

    def get_stoichiometric_number(self):
        equation=self.equation
        def equation_split(x):
            y={}
            for m in x.split('+'):
                if m.find('*')>=0:
                    value,key=m.split('*')
                    y[key]=int(value)
                else:
                    y[m]=1
            return y

        reacant,product=equation.split('=')
        return equation_split(reacant),equation_split(product)

    def reaction_enthalpies(self):
        def bolds_energy(molecule_dict):
            bolds_energy=0.0
            for i in molecule_dict:
                bolds_energy+=MOLECULES[i].bonds_energy*molecule_dict[i]
            return bolds_energy #kJ/mol

        return bolds_energy(self.reactant)-bolds_energy(self.product)

    def rate_equation(self,Temp,concentration):
        '''
        molecules,atoms or ions : float(mol/L), cover catalyst(催化剂),etc
        '''
        c = concentration  # mol/L

        R = constants.gas_constant/1000 # kJ/(mol*K)
        A,Ea=self.A,self.Ea #kJ/mol
        k = A*np.exp(-Ea/(R*Temp))
        tmp=1
        for i,j in self.reactant.items():
            tmp*=c.get(i,0)**j
        v=k*tmp # mol/(cm3*s)
        return v*1000 # mol/(L*s)

    def __repr__(self):
        return "ChemicalReaction({})".format(self.equation)


#list of reactions
CHEMICAL_REACTION=pd.read_excel('Init.xlsx',sheet_name='chemical_equation')
CHEMICAL_REACTION.set_index(keys='equation', inplace=True)
CHEMICAL_REACTION=CHEMICAL_REACTION.index.to_series().map(ChemicalReaction)
