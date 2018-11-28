#-*- coding:utf8 -*-
import numpy as np
import re
from scipy import constants
from Marvel import PERIODIC_TABLE,CHEMICAL_BOND,MOLECULE_TABLE,CHEMICAL_REACTION
from Marvel.Physics import *

#----------------atom
class Atom(Fermion):
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


ATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron']),axis=1)
antiATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron'],anti=True),axis=1)

#--------------------molecule
class Molecule(Fermion):
    def __init__(self,formula,isomeride='O',anti=False):
        '''
        :param formula: str
        '''
        self.formula=formula
        self.isomeride=isomeride
        self.name=get_property('name')
        self.standard_density=get_property('standard density')
        self.SHC=get_property('specific_heat_capacity kJ/(kg*K)')
        self.atoms_dict,self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge()
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy()  #kJ/mol

    def get_property(self,property):
        'get properties'
        return MOLECULE_TABLE.loc[self.formula,self.isomeride].loc[property]

    @staticmethod
    def get_atoms(formula):
        def split_formula(formula):
            split_form = re.findall('[A-Z][a-z]?|\(.+\)|\d+', formula)
            atoms_dict = {}
            for i, key in enumerate(split_form):
                if key.isnumeric() or key in list('^+-'):
                    next
                elif i + 1 == len(split_form):
                    atoms_dict[key] = atoms_dict.get(key, 0) + 1
                elif split_form[i + 1].isnumeric():
                    value = int(split_form[i + 1])
                    atoms_dict[key] = atoms_dict.get(key, 0) + value
                else:
                    atoms_dict[key] = atoms_dict.get(key, 0) + 1
            return atoms_dict

        atoms_dict = split_formula(formula)

        for key in list(atoms_dict):
            if not key.isalpha():
                value = atoms_dict.pop(key)
                key = key.replace('(', '').replace(')', '')
                sub_atoms = split_formula(key)
                for key in list(sub_atoms):
                    atoms_dict[key] = atoms_dict.get(key, 0) + sub_atoms[key]*value

        atoms=[ATOMS[i].copy() for i in atoms_dict.keys()]
        return atoms_dict,atoms

    def get_mass_charge(self,ionic_electron=0):
        anti=self.anti
        atoms_dict=self.atoms_dict

        electron=Electron(anti)
        ionic_charge=electron.charge*ionic_electron
        ionic_mass=electron.mass*ionic_electron

        if not set(atoms_dict.keys()).issubset(set(ATOMS.index)):
            return None,None

        mass=charge=0
        for i,j in atoms_dict.items():
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
        self.atoms_dict_dict,self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge(self.ionic_electron)
    def __ionic_electron(self):
        charge=re.findall('\^\d+[+-]|[+-]',self.formula)[0].replace('^','')
        sign=-1 if charge[-1]=='+' else 1
        n=1 if len(charge)==1 else int(charge[:-1])
        return n*sign

    def __repr__(self):
        return "Ion({})".format(self.formula)

MOLECULES=MOLECULE_TABLE.index.to_frame().apply(lambda x:Molecule(x['formula'],x['isomeride']),axis=1)
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

        reactant,product=equation.split('=')
        reactant={i:-j for i,j in equation_split(reactant).items()}
        return reactant,equation_split(product)

    def reaction_enthalpies(self):
        def bolds_energy(molecule_dict):
            bolds_energy=0.0
            for i in molecule_dict:
                bolds_energy+=MOLECULES[i].bonds_energy*molecule_dict[i]
            return bolds_energy #kJ/mol

        return -bolds_energy(self.reactant)-bolds_energy(self.product)

    def rate_constant(self,Temp):
        '''
        k = A*exp(-Ea/(RT))
        '''
        R = constants.gas_constant/1000 # kJ/(mol*K)
        A,Ea=self.A,self.Ea #kJ/mol
        return A*np.exp(-Ea/(R*Temp))

    def rate_equation(self,Temp,env):
        '''
        env(dict): reactant concentration, float(mol/L)
        aA+bB=cC+dD
        v = k([A]*[B])**2
        '''
        k=self.rate_constant(Temp)
        v=np.float64(k)
        for i,j in self.reactant.items():
            v=v*np.float64(env.get(i,0))**2
        return v # mol/(L*s)

    def __repr__(self):
        return "ChemicalReaction({})".format(self.equation)

CHEMICAL_REACTION=CHEMICAL_REACTION.index.to_series().map(ChemicalReaction)
