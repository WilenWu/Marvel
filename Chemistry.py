#-*- coding:utf8 -*-
import numpy as np
import re
from scipy import constants
from sympy import Symbol
from Marvel import PERIODIC_TABLE,CHEMICAL_BOND,MOLECULE_TABLE,CHEMICAL_REACTION
from Marvel.Physics import *
from Marvel.Geometry import *
import matplotlib.pyplot as plt

#----------------atom
class Atom(Fermion):
    def __init__(self,p,n,e=None,symbol=None,anti=False):
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

        self.element=symbol
        self.atomic_number=p
        self.electron_shells,self.electron_subshells=self.electronic_configuration()
        self.outermost_electron=list(self.electron_shells.values())[-1]

        self.name=PERIODIC_TABLE.loc[self.element,'name']
        self.period=PERIODIC_TABLE.loc[self.element,'period']
        self.group=PERIODIC_TABLE.loc[self.element,'group']

    def electronic_configuration(self):
        subshells=PERIODIC_TABLE.loc[self.element,'1s':'7f'].dropna().astype(int)
        g = subshells.index.map(lambda x: x[0])
        shells=subshells.groupby(g).sum()
        return shells.to_dict(),subshells.to_dict()
        
    def draw(self):
        n = len(self.electron_shells)
        subshells = [i + str(j).join(['^{', '}']) for i, j in self.electron_subshells.items()]
        subshells = '$' + ''.join(subshells) + '$'

        Z = self.protons
        A = self.protons + self.neutrons
        inner = '-' + str(Z) if self.anti else '+' + str(Z)
        formula = '$' + str(A).join(['^{', '}']) + str(Z).join(['_{', '}']) + str(self.element) + '$'

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.text(-3, -0.2, formula, bbox=dict(facecolor='white', edgecolor='white'), fontsize=20)
        ax.text(-3, -4, subshells, bbox=dict(facecolor='white', edgecolor='white'), fontsize=10)
        ax.text(-0.7, -0.2, inner, bbox=dict(alpha=0), fontsize=20)
        draw_circle([0, 0], 1,axes=ax,linewidth=3)

        for i, j in self.electron_shells.items():
            draw_circle([-1, 0], 2 + int(i), -np.pi / 6, np.pi / 6,axes=ax,linewidth=3)
            ax.text(0.8 + int(i), -0.2, str(j), bbox=dict(facecolor='white', edgecolor='white'), fontsize=20)

        ax.set_axis_off()
        ax.set_xlim(-3.5, 2 + n)
        ax.set_ylim(-1 - n, 1 + n)
        plt.show()

    def __repr__(self):
        return 'Atom({p},{n})'.format(p=self.protons,n=self.neutrons)


ATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['atomic number'],x['neutron'],symbol=x.name),axis=1)



#--------------------molecule
class Molecule(Fermion):
    def __init__(self,ID,anti=False):
        '''
        :param formula: str
        '''
        self.ID=ID
        self.formula,self.isomeride=self.extract_formula()
        self.name=self.get_property('name')
        self.standard_density=self.get_property('standard density')
        self.SHC=self.get_property('specific_heat_capacity kJ/(kg*K)')
        self.atoms_dict,self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge()
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy()  #kJ/mol

    def extract_formula(self):
        if self.ID.find('#')>=0:
            return self.ID.split('#')
        else:
            return self.ID,None

    def get_property(self,property):
        'get properties'
        return MOLECULE_TABLE.loc[self.ID,property]

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
        bonds=MOLECULE_TABLE.loc[self.ID,'bond']
        if bonds is None:
            return None,None
        bonds_energy=0
        for i,j in bonds.items():
            bonds_energy+=CHEMICAL_BOND[i]*j
        return bonds,bonds_energy

    def is_compound(self):
        return len(self.atoms_dict)>1

    def elementary_substance_group(self):
        if self.is_compound():
            return 'compound'
        elif self.atoms[0].outermost_electron<4:
            return 'metal'
        elif self.atoms[0].outermost_electron==7:
            return 'halogen'
        elif self.atoms[0].outermost_electron==8:
            return 'noble gas'
        else:
            return 'non-metal'

    def __repr__(self):
        return "Molecule({})".format(self.formula)

class Ion(Molecule):
    def __init__(self,formula,name=None,anti=False):
        '''
        :param formula: str
        '''
        self.ID=formula
        self.name=name
        self.formula,self.isomeride=self.extract_formula()
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

MOLECULES=MOLECULE_TABLE.index.to_series().map(Molecule)

#------------------------Enzymatic Reaction(酶促反应)
class EnzymaticReaction(ChemicalReaction):

    def __init__(self, equation):
        self.equation = equation
        self.reactant, self.product = self.get_stoichiometric_number()

        self.enzyme = self.chemical_property('enzyme')
        self.Km = self.chemical_property('Km(mmol/L)')
        self.Kcat = self.chemical_property('Kcat(1/s)')
        self.catalytic_efficiency = self.Kcat / self.Km

        self.reaction_enthalpies = self.reaction_enthalpies()  # kJ/mol
        self.ΔrHmΘ = self.reaction_enthalpies  # kJ/mol
        self.reaction_heat = self.reaction_enthalpies  # kJ/mol
        self.Qp = self.reaction_heat  # kJ/mol

    def rate_equation(self, Temp, **concentration):
        '''
        molecules,atoms or ions : float(mol/L), cover enzyme(酶),etc
        '''
        c = concentration  # mol/L

        enzyme = self.enzyme
        relative_conc = [c.get(i, 0) / j for i, j in self.reactant.items()]  # 相对浓度
        s = min(relative_conc)  # 底物浓度(mol/L)
        e = c.get(enzyme, 0)  # 酶浓度(mol/L)

        Km = self.Km  # mmol/L
        Kcat = self.catalytic_efficiency * Km  # 1/s
        Km = Km * 1000  # mol/L
        Vmax = Kcat * e  # mol/(L*s)

        v = Vmax * s / (Km + s)
        return v

#-----------------organic compound
'''
类别异构
有机物之间具有以下的类别异构关系：
1. 分子组成符合CnH2n(n≥3)的类别异构体: 烯烃和环烷烃;
2. 分子组成符合CnH2n-2(n≥4)的类别异构体: 炔烃和二烯烃;
3. 分子组成符合CnH2n+2O(n≥3)的类别异构体: 饱和一元醇和饱和醚;
4. 分子组成符合CnH2nO(n≥3)的类别异构体: 饱和一元醛和饱和一元酮;
5. 分子组成符合CnH2nO2(n≥2)的类别异构体: 饱和一元羧酸和饱和一元酯;
6. 分子组成符合CnH2n-6O(n≥7)的类别异构体: 苯酚的同系物,芳香醇及芳香醚;
如n=7,有以下五种: 邻甲苯酚,间甲苯酚,对甲苯酚;苯甲醇;苯甲醚.
7. 分子组成符合CnH2n+1O2N(n≥2)的类别异构体: 氨基酸和硝基化合物
结构特点编辑
有机化合物：种类繁多、数目庞大（已知有3000多万种、且还在以每年数百万种的速度增加）。
但组成元素少 有C、H、O、N 、P、 S、 X（卤素：F、Cl、Br、I ）等
'''
#alkyl group

class AlkylMolecule(Molecule):
    '烃'
    def __init__(self,n,IHD=0,anti=False):
        '''
        :param formula: str
        缺氢指数index of hydrogen deficiency (IHD)
        '''
        pass

#-----------------MacroMolecule(大分子)

class MacroMolecule(Molecule):
    pass


class PolysacCharideMolecule(Molecule):
    '''Carbohydrate(糖类)
    Monosaccharide(单糖,C6H12O6): glucose(葡萄糖),Fructose(果糖),Galactose(半乳糖),ribose(核糖,C5H12O5)
    Disaccharide(二糖,C12H22O11): sugar(蔗糖),maltose(麦芽糖)
    Polysaccharide(多糖,(C6H10O5)n): starch(淀粉),cellulose(纤维素,50000-2500000)
    '''
    def __init__(self, n, type='starch', anti=False):
        '''多糖
        :param n: int
        :param type: starch(s,淀粉),cellulose(c,纤维素,50000-2500000)
        isomeride: s,c
        :return: PolysacCharide
        '''
        self.name=type
        self.n=n if n else Symbol('n')
        self.general_formula = '(C6H10O5)n'+ '#' + self.name[0]
        self.ID=self.general_formula.replace('n',self.n) if n else self.general_formula
        self.formula,self.isomeride=self.extract_formula()

        self.standard_density=self.get_property('standard density')
        self.SHC=self.get_property('specific_heat_capacity kJ/(kg*K)')

        self.atoms_dict, self.atoms = self.get_atoms(self.formula)
        self.atoms_dict={i:j*self.n for i,j in self.atoms_dict.items()}
        self.anti = anti
        self.mass, self.charge = self.get_mass_charge()
        self.chemical_bonds, self.bonds_energy = self.get_bonds_energy()  # kJ/mol

    def get_bonds_energy(self):
        bonds=MOLECULE_TABLE.loc[self.ID,'bond']
        bonds={i:re.search('\d+',str(j)).group() for i,j in bonds.items()}
        bonds={i:float(j)*self.n for i,j in bonds.items()}
        if bonds is None:
            return None,None
        bonds_energy=0
        for i,j in bonds.items():
            bonds_energy+=CHEMICAL_BOND[i]*j
        return bonds,bonds_energy

    @classmethod
    def draw_glucose(cls, center=[0, 0], axes=plt, dehydration=False):
        ax = axes
        x, y = draw_polygon(center, r=1, n=6, axes=ax, linewidth=3)
        for i, j in zip(x[-1:-6:-1], y[-1:-6:-1]):
            ax.plot([i, i], [j - 0.5, j + 0.5])
        ax.plot([x[0], x[0] + 1], [y[0] + 0.5, y[0] + 0.5])
        ax.plot([x[3], x[3] - 1.5], [y[3] + 0.5, y[3] + 0.5])
        ax.text(x[1], y[1], 'O', fontsize=25, ha='center', va='center',
                bbox=dict(boxstyle='circle,pad=0', fc='w', ec='w'))
        ax.text(x[3] - 0.5, y[3] + 0.5, 'O', fontsize=25, ha='center', va='center',
                bbox=dict(boxstyle='circle,pad=0', fc='w', ec='w'))
        sx = [x[2] - 0.2, x[4] - 0.2, x[5] - 0.2]
        sy = [y[2] + 0.5, y[4] + 0.5, y[5] - 0.7]
        s = ['CH$_2$OH', 'OH', 'OH']
        for i, j, k in zip(sx, sy, s):
            ax.text(i, j, k, fontsize=25, bbox=dict(boxstyle='square,pad=0', fc='w', ec='w'))
        if not dehydration:
            ax.text(x[0] + 0.5, y[0] + 0.5, 'OH', fontsize=25, va='center',
                    bbox=dict(boxstyle='square', fc='w', ec='w'))
            ax.text(x[3] - 1.5, y[3] + 0.5, 'H', fontsize=25, va='center', bbox=dict(boxstyle='square', fc='w', ec='w'))
        return x, y

    @classmethod
    def draw_structure(cls, n=None):
        m = n - 1 if n else 'n-1'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x1, y1 = cls.draw_glucose(axes=ax, dehydration=True)
        x2, y2 = cls.draw_glucose(center=[4.5, 0], axes=ax, dehydration=True)
        ax.plot([x1[0] + 0.5, x1[0] + 0.6, x1[0] + 0.6, x1[0] + 0.5], [y1[0] + 1, y1[0] + 1, y1[0] - 1, y1[0] - 1],
                linewidth=2, color='black')
        ax.plot([x1[3] - 0.9, x1[3] - 1, x1[3] - 1, x1[3] - 0.9], [y1[3] + 1, y1[3] + 1, y1[3] - 1, y1[3] - 1],
                linewidth=2, color='black')
        ax.text(x1[0] + 0.6, y1[0] - 1, m, fontsize=25)
        ax.text(x1[3] - 2, y1[3] + 0.5, 'H', va='center', fontsize=25,
                bbox=dict(boxstyle='square,pad=0', fc='w', ec='w'))
        ax.text(x2[0] + 1, y2[0] + 0.5, 'OH', va='center', fontsize=25,
                bbox=dict(boxstyle='square,pad=0', fc='w', ec='w'))
        ax.set_axis_off()
        plt.show()

    def __repr__(self):
        return "{0}: {1}".format(self.name.title(), self.formula)

'''
Protein(蛋白质)
amino acid(氨基酸,R-CH(NH2)-COOH)=amino(氨基,-NH2)+carboxyl(羧基,-COOH)+radical group(基团,R)
'''

class AminoAcidMolecule(Molecule):
    def __init__(self,abbreviation):
        self.abbreviation=abbreviation
        self.ID=self.get_property(self, formula)
        Molecule.__init__(ID)
        self.radical_group=self.formula.split('-')
        self.codon=self.get_property(self, 'codon')


    def get_property(self,property):
        'get physical properties'
        return MOLECULE_TABLE[property].get(self.abbreviation,None)

    def draw_structure(self):
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot([-2,2],[0,0])
        ax.plot([0,0],[-1,1])
        ax.plot([-1,-1],[0,-1])
        ax.plot([0.95,0.95],[0,1])
        ax.plot([1.05,1.05],[0,1])
        x=[-2,-1,0,1,2, 0,1, -1,0]
        y=[0,0,0,0,0, 1,1, -1,-1]
        s=['H','N','C','C','OH', 'R','O', 'H','H']
        for i,j,k in zip(x,y,s):
            ax.text(i,j,k,fontsize=50,ha='center',va='center',bbox=dict(fc='w',ec='w'))
        ax.set_axis_off()
        plt.show()

AMINO_ACID=AMINO_ACID.apply(lambda x:AminoAcidMolecule())



class ProteinMolecule(Molecule):
    pass

'''Lipid(脂质)
Fat(,脂肪): glycerol(甘油) + fatty acids(脂肪酸)
Phospholipid(磷脂)
Sterol(固醇)
'''


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





class OrganicReaction(ChemicalReaction):
    '有机物通用的反应'
    def __init__(self,equation):
        self.equation=equation
        self.reactant,self.product=self.get_stoichiometric_number()

        self.A=self.chemical_property('A')
        self.Ea=self.chemical_property('Ea(kJ/mol)')#kJ/mol

        self.reaction_enthalpies=self.reaction_enthalpies() #kJ/mol
        self.ΔrHmΘ=self.reaction_enthalpies #kJ/mol
        self.reaction_heat=self.reaction_enthalpies #kJ/mol
        self.Qp=self.reaction_heat #kJ/mol



