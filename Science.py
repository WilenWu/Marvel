#-*- coding:utf8 -*-

import pandas as pd
import numpy as np
import re

from scipy import constants,linalg
from sympy import symbols
from math import pi,e,exp,inf
from itertools import combinations,permutations,product

# System International(SI)
class Units:
    'unit: a standard of measurement'
    base=pd.DataFrame(
       columns=['Unit','Symbol','Physical quantities'],
        data=[
            ['meter', 'm', 'length'],
            ['kilogram', 'kg', 'mass'],
            ['second', 's', 'time'],
            ['ampere', 'A', 'electron flow'],
            ['kelvin', 'K', 'thermodynamic temperature'],
            ['mole', 'mol', 'amount of substance'],
            ['candela', 'cd', 'luminousintensity']
        ],dtype=str)
    derived = pd.DataFrame(
        columns=['Unit', 'Symbol', 'Physical quantities'],
        data=[
            ['hertz', 'Hz', 'frequency'],
            ['newton', 'N', 'force'],
            ['pascal', 'Pa', 'pressure'],
            ['joule', 'J', 'energy'],
            ['watt', 'W', 'power'],
            ['coulomb', 'C', 'electric charge'],
            ['volt', 'V', 'electric potential'],
            ['ohm', 'Ω', 'resistor'],
            ['henry', 'H', 'nductance'],
            ['farad', 'F', 'capacitance'],
            ['tesla', 'T', 'magnetic flux density'],
            ['siemens', 'S', 'conductance'],
            ['lux', 'lx', 'lux'],
            ['lumen', 'lm', 'luminous flux'],
            ['weber', 'Wb', 'magnetic flux'],
            ['radius', 'rad', 'plane angle'],
            ['steradian', 'sr', 'solid angle'],
            ['becquerel','Bq','disintegration rate'],
            ['gray','Gy','radiation absorbed dose'],
            ['sievert','Sv','dose equivalent'],
            ['katal','kat','catalytic activity'],
            ['degree_Celsius','℃','Celsius temperature']
        ],dtype=str)

    m,kg,s,A,K,mol,cd=symbols('m,kg,s,A,K,mol,cd')

    Hz=1/s
    N=kg*m/s**2
    Pa=N/m**2
    J=N*m
    W=J/s
    C=A*s
    V=W/A
    Ω=V/A
    S=1/Ω
    F=C/V
    T=N/(A*m)
    Wb=T*m**2
    H=V/(A*s)
    rad=symbols('rad')
    sr=symbols('sr')
    lm=cd*sr
    lx=lm/m**2
    Bq=1/s
    Gy=J/kg
    Sv=J/kg
    kat=mol/s
    zero_Celsius=constants.zero_Celsius*K

    convert_temperature=constants.convert_temperature
    physical_constants=constants.physical_constants

    derived['Derived formula']=list(map(eval,derived.Symbol[:-1])).append('n℃=n+273.15K')
    SI_prefixes={
        'Y': (1e+24,'yotta'),
        'Z': (1e+21,'zetta'),
        'E': (1e+18,'exa'),
        'P': (1e+15,'peta'),
        'T': (1e-12,'tera'),
        'G': (1e+9, 'giga'),
        'M': (1e+6, 'mega'),
        'ma':(1e+4, 'myria'),
        'k': (1e+3, 'kilo'),
        'h': (1e+2,'hecto'),
        'da':(10, 'deca'),

        'd': (0.1, 'deci'),
        'c': (1e-2, 'centi'),
        'm': (1e-3, 'milli'),
        'mo': (1e-4, 'myrio'),
        'μ': (1e-6, 'micro'),
        'n': (1e-9, 'nano'),
        'p': (1e-12, 'pico'),
        'f': (1e-15, 'femto'),
        'a': (1e-18, 'atto'),
        'z': (1e-21, 'zepto'),
        'y': (1e-24, 'yocto')
    }


    loaded={}
    def load_prefix_unit(self,symbols):
        'loading SI prefix units'
        for symbol in symbols:
            symbol=str(symbol)

            prefix=self.SI_prefixes.get(symbol[0],(0,None))
            constant=prefix[0]
            SI_unit=getattr(self,symbol[1:])

            setattr(self, symbol, constant*SI_unit)
            self.loaded[symbol]=(symbol,'{}*{}'.format(constant,symbol[1:]))

    def load_unit(self,symbol,constant=None,description=None,**kwargs):
        'loading other unit'
        symbol=str(symbol)
        if not constant:
            constant=1

        expr=str(constant)
        for i,j in kwargs.items():
            if j==1:
                expr += '*self.{}'.format(i)
            else:
                expr += '*self.{}**{}'.format(i, j)

        setattr(self, symbol, eval(expr))
        self.loaded[symbol]=(expr.replace('self.',''),description)

    def __init__(self):
        self.load_unit('g', constants.gram, kg=1, description='gram')
        self.load_unit('u', constants.atomic_mass,kg=1,description='Relative atomic mass')
        self.load_unit('eV', constants.eV, J=1,description='Electron-Volt')
        self.load_unit('e', constants.elementary_charge, C=1,description='electron charge')
        self.load_unit('au', constants.au, m=1, description='astronomical_unit')
        self.load_unit('atm', constants.atm, Pa=1, description='standard atmosphere in pascals')
        self.load_unit('cal', constants.calorie, J=1,description='calorie')
        self.load_prefix_unit(['km', 'dm','cm','mm', 'nm', 'μm'])
        self.load_prefix_unit(['ms','μs', 'ns','kHz','MHz'])

        self.load_unit('L', dm=3)
        self.load_prefix_unit(['mg','mL','kJ','kcal', 'kPa','MPa','kmol','kN'])

    def __repr__(self):
        return 'Units'

UNITS=Units()

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

#-----------------Physical law

class PhysicalLaw:
    '''
    物理定律接口
    '''
    def __init__(self):
        pass

class NewtonLaw(PhysicalLaw):

    @classmethod
    def gravity(self, mass):
        gravity = mass * self.gravity_constant

    @classmethod
    def first_law(self):
        pass

    @classmethod
    def second_law(self,force=None,mass=None,a=None):
        force=mass*a

    @classmethod
    def third_law(self, force,):
        force=-force

    @classmethod
    def motion(self, x0, y0, t0, t, speed=[0, 0]):
        X = x0 + speed[0] * (t - t0)
        Y = y0 + speed[1] * (t - t0)

    @classmethod
    def friction(self, mass, coefficient):
        '''摩擦力'''
        friction = mass * coefficient * self.gravity_constant

class ThermodynamicLaw(PhysicalLaw):

    energy_conversion_efficiency = 0.8 #能量利用率

    '''温度'''
    @classmethod
    def temperature(cls, photon):
        temperature = photon.density * photon.energy
        
    @classmethod
    def thermal_equilibrium(cls,composition):
        '''热平衡
        linalg matrix:
        SHC*mass*temp-Q=SHC*mass*t0
        Q1+Q2+Q3+...=0

        系数矩阵
        [       a             ]      b
        temp         Q1  Q2  Q3  SHC*mass*t0
        SHC1*mass1   -1   0   0  SHC1*mass1*t0
        SHC2*mass2   0   -1   0  SHC2*mass2*t0
        SHC3*mass3   0   0   -1  SHC3*mass3*t0
            0        1    1   1       0
        '''
        num=len(composition)
        a=np.zeros((num+1,num+1))
        b=np.zeros((num+1,1))

        a[:-1,1:]=np.eye(num)*(-1)
        a[num,1:]=1
        a[num,0]=0
        a[:num,0]=[i.SHC*i.mass for i in composition]
        b[:num,0]=[i.SHC*i.mass*i.degree_Celsius for i in composition]
        b[num,0]=0

        temp=linalg.solve(a,b)[0]
        return temp

    @classmethod
    def energy_conservation(cls):
        dQ=dU+pdV

    '''动量守恒'''
    @classmethod
    def momentum(self,mass,velocity):
        momentum=mass*velocity

    '''熵增加定律'''

    '''做功与能量'''
    @classmethod
    def work(self,force,dist):
        energy=force*dist

    '''动能'''
    @classmethod
    def kinetic_energy(cls,mass,speed):
        return 0.5*mass*speed**2

    '''重力势能'''
    @classmethod
    def gravitational_potential_energy(cls,mass,height=None,M=None,r=None,on_earth=True):
        if on_earth:
            return mass*constants.g*height
        elif r>0:
            G=constants.gravity_constant
            return -G*M*mass/r

# -------------------------------Pure Substance
class PureSubstance:

    def __init__(self,molecule,mass,volume=None,temp=15):
        '''
         PureSubstance: Compound and Elementary Substance
         molecule:  Molecule
         temp: ℃
         amount of substance: mol
        '''
        self.molecule=molecule
        self.name=self.molecule.name
        self.anti=self.molecule.anti
        self.molecular_formula=self.molecule.formula
        self.relative_molecular_mass=self.molecule.mass

        self.set_mass_charge(mass)
        self.volume=volume  # cubic meter
        self.density=self.mass/self.volume

        self.standard_density=self.physical_property('standard density') #标准密度
        self.standard_concentration=self.get_standard_concentration()
        self.melting_point=self.physical_property('melting_point')  #熔点
        self.boiling_point=self.physical_property('boiling_point')  #沸点
        self.solubility=self.physical_property('solubility')  #溶解度(g/100g H2O)
        self.SHC=self.physical_property('specific_heat_capacity kJ/(mol*K)')  #比热容
        self.ionization=self.physical_property('ionization') #离子集

        self.set_temp(temp)
        self.state=self.get_state()
        self.pressure=self.get_pressure() #Pa

    def set_temp(self,temp):
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp
        
    def get_amount_of_substance(self,mass):
        return mass/(self.relative_molecular_mass*constants.atomic_mass*constants.Avogadro)  #mol

    def set_mass_charge(self,mass):
        self.mass=mass  #kg
        self.amount_of_substance=self.get_amount_of_substance(mass)
        self.charge=self.molecule.charge*self.amount_of_substance*constants.e*constants.Avogadro   # C

    def physical_property(self,property):
        'get physical properties'
        return MOLECULE_TABLE[property].get(self.molecular_formula,None)

    def get_state(self):
        '''There are four states of matter: solid,liquid and gas'''
        temp=self.degree_Celsius
        if temp<=self.melting_point:
            state='solid'
        elif temp<=self.boiling_point:
            state='liquid'
        elif temp>self.boiling_point:
            state='gas'

        return state

    def heat_transfer(self,heat):
        'heat:kJ'
        SHC=self.SHC
        delta_t=heat/(SHC*self.mass)
        self.set_temp(self.degree_Celsius+delta_t)
        self.state=self.get_state()  # phase transition(相变)

    def materia_exchange(self,mass=None,mol=None):
        mass=mass if mass else mol*self.relative_molecular_mass*1e-3
        self.set_mass_charge(mass)

    def solubility(self,solvent):
        '溶解'
        s=self.solubility/100*solvent.mass
        mass=self.mass if s>=self.mass else s
        return mass #kg

    def diffusion(self,volume):
        if self.state=='gas':
            self.volume=volume
            self.density=self.mass/self.volume
            self.pressure=self.get_pressure()

    def get_pressure(self):
        if self.volume is not None and self.state=='gas':
            return self.amount_of_substance*constants.gas_constant*self.degree_Kelvin/self.volume

    def annihilate(self,others):
        '湮灭'
        if (not self.anti)==others.anti:
            Δm=math.abs(self.mass-others.mass)
            E=2*Δm*constants.c**2
            return 2*Δm,E
        
    def get_standard_concentration(self):
        mass=self.standard_density/1000  #kg/dm3
        mol=self.get_amount_of_substance(mass)
        return mol #mol/L

    def __repr__(self):
        '显示单质或化合物'
        name='Elementary Substance' if len(self.molecule.atoms)==1 else 'Compound'
        return '{}({},mass={},volume={},{}℃)'.\
            format(name,self.molecular_formula,self.mass,self.volume,self.degree_Celsius)

    def __add__(self,other):
        if type(self.molecule)==Molecule and self.molecular_formula==other.molecular_formula:
            temp = ThermodynamicLaw.thermal_equilibrium([self,other])
            return PureSubstance(self.molecule,self.mass+other.mass,self.volume+other.volume,temp)


unitPURE=MOLECULES.map(lambda x:PureSubstance(x,mass=1,volume=1/x.standard_density))

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
        equation=self.equation
        def bolds_energy(molecule_dict):
            bolds_energy=0.0
            for i in molecule_dict:
                bolds_energy+=MOLECULES[i].bonds_energy*molecule_dict[i]
            return bolds_energy #kJ/mol

        return bolds_energy(self.reactant)-bolds_energy(self.product)

    def rate_equation(self,Temp,**concentration):
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

#-----------------------------------Matter
class Mixture:
    def __init__(self,volume,composition):
        '''
        composition: tuple or list of PureSubstance
        '''
        self.composition={i.molecular_formula:i for i in composition}
        self.set_mass_charge()
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.items()}
        self.volume=volume
        self.density=self.mass/self.volume
        
        self.diffusion()
        self.__temp_init()


    def __temp_init(self):
        temp=ThermodynamicLaw.thermal_equilibrium(self.composition.values())
        self.set_temp(temp)
        for i in self.composition.values():
            i.set_temp(temp)

    def set_temp(self,temp):
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp

    def set_mass_charge(self):
        mass=charge=0
        for i in self.composition.values():
            mass+=i.mass
            charge+=i.charge
        self.mass=mass
        self.charge=charge

    def solubility(self):
        '溶解'
        solvent=self.composition.get('H2O',None)
        solute={}
        if solvent is not None:
            for i,j in self.composition.items():
                solute[i]=j.solubility(solvent)
        return solvent,solute  #溶剂和溶质

    def diffusion(self):
        '扩散'
        used_volume=[i.volume for i in self.composition.values() if i.state!='gas']
        used_volume=sum(used_volume)
        residual_volume=self.volume-used_volume
        for i in self.composition.values():
            i.diffusion(residual_volume)

    def get_ions(self):
        ions={}
        if self.composition.get('H2O',None) is not None:
            for i in self.composition.values():
                ions.update(i.ionization)
        return ions

    def get_SHC(self):
        '''
        精简计算混合比热容 kJ/(mol*K)
        Q=CmΔt
        '''
        Q=0
        for i in self.composition.values():
            Q+=i.SHC*i.mass*1
        return Q/self.mass

    def heat_transfer(self,heat,SHC=None):
        'heat:kJ'
        SHC=self.get_SHC() if SHC is None else SHC
        delta_t=heat/(SHC*self.mass)
        self.set_temp(self.degree_Celsius+delta_t)
        for i in self.composition.values():
            i.set_temp(i.degree_Celsius+delta_t)
            i.state=i.get_state()

    def materia_exchange(self,PureSubstance):
        p=PureSubstance
        p=pd.Series(p,index=[p.molecular_formula])
        p0=pd.Series(self.composition)
        p1=p0+p
        p1.fillna(0)
        self.composition=dict(p1)
        self.set_mass_charge()


    def __repr__(self):
        return "Mixture: mass={},volume={},{}℃".\
            format(self.mass,self.volume,self.degree_Celsius)

    def physical_law(self,law):
        '物理定律接口'
        pass

class Matter(Mixture):
    def __reaction_rate(self,T,**env):
        rate={}
        power_total=0
        for k in CHEMICAL_REACTION:
            reaction=k
            reactant,product=reaction.reactant,reaction.product
            v=reaction.rate_equation(T,**env)
            power=reaction.Qp*v
            power_total+=power
            for i,j in reacant.items():
                rate[i]=rate.get(i,0)-v*j
            for i,j in product.items():
                rate[i]=rate.get(i,0)+v*j
        return rate,power_total

    def chemical_reaction_process(self,Δt,accuracy=1e-3):

        dt=accuracy*Δt
        energy=0

        for m in np.arange(0,Δt,dt):
            T=self.Temp
            env={i:j.amount_of_substance/self.volume for i,j in self.composition.items()}

            rate,power=__reaction_rate(T,env)
            t=dt if Δt-m>=dt else Δt-m
            Q=power*t
            energy+=Q

            p=[PureSubstance(i,MOLECULES[i].mass*j*1e-3) for i,j in rate.items()]
        return p,energy
