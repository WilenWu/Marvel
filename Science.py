#-*- coding:utf8 -*-

import pandas as pd
import numpy as np
import re

from scipy import constants,linalg
from sympy import symbols
from math import pi,e,exp,inf

# System International(SI)
class units:
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
        self.load_unit('c', constants.c, m=1,s=-1,description='speed_of_light')
        self.load_unit('au', constants.au, m=1, description='astronomical_unit')
        self.load_unit('atm', constants.atm, Pa=1, description='standard atmosphere in pascals')
        self.load_prefix_unit(['km', 'dm','cm','mm', 'nm', 'μm'])

        self.load_unit('L', dm=3)
        self.load_unit('cal', constants.calorie, J=1)
        self.load_prefix_unit(['kJ', 'kcal', 'kPa','kmol'])


UNITS=units()

#-----------------subatomic particles

class Electron:
    def __init__(self,anti=False):
        self.charge=1 if anti else -1
        self.mass=constants.m_e/constants.atomic_mass
    def __str__(self):
        return 'Electron'

class Proton:
    def __init__(self, anti=False):
        self.charge = -1 if anti else 1
        self.mass=constants.m_p/constants.atomic_mass
    def __str__(self):
        return 'Proton'

class Neutron:
    def __init__(self, anti=False):
        self.charge=0
        self.mass=constants.m_u/constants.atomic_mass
    def __str__(self):
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

    def __str__(self):
        return 'Photon'

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

    def nuclear_reaction():
        'fission and fusion(裂变和聚变)'
        pass

    def __str__(self):
        return 'Atom({p},{n})'.format(p=self.protons,n=self.neutrons)


#----------------periodic_table
PERIODIC_TABLE=pd.read_excel('Marvel/Science.xlsx',sheet_name='periodic_table')
PERIODIC_TABLE.set_index(keys='element', inplace=True)
ATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron']),axis=1)
antiATOMS=PERIODIC_TABLE.apply(lambda x:Atom(x['proton'],x['neutron'],anti=True),axis=1)

#--------------------molecule
class Molecule:
    def __init__(self,formula,ionic_electron=0,anti=False):
        '''
        :param formula: int
        :param ionic_electron: int
        '''
        self.formula=formula
        self.name=MOLECULE_TABLE.loc[self.formula, 'name']
        self.atoms=self.get_atoms(self.formula)
        self.anti=anti
        self.mass,self.charge=self.get_mass_charge(ionic_electron)
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy()  #kJ/mol

    @staticmethod
    def get_atoms(formula):
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
        for i in bonds:
            bonds_energy+=CHEMICAL_BOND.loc[i,'energy(KJ/mol)']*bonds[i]
        return bonds,bonds_energy

    def __str__(self):
        return 'Molecule({})'.format(self.formula)

# chemical bond: ionic bond,covalent bond and metallic bond
CHEMICAL_BOND=pd.read_excel('Marvel/Science.xlsx',sheet_name='chemical_bond')
def __reverse(x):
    if x in ['N-H..O','N..H-O']:
        return None
    else:
        x=re.split('[-=3]',x)
        x.reverse()
        return '-'.join(x)
__CB=CHEMICAL_BOND.copy()
__CB['bond']=__CB['bond'].map(__reverse)
CHEMICAL_BOND=CHEMICAL_BOND.append(__CB.dropna()).drop_duplicates()
CHEMICAL_BOND.set_index(keys='bond', inplace=True)

MOLECULE_TABLE=pd.read_excel('Marvel/Science.xlsx',sheet_name='molecule')
MOLECULE_TABLE.loc[:,['bond','ionization']]=MOLECULE_TABLE.loc[:,['bond','ionization']].applymap(eval)
MOLECULE_TABLE.replace('Inf', inf, inplace=True)
MOLECULE_TABLE.set_index(keys='formula', inplace=True)

MOLECULES=MOLECULE_TABLE.index.to_series().map(Molecule)
antiMOLECULES=MOLECULE_TABLE.index.to_series().map(lambda x:Molecule(x,anti=True))

# -------------------------------Pure Substance
class PureSubstance:
    '''
     PureSubstance: Compound and Elementary Substance
     molecule:  molecule
     temp: ℃
     amount of substance: mol
    '''
    def __init__(self,molecule,mass,volume=None,temp=15):

        self.molecule=molecule
        self.name=self.molecule.name
        self.anti=self.molecule.anti
        self.molecular_formula=self.molecule.formula
        self.relative_molecular_mass=self.molecule.mass

        self.set_mass_charge(mass)
        self.volume=volume  # cubic meter
        self.density=self.mass/self.volume
        self.pressure=self.get_pressure() #Pa

        self.set_temp(temp)
        self.state=self.get_state()
        self.melting_point=self.physical_property('melting_point')  #熔点
        self.boiling_point=self.physical_property('boiling_point')  #沸点
        self.solubility=self.physical_property('solubility')  #溶解度(g/100g H2O)
        self.SHC=self.physical_property('specific_heat_capacity kJ/(mol*K)')  #比热容
        self.ionization=self.physical_property('ionization') #离子集

    def set_temp(self,temp):
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp

    def set_mass_charge(self,mass):
        self.mass=mass  #kg
        self.amount_of_substance=self.mass*1000/self.relative_molecular_mass  #mol
        self.charge=self.self.molecule.charge*self.amount_of_substance*constants.e   # C

    def physical_property(self,property):
        'get physical properties'
        return MOLECULE_TABLE.loc[self.molecular_formula, property]

    def get_state(self):
        '''There are four states of matter: solid,liquid and gas'''
        temp=self.degree_Celsius
        if temp<=self.melting_point:
            state='solid'
        elif temp<=self.boiling_point:
            state='liquid'
        elif temp>self.boiling_point:
            state='gas'
        else:
            state=None
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

    def __str__(self):
        '显示单质或化合物'
        name='Elementary Substance' if len(self.molecule.atoms)==1 else 'Compound'
        return '{}:{}'.format(name,self.molecular_formula)

    def __add__(self,other):
        if type(self)==type(other) and self.molecular_formula==other.molecular_formula:
            return PureSubstance(self.molecule,self.mass+other.mass,self.volume+other.volume)


#----------------Chemical_reaction
class ChemicalReaction:
    def __init__(self,equation):
        self.equation=equation
        self.reactant,self.product=self.get_stoichiometric_number()
        self.A,self.Ea=CHEMICAL_REACTION.loc[equation,['A','Ea(kJ/mol)']] #kJ/mol
        self.enzyme=CHEMICAL_REACTION.loc[equation,'enzyme']
        self.Km=CHEMICAL_REACTION.loc[equation,'Km(mmol/L)']
        self.Kcat = CHEMICAL_REACTION.loc[equation, 'Kcat(1/s)']
        self.catalytic_efficiency=self.Kcat/self.Km

        self.reaction_enthalpies=self.reaction_enthalpies() #kJ/mol
        self.ΔrHmΘ=self.reaction_enthalpies #kJ/mol
        self.reaction_heat=self.reaction_enthalpies #kJ/mol
        self.Qp=self.reaction_heat #kJ/mol

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
        def bold_energy(molecule_dict):
            bold_energy=0.0
            for i in molecule_dict:
                bold_energy+=MOLECULES[i].bond_energy*molecule_dict[i]
            return bold_energy #kJ/mol

        return bold_energy(self.reactant)-bold_energy(self.product)

    def rate_equation(self,Temp,**concentration):
        '''
        molecules,atoms or ions : float(mol/L), cover catalyst(催化剂),enzyme(酶),etc
        '''

        if self.enzyme is not None:
            c=concentration #mol/L
            enzyme=self.enzyme
            relative_conc=[c.get(i,0)/j for i,j in self.reactant.items()] #相对浓度
            s=min(relative_conc) #底物浓度(mol/L)
            e=c.get(enzyme['enzyme'],0) #酶浓度(mol/L)

            Km=self.Km # mmol/L
            Kcat=self.catalytic_efficiency*Km  #1/s
            Km=Km*1000  #mol/L
            Vmax=Kcat*e # mol/(L*s)

            v=Vmax*s/(Km+s)
            return v
        else:
            R = constants.gas_constant * 0.001 # kJ/(mol*K)
            c = {i: j * 0.001 for i, j in concentration.items()} # mmol/L
            A,Ea=self.A,self.Ea #kJ/mol
            k = A*np.exp(-Ea/(R*Temp))
            tmp=1
            for i,j in self.reactant.iterms():
                tmp*=c.get(i,0)**j
            v=k*tmp # mmol/(L*s)
            return v*1000 # mol/(L*s)

    def __str__(self):
        return self.equation



#list of reactions
CHEMICAL_REACTION=pd.read_excel('Marvel/Science.xlsx',sheet_name='chemical_equation')
CHEMICAL_REACTION.loc[:,'enzyme']=CHEMICAL_REACTION.loc[:,'enzyme'].map(eval)
CHEMICAL_REACTION.set_index(keys='equation', inplace=True)

CHEMICAL_REACTION=CHEMICAL_REACTION.index.to_series().map(Chemical_reaction)

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
    '''
    energy_conversion_efficiency(能量利用率) = 80 %
    '''
    energy_conversion_efficiency = 0.8

    '''温度'''
    @classmethod
    def temperature(self, photon):
        temperature = photon.density * photon.energy

    @classmethod
    def energy_conservation(self):
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


#-----------------------------------Mixture
class Matter:
    def __init__(self,volume,*composition):
        '''
        composition: PureSubstance
        '''
        self.composition={i.molecular_formula:i for i in composition}
        self.set_mass_charge()
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.iterms()}
        self.diffusion()
        self.volume=volume

    def __temp_init(self):
        '''
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
        num=len(self.composition)
        a=np.zeros((num+1,num+1))
        b=np.zeros((num+1,1))

        a[:-1,1:]=np.eye(num)*(-1)
        a[num,1:]=1
        a[num,0]=0
        a[:num,0]=[i.SHC*i.mass for i in self.composition.values]
        b[:num]=[i.SHC*i.mass*i.degree_Celsius for i in self.composition.values]
        b[num]=0

        temp=linalg.solve(a,b)[0]
        self.set_temp(temp)
        for i in self.composition.values:
            i.set_temp(temp)

    def set_temp(self,temp):
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp

    def set_mass_charge(self):
        mass=charge=0
        for i in self.composition.values:
            mass+=i.mass
            charge+=i.charge
        self.mass=mass
        self.charge=charge

    def solubility(self):
        '溶解'
        solvent=self.composition.get('H2O',None)
        solute={}
        if solvent is not None:
            for i,j in self.composition.iterms():
                solute[i]=j.solubility(solvent)
        return solvent,solute

    def diffusion(self):
        '扩散'
        used_volume=[i.volume for i in self.composition.values if i.state!='gas']
        used_volume=sum(used_volume)
        residual_volume=self.volume-used_volume
        for i in self.composition.values:
            i.diffusion(residual_volume)

    def get_ions(self):
        ions={}
        if self.composition.get('H2O',None) is not None:
            for i in self.composition.values:
                ions.update(i.ionization)
        return ions

    def get_SHC(self):
        '''
        精简计算混合比热容 kJ/(mol*K)
        Q=CmΔt
        '''
        Q=0
        for i in self.composition.values:
            Q+=i.SHC*i.mass*1
        return Q/self.mass

    def heat_transfer(self,heat):
        'heat:kJ'
        SHC=self.get_SHC()
        delta_t=heat/(SHC*self.mass)
        self.set_temp(self.degree_Celsius+delta_t)
        for i in self.composition.values:
            i.set_temp(i.degree_Celsius+delta_t)
            i.state=i.get_state()

    def materia_exchange(self,PureSubstance):
        m={i.molecular_formula:i for i in PureSubstance}
        m=pd.Series(m)
        m1=pd.Series(self.composition)
        n=m+m1
        n.fillna(0)
        self.composition=dict(n)
        self.set_mass_charge(mass)


    def __reaction_rate(self,T,env):
        rate={}
        power_total=0
        for k in CHEMICAL_REACTION:
            reaction=k
            reactant,product=reaction.reactant,reaction.product
            v=reaction.rate_equation(T,env)
            power=reaction.Qp*v
            power_total+=power
            for i,j in reacant.iterms():
                rate[i]=rate.get(i,0)-v*j
            for i,j in product.iterms():
                rate[i]=rate.get(i,0)+v*j
        return rate,power_total

    def chemical_reaction_process(self,Δt,accuracy=1e-3):

        dt=accuracy*Δt
        energy=0

        for m in np.arange(0,Δt,dt):
            T=self.Temp
            env={i:j.amount_of_substance/self.volume for i,j in self.composition.iterms()}

            rate,power=__reaction_rate(T,env)
            t=dt if Δt-m>=dt else Δt-m
            Q=power*t
            energy+=Q

            p=[PureSubstance(i,MOLECULES[i].mass*j*1e-3) for i,j in rate.iterms()]
        return p,energy

    def physical_law(self,law):
        '物理定律接口'
        pass



THEORY={}
THEORY['chemistry']=r'''
化学反应的本质是旧化学键断裂和新化学键形成的过程
chemical equation:
s: solid(固体), l: liquid(液体), g: gas(气体), aq: Aqueous solution(溶液)

反应是否进行由体系的吉布斯自由能(Gibbs free energy)变化确定
反应热(reaction heat)为吸热或放热反应，完全由体系的反应焓变(reaction enthalpies)决定
键能(Bond Energy)是从能量因素衡量化学键强弱的物理量。

ΔG(Gibbs free energy change) = ΔH - TΔS
Qp(reaction heat) = ΔU + pΔV
Qp = ΔH
ΔH(reaction enthalpies) = ΣΕ(reactant) — ΣΕ(product) (利用键能计算反应热)

ΔS: entropy change
T: Kelvin temperature, E: chemical bold energy

ΔG < 0 : 自发反应
ΔG = 0 : 不能反应
ΔG > 0 : 逆反应自发进行

ΔH < 0 : exothermic reaction(放热反应)
ΔH > 0 : endothermic reaction(吸热反应)

化学计量数(stoichiometric number)：化学反应方程式中,参与反应的物质前的系数,称化学计量数
化学反应速率定义为单位时间内反应物或生成物浓度的变化量

-----------rR=nM{m+}+mN{n-}
离子积常数是化学平衡常数的一种形式，多用于纯液体和难溶电解质的电离
K=[M{m+}]**n*[N{n-}]**m
Kw=[H{+}]*[OH{-}]=1e-14 mol/L
pH=-log[H+]
K: ion product(离子积), Kw: 水的离子积

Ka=[M{m+}]**n*[N{n-}]**m/[R]**r
pKa=-logKa
Ka: ionization constant(电离平衡常数)

----------- aA+bB=cC+dD
国际单位制建议反应速率
v = -1/a*d[A]/dt = -1/b*d[B]/dt = 1/c*d[C]/dt = 1/d*d[D]/dt
化学反应速率方程(the rate law or rate equation for a chemical reaction)
v = k*[A]**m*[B]**n
Arrhenius equation： k = A*exp(-Ea/(RT))
k: 反应速率常数(reaction rate constant), [X]: 表示物质X的浓度, m,n: 反应级数(order of reaction)
A: 指前因子(Pre-exponential factor), Ea: 活化能(Activation energy)
R: 摩尔气体常数(Molar gas constant), T: 开尔文温度(kelvin temperature)

Enzyme catalysis(酶促反应)
Michaelis-Menten equation：
v = Vmax[S]/(Km+[S])
Vmax=Kcat[E]
[S]: substrate concentration
[E]: enzyme concentration
Kcat: turnover number(转化数)
Km: Michaelis constant(米氏常数)
Soil enzyme activity was mostly affected by substrate concentration

|   -- reaction process --
|      /          \
|  Ea(forward)     \
|    /           Ea(reverse)
|-------------|      \
| A+B        ΔU       \
|           --|----------
|                    C+D
energy

'''
THEORY['physics']=r'''
physical constants:

'''
