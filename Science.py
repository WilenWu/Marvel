#-*- coding:utf8 -*-

import pandas as pd
import numpy as np

from scipy import constants,linalg
from sympy import symbols
from math import pi,e,exp

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

class Proton:
    def __init__(self, anti=False):
        self.charge = -1 if anti else 1
        self.mass=constants.m_p/constants.atomic_mass

class Neutron:
    def __init__(self, anti=False):
        self.charge=0
        self.mass=constants.m_u/constants.atomic_mass

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
        self.chemical_bonds,self.bonds_energy=self.get_bonds_energy()

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
        bonds_energy=0
        for i in bonds:
            bonds_energy+=CHEMICAL_BOND.loc[i,'energy(KJ/mol)']*bonds[i]

        return bonds,bonds_energy

# chemical bond: ionic bond,covalent bond and metallic bond
CHEMICAL_BOND=pd.read_excel('Marvel/Science.xlsx',sheet_name='chemical_bond')
CHEMICAL_BOND.set_index(keys='bond', inplace=True)

MOLECULE_TABLE=pd.read_excel('Marvel/Science.xlsx',sheet_name='molecule')
MOLECULE_TABLE.loc[:,'bond']=MOLECULE_TABLE.loc[:,'bond'].map(eval)
MOLECULE_TABLE.set_index(keys='formula', inplace=True)

MOLECULES=MOLECULE_TABLE.index.to_series().map(Molecule)
antiMOLECULES=MOLECULE_TABLE.index.to_series().map(lambda x:Molecule(x,anti=True))

# -------------------------------matter
class PureSubstance:
    '''
     PureSubstance: Compound and Elementary Substance
     molecule:  molecule
     temp: ℃
     amount of substance: mol
    '''
    def __init__(self,molecule,amount_of_substance,temp=15,volume=None):

        self.molecule=molecule
        self.name=self.molecule.name
        self.anti=self.molecule.anti
        self.molecular_formula=self.molecule.formula
        self.amount_of_substance=amount_of_substance

        self.relative_molecular_mass=self.molecule.mass
        self.mass=self.molecule.mass*self.amount_of_substance*1e-3   #kg
        self.charge=self.self.molecule.charge*self.amount_of_substance*constants.e   # C

        self.volume=volume
        self.density=self.mass/self.volume

        self.__set_temp(temp)
        self.state=self.get_state()
        self.melting_point=self.physical_property('melting_point')  #熔点
        self.boiling_point=self.physical_property('boiling_point')  #沸点
        self.solubility=self.physical_property('solubility')  #溶解度(g/100g H2O)
        self.SHC=self.physical_property('specific_heat_capacity')  #比热容
        self.ionization_constant=self.physical_property('ionization_constant')  #电离平衡常数

    def __set_temp(self,temp):
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp

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
        SHC=self.SHC
        delta_t=heat/(SHC*self.mass)
        self.__set_temp(self.degree_Celsius+delta_t)
        self.state=self.get_state()  # phase transition(相变)

    def materia_exchange(self,mass=None,mol=None):
        mass=mass if mass else mol*self.relative_molecular_mass*1e-3
        self.mass+=mass
        self.density=self.mass/self.volume

    def diffusion(self,volume):
        if self.state=='gas':
            self.volume=Volume
            self.density=self.mass/self.volume

ElementarySubstance=PureSubstance  #单质
Compound=PureSubstance  # 化合物

class Matter:
    def __init__(self,volume,*composition):
        self.composition={i.molecular_formula:i for i in composition}
        self.mass,self.charge=self.get_mass_charge()
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.iterms()}

        self.__temp_init()
        self.__physics_init() #溶解，电离等

        self.AqueousSolution
        self.SHC=self.get_SHC()  #混合比热容

    def __temp_init(self):
        '''
        linalg matrix:
        SHC*mass*temp-Q1=SHC*mass*t0
        Q1+Q2+Q3+...=0
        temp,Q1,Q2,Q3...,SHC*mass*t0
        '''
        num=len(self.composition)
        a=np.zeros((num+1,num+1))
        b=np.zeros((num+1,1))

        '系数矩阵'
        a[:-1,1:]=np.eye(num)*(-1)
        a[num,1:]=1
        a[num,0]=0
        a[:num,0]=[i.SHC*i.mass for i in self.composition.values]
        b[:num]=[i.SHC*i.mass*i.degree_Celsius for i in self.composition.values]
        b[num]=0

        temp=linalg.solve(a,b)[0]
        self.Temp=temp+273.15
        self.degree_Kelvin=self.Temp
        self.degree_Celsius=temp

    def __physics__init(self):
        pass

    def get_mass_charge(self):
        mass=charge=0
        for i in self.composition.values:
            mass+=i.mass
            charge+=i.charge
        return mass,charge

    def ionization(self,concentration,solvent='H2O'):
        -----------rR=nM{m+}+mN{n-}
        离子积常数是化学平衡常数的一种形式，多用于纯液体和难溶电解质的电离
        K=[M{m+}]**n*[N{n-}]**m
        Kw=[H{+}]*[OH{-}]=1e-14 mol/L
        pH=-log[H+]
        K: ion product(离子积), Kw: 水的离子积

        Ka=[M{m+}]**n*[N{n-}]**m/[R]**r
        pKa=-logKa
        Ka: ionization constant(电离平衡常数)




        def energy_exchange(self):
            pass
        def materia_exchange(self):
            pass

        def heat_transfer(self,heat):
            heat=SHC*Mass*Temp

        def diffusion(self):
            pass

        def pressure(self):
            return n*R*T/v

        def __reaction_rate(T,environment):
            c={}
            power_total=0
            for k in CHEMICAL_REACTION.index:
                reaction=CHEMICAL_REACTION.loc[k,'Chemical_reaction']
                reactant,product=reaction.reactant,reaction.product
                v=reaction.rate_equation(T,environment)
                power=reaction.Qp*v
                power_total+=power
                for i,j in reacant.iterms():
                    c[i]=c.get(i,0)-v*j
                for i,j in product.iterms():
                    c[i]=c.get(i,0)+v*j
            return c,power_total

        def reaction(reaction_rate,Δt,Accuracy=1,energy_conversion_efficiency=0.8,**environment):
            energy=0
            env=environment
            for m in range(0,dt,Δt):
                rate,power=reaction_rate(T,env)
                t=dt if Δt-m>=dt else Δt-m
                Q=power*t
                energy+=Q
                for i,j in rate.iterms():
                    env[i]=env.get(i,0)+j*t
            return env,energy

        def chemical_equilibrium(self):
            '化学平衡'
            pass

        def ionization(self):
            '电离'
            pass

        def ion_balance(self):
            '离子平衡'
            pass

        def annihilate(self,others):
            '湮灭'
            if (not self.anti)==others.anti:
                Δm=math.abs(self.mass-others.mass)
                E=2*Δm*constants.c**2
                return 2*Δm,E

        def physical_law(self,law):
            '物理定律接口'
            pass


class AqueousSolution(Matter):
    pass
#----------------Chemical_reaction
THEORY={}
THEORY['chemical']=r'''
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
k = A*exp(-Ea/(RT))
k: 反应速率常数(reaction rate constant), [X]: 表示物质X的浓度, m,n: 反应级数(order of reaction)
A: 指前因子(Pre-exponential factor), Ea: 活化能(Activation energy)
R: 摩尔气体常数(Molar gas constant), T: 开尔文温度(kelvin temperature)

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

class ChemicalReaction:
    def __init__(self,equation):
        self.equation=equation
        self.reactant,self.product=self.get_stoichiometric_number()
        self.A,self.Ea=CHEMICAL_REACTION.loc[equation,['pre-exponential factor','activation energy']]
        self.enzyme=CHEMICAL_REACTION.loc[equation,'enzyme']

        self.reaction_enthalpies=self.reaction_enthalpies()
        self.ΔrHmΘ=self.reaction_enthalpies
        self.reaction_heat=self.reaction_enthalpies
        self.Qp=self.reaction_heat

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
            return bold_energy

        return bold_energy(self.reactant)-bold_energy(self.product)

    def rate_equation(self,Temp,**concentration):
        '''
        molecules,atoms or ions : float(mol/L), cover catalyst(催化剂),enzyme(酶),etc

        aA+bB=cC+dD

        v = k*[A]**a*[B]**b
        Arrhenius equation： k = A*exp(-Ea/(RT))
        Energy Units	 	kJ	 	Molecular Units	 	Molecule
        Pressure Units	 	Pa	 	Temperature Units	 	K
        Base Volume Unit	cm

        Enzyme catalysis(酶促反应)
        Michaelis-Menten equation：
        v = Vmax[S]/(Km+[S])
        [S]: substrate concentration
        Soil enzyme activity was mostly affected by substrate concentration
        '''
        R=constants.gas_constant*0.001
        c={i:j*0.001 for i,j in concentration.items()}
        A,Ea,order=self.A,self.Ea,self.reactant
        k = A*np.exp(-Ea/(R*Temp))

        # 催化
        enzyme=self.enzyme
        if len(enzyme)>=2:
            enzyme,active=enzyme[:1]
            Ea1=Ea/active
            k1=A*exp(-Ea1/(R*Temp))

            relative_conc=[c.get(i,0)/j for i,j in self.reactant.items()] #相对浓度
            saturation=c.get(enzyme,0)/min(relative_conc)*100  #相对饱和度
            k=atan(saturation)*2/pi*k1

        v=1
        for i in self.reactant:
            v*=c.get(i,0)**order[i]

        return k*v*1000 # mol/(L*s)


#list of reactions
CHEMICAL_REACTION=pd.read_excel('Marvel/Science.xlsx',sheet_name='chemical_equation')
CHEMICAL_REACTION.loc[:,'enzyme']=CHEMICAL_REACTION.loc[:,'enzyme'].map(eval)
CHEMICAL_REACTION.set_index(keys='equation', inplace=True)

CHEMICAL_REACTION=CHEMICAL_REACTION.index.to_series().map(Chemical_reaction)

#--------------------------- matter

#-----------------Physical law
THEORY['physical']=r'''
physical constants:

'''

class PhysicalLaw:
    '''
    物理定律接口
    '''
    def __init__(self):
        pass

class NewtonLaw(Physical_law):
    def gravity(self, mass):
        gravity = mass * self.gravity_constant

    def first_law(self):
        pass

    def second_law(self,force=None,mass=None,a=None):
        force=mass*a

    def third_law(self, force,):
        force=-force

    def motion(self, x0, y0, t0, t, speed=[0, 0]):
        X = x0 + speed[0] * (t - t0)
        Y = y0 + speed[1] * (t - t0)

    def friction(self, mass, coefficient):
        '''摩擦力'''
        friction = mass * coefficient * self.gravity_constant

'''温度'''
def temperature(self, photon):
    temperature = photon.density * photon.energy


class ThermodynamicLaw(hysical_law):
    def energy_conservation(self):
        dQ=dU+pdV

    '''做功与能量'''
    def work(self,force,dist):
        energy=force*dist

    '''动量守恒'''
    def momentum(self,mass,speed):
        momentum=mass*speed
        #velocity
    '''熵增加定律'''
