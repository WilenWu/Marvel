#-*- coding:utf8 -*-

import pandas as pd
from scipy import constants
# System International(SI)
class units:
    'unit: a standard of measurement'
    from sympy import symbols
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
        ])
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
    Avogadro=constants.Avogadro  #Avogadro’s number,1/mol

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

    derived['Derived formula']=list(map(eval,derived.Symbol[:-1])).append('n℃=n+273.15K')
    SI_prefixes={
        'G': (1e+9, 'giga'),
        'M': (1e+6, 'mega'),
        'k': (1e+3, 'kilo'),
        'c': (1e+2, 'centi'),
        'm': (1e-3, 'milli'),
        'μ': (1e-6, 'micro'),
        'n': (1e-9, 'nano'),
        'p': (1e-12, 'pico')
    }

    physical_constants=constants.physical_constants
    convert_temperature=constants.convert_temperature

    loaded={}
    def load_prefix_unit(self,symbols):
        'loading SI prefix units'
        for symbol in symbols:
            symbol=str(symbol)

            constant=self.SI_prefixes[symbol[0]][0]
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
                expr+='*self.{}'.format(i)
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


UNITS=units()

#electron,proton and neutron
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

#----------------atom
class Atom:
    def __init__(self,p,n,e=None,anti=False):
        e = p if e is None else e
        self.extra_nuclear_electron=e #核外电子数

        proton=Proton(anti)
        neutron=Neutron(anti)
        electron=Electron(anti)

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



