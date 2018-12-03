#-*- coding:utf8 -*-

import pandas as pd
import numpy as np
import re
from scipy import constants,linalg
from sympy import symbols
from math import pi,e,exp,inf
from copy import deepcopy
from itertools import combinations,permutations,product
import matplotlib.pyplot as plt

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

#------------------------------- Standard Model(SM)
# fermion(费米子): quark(夸克) and lepton(轻子)
# boson(玻色字): gluon(胶子), W and Z boson,photon,graviton,higgs(希格斯)

class Fermion:
    def copy(self):
        return deepcopy(self)
    def physical_law(self,law):
        '物理定律接口'
        pass
    
class Boson:
    def copy(self):
        return deepcopy(self)
    def physical_law(self,law):
        '物理定律接口'
        pass
    
#-----------------subatomic particles

class Electron(Fermion):
    def __init__(self,anti=False):
        self.symbol='e'
        self.spin=1/2
        self.magnetic_moment=None
        self.half_life=None
        self.charge=1 if anti else -1
        self.mass=constants.m_e/constants.atomic_mass
    def __repr__(self):
        return 'Electron'
    
class Neutrino(Fermion):
    '中微子'
    def __init__(self, anti=False):
        self.symbol='v'
        self.spin=1/2
        self.magnetic_moment=-1838.28197234
        self.half_life=None
        self.charge=0
        self.mass=None  

class Proton(Fermion):
    def __init__(self, anti=False):
        self.symbol='p'
        self.spin=1/2
        self.magnetic_moment=2.7928473508
        self.half_life=32661417598.559998
        self.charge = -1 if anti else 1
        self.mass=constants.m_p/constants.atomic_mass
    def __repr__(self):
        return 'Proton'

class Neutron(Fermion):
    def __init__(self, anti=False,free=False):
        self.symbol='n'
        self.spin=1/2
        self.magnetic_moment=-1.91304273
        self.half_life=611
        self.charge=0
        self.mass=constants.m_u/constants.atomic_mass
        
    def radioactive_decay(self):
        anti=self.anti
        if self.free:
            return {Proton(anti),Electron(anti),Neutrino(-anti)}
    
    def __repr__(self):
        return 'Neutron'

class Photon(Boson):
    def __init__(self,frequency):
        self.spin = 1
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


def annihilate(self,others):
    '湮灭'
    if (not self.anti)==others.anti:
        Δm=math.abs(self.mass-others.mass)
        E=2*Δm*constants.c**2
        return 2*Δm,E

def nuclear_reaction(self):
    'fission and fusion(裂变和聚变)'
    pass

def radioactive_decay(N0,half_life,t):
    N=No*(1/2)**(t/half_life)
    return N

class Field:
    pass

class ElectricField(Field):
    pass

class PhotonField(Field):
    pass
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
        return float(temp)

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

