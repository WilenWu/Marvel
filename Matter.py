#-*- coding:utf8 -*-

import numpy as np
from scipy import constants
from Physics import *
from Chemistry import *

# -------------------------------Pure Substance
def kg2mol(kg,relative_molecular_mass):
    return kg / (relative_molecular_mass * constants.atomic_mass * constants.Avogadro)

def mol2kg(mol,relative_molecular_mass):
    return relative_molecular_mass * constants.atomic_mass * constants.Avogadro * mol

class PureSubstance(Fermion):

    def __init__(self,molecule,mass=None,amount=None,volume=None,temp=15):
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

        self.set_mass_charge(mass,amount)
        self.volume=volume  # cubic meter
        self.density=self.mass/self.volume

        self.standard_density=self.physical_property('standard density') #标准密度
        self.standard_concentration=self.get_standard_concentration()
        self.melting_point=self.physical_property('melting_point')  #熔点
        self.boiling_point=self.physical_property('boiling_point')  #沸点
        self.solubility_H2O=self.physical_property('solubility')  #溶解度(g/100g H2O)
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
        return kg2mol(mass,self.relative_molecular_mass)  #mol

    def set_mass_charge(self,mass=None,amount=None):
        if mass is not None:
            self.mass=mass
            self.amount_of_substance=kg2mol(mass,self.relative_molecular_mass)  #mol
        elif amount is not None:
            self.amount_of_substance=amount
            self.mass=mol2kg(amount,self.relative_molecular_mass) #kg
        else:
            self.mass,self.amount_of_substance=(0,0)
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

    def materia_exchange(self,other):
        if self.molecular_formula==other.molecular_formula:
            temp = ThermodynamicLaw.thermal_equilibrium([self,other])
            self.set_mass_charge(self.mass+other.mass)
            self.volume=self.volume+other.volume
            self.set_temp(temp)
            
    def solubility(self,solvent):
        '溶解'
        if solvent.molecular_formula=='H2O':
            s=self.solubility_H2O/100*solvent.mass
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

    def get_standard_concentration(self):
        mass=self.standard_density/1000  #kg/dm3
        mol=self.get_amount_of_substance(mass)
        return mol #mol/L

    def __repr__(self):
        '显示单质或化合物'
        name='Elementary Substance' if len(self.molecule.atoms)==1 else 'Compound'
        return '{}({}({}),{}kg,{:.1e}m3,{:.2f}℃)'\
        .format(name,self.molecular_formula,self.state[0],self.mass,self.volume,self.degree_Celsius)
        

unitPURE=MOLECULES.map(lambda x:PureSubstance(x.copy(),mass=1,volume=1/x.standard_density))

#-----------------------------------Matter
class Mixture(Fermion):
    def __init__(self,volume,composition):
        '''
        composition: tuple or list of PureSubstance
        '''
        self.composition={i.molecular_formula:i for i in composition}
        self.volume=volume  #混合体积
        self.property_update()
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

    def property_update(self):
        self.set_mass_charge()
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.items()}
        self.density = self.mass / self.volume
        self.diffusion()

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
        ions=set()
        if self.composition.get('H2O',None) is not None:
            for i in self.composition.values():
                if i.ionization is not None:
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
        self.diffusion()

    def materia_exchange(self,list_of_PureSubstance):
        for i in list_of_PureSubstance:
            if i.molecular_formula in self.composition.keys():
                self.composition[i.molecular_formula].materia_exchange(i)
            else:
                self.composition[i.molecular_formula]=i
        self.property_update()
        self.__temp_init()


    def __repr__(self):
        return "Mixture({}kg,{:.1e}m3,{:.2f}℃)".\
            format(self.mass,self.volume,self.degree_Celsius)

    def physical_law(self,law):
        '物理定律接口'
        pass

class Matter(Mixture):
    @classmethod
    def __reaction_rate(cls,Temp,env,volume):
        rate={}
        power_total=0
        for k in CHEMICAL_REACTION:
            reaction=k
            reactant,product=reaction.reactant,reaction.product
            v=reaction.rate_equation(Temp,env) ## mol/(L*s)
            power=reaction.Qp*v*volume*1000 #kJ/s
            power_total+=power
            for i,j in reactant.items():
                rate[i]=rate.get(i,0)-v*j*volume*1000
            for i,j in product.items():
                rate[i]=rate.get(i,0)+v*j*volume*1000
        return rate,power_total #mol/s, kJ/s

    def chemical_reaction_process(self,Δt,accuracy=1e-3):
        'Δt: sec'
        dt=accuracy*Δt
        
        for m in np.arange(0,Δt,dt):
            T=self.Temp
            env={i:j.amount_of_substance/(self.volume*1000) for i,j in self.composition.items()}
            rate,power=self.__reaction_rate(T,env,self.volume)
            t=dt if Δt-m>=dt else Δt-m
            
            p=[]
            subs={i:j*t for i,j in rate.items()}
            for i,j in subs.items():
                m=MOLECULES[i].copy()
                mass=mol2kg(j,m.mass)
                v=mass/m.standard_density
                p.append(PureSubstance(m,mass=mass,volume=v,temp=self.degree_Celsius))
                
            #self.materia_exchange(p) 
            #self.heat_transfer(power*t)

x,y,z=(i.copy() for i in unitPURE[['H2','O2','H2O']])
a=Matter(2,[x,y,z])
a.heat_transfer(100)

a.materia_exchange([x,z])

a.chemical_reaction_process(30)

