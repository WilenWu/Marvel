#-*- coding:utf8 -*-

import numpy as np
from scipy import constants
from scipy.integrate import odeint
from Marvel.Physics import *
from Marvel.Chemistry import *

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
        else:
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
        return '{}({}({}),{:.4f}kg,{:.4f}m3,{:.2f}℃)'\
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
        self.set_residual_volume()
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

    def set_residual_volume(self):
        used_volume=[i.volume for i in self.composition.values() if i.state!='gas']
        used_volume=sum(used_volume)
        
        if self.volume-used_volume<=0:
            self.volume=used_volume+self.volume*0.01
        self.residual_volume=self.volume-used_volume        
        
    def diffusion(self):
        '扩散'
        residual_volume=self.residual_volume
        for i in self.composition.values():
            i.diffusion(residual_volume)
            
    def get_gas_pressure(self):
        gas_list=[i.amount_of_substance for i in self.composition.values() if i.state=='gas']
        residual_volume=self.residual_volume
        if len(gas_list)>0 and residual_volume>0:
            gas_mol=sum(gas_list)
            return gas_mol*constants.gas_constant*self.degree_Kelvin/residual_volume
        
    def get_ions(self):
        ions=set()
        if self.composition.get('H2O',None) is not None:
            for i in self.composition.values():
                if i.ionization is not None:
                    ions.update(i.ionization)
        return ions

    def get_SHC(self):
        '''
        精简计算混合比热容 kJ/(kg*K)
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
        
        self.set_residual_volume()
        self.property_update()
        self.__temp_init()


    def __repr__(self):
        return "Mixture({:.4f}kg,{:.4f}m3,{:.2f}℃)".\
            format(self.mass,self.volume,self.degree_Celsius)

    def physical_law(self,law):
        '物理定律接口'
        pass

class Matter(Mixture):
    
    def _reaction_rate(self,Y0,t,keys,reactions):
        '''for Differential Equations
        env: Q,T,*X
        '''
        Q,T,*X=Y0
        X=dict(zip(keys,X))
        Q_T={i:mol2kg(j,MOLECULES[i].mass) for i,j in X.items()} #kg/L
        Q_T=[j*MOLECULES[i].SHC for i,j in Q_T.items()] #kJ/(L*T) 每L物质上升1K所需的热量kJ
        Q_T=sum(Q_T)
        
        dQ=0   #kJ/(L*s)
        dX={}  #mol/(L*s)
        for i in reactions:
            v=i.rate_equation(T,X) # mol/(L*s)
            dQ+=-i.Qp*v # kJ/(L*s)
            dX[i.equation]={m:n*v for m,n in dict(i.reactant,**i.product).items()}
        dT=dQ/Q_T #T/s
        
        dX=pd.DataFrame(dX)
        dX.fillna(0)
        dX=dX.sum(axis=1)
        return np.concatenate(([dQ],[dT],dX.values))

    def _reaction_process(self,Δt):
        '''
        Δt(sec): tuple, list, array..
        
        Linear Differential Equations(线性微分方程组):
        Reaction:R1,R2...
        d[Z]/dt=Σzv
        ρ=m/v
        
        dQ/dt=(Qp1*v1+Qp2*v2+...)  kJ/(L*s)
        dT/dt=1/(1000Cρ)*dQ/dt=1/(1000Cρ)*(Qp1*v1+Qp2*v2+...) K/s
        
        d[A]/dt=(a1*v1+a2*v2+...)
        d[B]/dt=(b1*v1+b2*v2+...)
        d[C]/dt=...       
              
        init: Q=0,T=Temp,[reactant],[product]=0
        '''
        reactants=CHEMICAL_REACTION.map(lambda x:set(x.reactant))
        products=CHEMICAL_REACTION.map(lambda x:set(x.product))

        volume=self.volume
 
        reactants=CHEMICAL_REACTION.map(lambda x:set(x.reactant))
        products=CHEMICAL_REACTION.map(lambda x:set(x.product))
        
        env=set(self.composition)       
        while True:
            reactions=reactants.map(lambda x:x.issubset(env))
            for i in products[reactions]:
                env.update(i)       
            if reactions.sum()==reactants.map(lambda x:x.issubset(env)).sum():
                break
        
        reactions=CHEMICAL_REACTION[reactions.values]
        amount={i:j.amount_of_substance for i,j in self.composition.items()} # mol
        amount={i:round(j,6) for i,j in amount.items()} # 修正函数返回值微误差
        env={i:amount.get(i,0)/(volume*1000) for i in env} # mol/L
        
        Y0=(0,)+(a.Temp,)+tuple(env.values()) #初始
        Y=odeint(self._reaction_rate,Y0,Δt,args=(env.keys(),reactions))
        heat=pd.Series(Y[:,0]*volume,index=Δt,name='reaction heat') #kJ
        Temp=pd.Series(Y[:,1],index=Δt,name='Temp') #K
        comp=pd.DataFrame(Y[:,2:].T,index=env.keys(),columns=Δt) #mol/L

        return heat,Temp,comp

    def chemical_reaction(self,Δt):
        'Δt(sec): a collection of some kind'
        heat,Temp,comp=self._reaction_process(Δt)
        heat=float(heat[-1])
        temp=float(Temp[-1])-273.15
        comp=comp[Δt]

        comp_density={i:j.density for i,j in self.composition}        
        self.composition={}
        for i,j in comp.items():
            m=MOLECULES[i].copy()
            density=comp_density.get(i,m.standard_density)
            amount=j*self.volume
            volume=mol2kg(amount,m.mass)/density 
            pure=PureSubstance(m,amount=amount,volume=volume,temp=temp)
            self.composition[i]=pure
            
        self.set_temp(temp)
        self.set_residual_volume()
        self.property_update()

            

x,y,z=(i.copy() for i in unitPURE[['H2','O2','H2O']])
a=Matter(0.01,[x,y,z])
a.heat_transfer(2000)
a.materia_exchange([x,z])

print(a.composition)
t=[100]
heat,temp,comp=a._reaction_process(t)
print(temp)
print(comp)

