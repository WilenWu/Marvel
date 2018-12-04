时#-*- coding:utf8 -*-

import pandas as pd
import numpy as np

from Marvel import Science
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib import patches,collections
from Science import Molecule,ChemicalReaction
UNITS.load_unit('ATP', 30.54, kJ=1, mol=-1)

# 元素
# O Si C N H
# 地壳[48.6,26.3,0.087,0.03,0.76]
# 细胞[65,0,18,3,10]
# 大量元素macroelement=[C,H,O,N,P,S,K,Ca,Mg]
# 微量元素microelement=[Fe,Mn,Zn,Cu,B,Mo,...]
#
# # 组成细胞的化合物
# 无机
# 水85-90
# 结合水(全部水分的4.5)和自由水
#
# 无机盐1-1.5
# 阳离子Na+ K+ Ca2+ Mg2+ Fe2+ Fe3+
# 阴离子Cl- SO4^2- PO4^3- HCO3-
#
# 有机
# protein蛋白质7-10
# 脂质1-2
# 糖类和核酸1-1.5
#



# R:
# H甘氨酸
# CH3丙氨酸

#
# 脱水缩合
#
# # 核酸 nucleic acid
# DNA and RNA 细胞质
# 磷酸PO4+核糖+碱基
# ACG+T|U
# A=T,G=C
# tRNA
# transcription转录 and translation翻译
# rRNA
# mRNA
#
# 密码子codon
# gene mutation基因突变
# gene recombination基因重组
# chromosome variations染色体变异 结构和数目

#密码子
#密码子，基因矩阵

#-----------------有机分子(organic molecule)
class Carbohydrate(Molecule):
    '''糖类
    Monosaccharide(单糖,C6H12O6): glucose(葡萄糖),ribose(核糖,C5H12O5)
    Disaccharide(二糖,C12H22O11): sugar(蔗糖),maltose(麦芽糖)
    Polysaccharide(): starch(淀粉),cellulose(纤维素)
    '''
    pass

class Lipid(Molecule):
    '''脂质
    Fat(脂肪)
    Phospholipid(磷脂)
    Sterol(固醇)
    '''
    pass

class Protein(Molecule):
    '''蛋白质
# amino acid 氨基酸 20多种
# H2N-CHR-COOH
    '''
    pass

class OrganicCompounds(PureSubstance):
    pass

#-----------------生化反应和酶促反应
class BiochemicalReaction(ChemicalReaction):
    '反应物生成物可以不用化学式'
    pass


class EnzymaticReaction(BiochemicalReaction):

    def rate_equation(self,Temp,**concentration):
        '''
        molecules,atoms or ions : float(mol/L), cover enzyme(酶),etc
        '''
        c = concentration  # mol/L

        enzyme=self.enzyme
        relative_conc=[c.get(i,0)/j for i,j in self.reactant.items()] #相对浓度
        s=min(relative_conc) #底物浓度(mol/L)
        e=c.get(enzyme,0) #酶浓度(mol/L)

        Km=self.Km # mmol/L
        Kcat=self.catalytic_efficiency*Km  #1/s
        Km=Km*1000  #mol/L
        Vmax=Kcat*e # mol/(L*s)

        v=Vmax*s/(Km+s)
        return v
#--------------------细胞结构

# cellular metabolism 细胞代谢
# ADP+PI=AT
# 主要来源
# cell respiration 细胞呼吸
# aerobic and anaerobic 有氧和无氧
# 1细胞质基质 c6h12o6=2c3h4o3+4[H]+2ATP+热能
# 2 线粒体 丙酮酸2c3h4o3+6h20=6co2+20[H]+2ATP
# 3 线粒体 24[H]+6o2=12H2o+34ATP
#
# 无氧呼吸只有1和2
# #有氧呼吸公式：C6H12O6+6H2O+6O2酶→6CO2+12H2O+38ATP
# #无氧呼吸公式：C6H12O6→2C2H5OH+2CO2+2ATP
# #C6H12O6→2C3H6O3+2ATP
#
# 光合作用photosynthesis
# 叶绿素3/4和类胡萝卜素1/4
# 6CO2+6H2O（ 光照、酶、 叶绿体）→C6H12O6（CH2O）+6O2
# 水的光解：12h20=24[h]+6o2
# co2的固定：6co2+6c5=12c3
# c3的还原：12c3+24[h]+ATP=C6H12O6+6h2o+6c5
#
# [H]:NADPH
#
#
# 细胞增殖cell proliferation
# meiosis减数分裂
# fertilization受精
# 细胞分化 cell differentiation
# 细胞凋亡apoptosis
#------------cell
#------------cell
class Biomembrance(Mixture):
    '''生物膜系统
    分割+控制物质流+信息交流
    '''
    def __init__(self,radius,thickness=0.01,seed=None):
        '''
        biomembrance thickness:8-10nm
        phospholipid molecules:5e+6/μm2

        transport protein:carrier protein and channel protein
        '''
        self.thickness=thickness #μm
        self.radius=radius #μm
        self.surface_aera=4*pi*radius**2  #μm2

        self.composition=self.__init_composition()
        self.volume=self.surface_aera*self.thickness  #μm3
        self.set_mass_charge()
        self.density=self.mass/self.volume
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.items()}

        self.carrier_protein=None
        self.channel_protein=None

    def __init_composition(self):
        '''磷脂双分子层（镶嵌蛋白）
        脂质（磷脂）50%+蛋白质40%+糖2-10%
        '''
        lm=LipidMolecule()
        cm=CarbohydrateMolecule()
        pm=ProteinMolecule()

        lcnt=np.ceil(self.surface_aera*5e+6)

        lmass=lm.mass*constants.atomic_mass*pcnt
        pmass=lmass*abs(normal(0.8,0.05))
        cmass=lmass*abs(normal(0.12,0.5))

        lipid=Lipid(lm,lmass)
        carbohydrate=Carbohydrate()
        protein=Protein()

        return {'phospholipid':lipid,'protein':protein,'carbohydrate':carbohydrate}

    def passive_transport(self):
        '''被动运输
        free diffusion(自由扩散): H2O, CO2
        facilitated diffusion(协助扩散): 葡萄糖, 载体蛋白
        '''
        pass

    def active_transport(self,ATP):
        '''主动运输
        Na+, K+, Ca^2+
        '''
        pass

    def exocytosis(self,ATP,Ca):
        '胞吞和胞吐'
        pass

    def phagocytosis(self,ATP,Ca):
        '吞噬作用'
        pass

    def create(self):
        pass


class CellWall(Mixture):
    '''细胞壁
    bacteria(细菌): peptidoglycan(肽聚糖,C3H5NO4), 0.1-0.3μm
    plant: cellulose(纤维素) and lignin(果胶), 5-10μm
    fungus(真菌): chitin(可壳糖), 0.1-0.25μm
    '''
    def __init__(self,radius,thickness,species='bacteria'):
        self.thickness=thickness #μm
        self.radius=radius #μm
        self.surface_aera=4*pi*radius**2  #μm2

        self.composition=self.__init_composition(species)
        self.volume=self.surface_aera*self.thickness  #μm3
        self.set_mass_charge()
        self.density=self.mass/self.volume
        self.mass_percent={i:j.mass/self.mass for i,j in self.composition.items()}

    def __init_composition(self,species):
        if species=='bacteria':
            pass
        elif species=='fungus':
            pass
        elif species=='plant':
            pass
        else:
            return None

    def hold(self):
        '保护细胞，避免变形破裂'
        pass

    def protect(self):
        pass

    def create(self):
        pass


class Photosynthesis:
    '''光合作用
    6CO2+6H2O=C6H12O6+6O2
    '''
    def __init__(self,photo,enzyme):
        pass

    @classmethod
    def light_reaction(cls):
        '''水的光解
        24*NADP+12*H2O=24*NADPH+6*O2
        ADP+Pi=ATP+H2O
        '''
        pass

    @classmethod
    def carbon_fixation_reaction(cls):
        '''固碳反应
        CO2的固定：6*CO2+6*C5=12*C3
        C3的还原：12*C3+24*NADPH=C6H12O6+6H2O+6C5
        ATP+H2O= ADP+Pi
        '''
        pass


class Thylakoid(Mixture):
    '''类囊体
    叶绿素3/4和类胡萝卜素1/4
    类囊体中和光合成反应有直接关系的有：类胡萝卜素、质体醌、叶绿素、质体蓝素、铁氧化还原蛋白等
    类囊体膜由蛋白质与脂质所组成，其重量比约1∶1。
    脂质的成分主要为糖脂质（半乳糖脂质和硫脂质）
    蛋白主要有细胞色素b6/f复合体、质体醌（PQ）、质体蓝素（PC）、铁氧化还原蛋白、黄素蛋白、光系统Ⅰ、光系统Ⅱ复合物等。
    '''
    def __init__(self,diameter=0.25,height=0.05):
        ''' ellipsoid(椭球体)
        0.25-0.8μm
        chlorophyl a(叶绿素a):C55H72O5N4Mg
        chlorophyl b(叶绿素b):C55H70O6N4Mg
        carotenoids(类胡萝卜素):C40H56O2
        '''
        self.diameter=diameter
        self.radius=self.diameter/2
        self.height=height
        self.volume=4/3*pi*self.radius**2*self.height/2
        self.surface_aera=4/3*pi*(self.radius**2+2*self.radius*self.height/2)

        self.light_reaction=Photosynthesis.light_reaction #光反应



class Chloroplasts(Mixture):
    '''叶绿体，40-60个基粒
    diameter=5-100um
    grana(叶绿体基粒): diameter=0.15-0.2μm, length=0.5μm, 10-100个类囊体
    stroma(基质)主要成分包括：
    碳同化相关的酶类：如RuBP羧化酶占基质可溶性蛋白总量的60%。
    叶绿体DNA、蛋白质合成体系：如，ctDNA、各类RNA、核糖体等。
    一些颗粒成分：如淀粉粒、质体小球和植物铁蛋白等
    '''

    def __init__(self,diameter,height,grana_cnt=50):
        self.diameter=diameter
        self.radius=self.diameter/2
        self.height=height
        self.volume=4/3*pi*self.radius**2*self.height/2
        self.surface_aera=4/3*pi*(self.radius**2+2*self.radius*self.height/2)

        self.membrance=Biomembrance(radius)

        self.grana_cnt=grana_cnt
        self.thylakoid_cnt=self.grana_cnt*randint(10,100)
        self.thylakoid=Thylakoid()
        self.thylakoid_aera=None

        self.stroma=self.__get_stroma()

    def photosynthesis(self,photo):
        self.thylakoid.light_reaction()
        Photosynthesis.carbon_fixation_reaction(self.stroma['c'])


class Vacuole:
    def __init__(self,radius):
        self.membrance=Biomembrance(radius)
        self.cellsap=None

# ADP+PI=ATP
# 主要来源
# cell respiration 细胞呼吸
# aerobic and anaerobic 有氧和无氧
# 1细胞质基质 c6h12o6=2c3h4o3+4[H]+2ATP+热能
# 2 线粒体 丙酮酸2c3h4o3+6h20=6co2+20[H]+2ATP
# 3 线粒体 24[H]+6o2=12H2o+34ATP
#
# 无氧呼吸只有1和2
# #有氧呼吸公式：C6H12O6+6H2O+6O2酶→6CO2+12H2O+38ATP
# #无氧呼吸公式：C6H12O6→2C2H5OH+2CO2+2ATP
# #C6H12O6→2C3H6O3+2ATP

class Mitochondrion:
    '''
    diameter:0.5-1um
    length:1-2um
    count:1000-2000
    '''
    def respiration(self):
        '呼吸作用'
        pass

# 细胞核nucleus
# 染色体chromosome

class Cell:
    '''
    # 细胞器organelle
# 线粒体 叶绿体
# 内质网+核糖体->蛋白质和脂质
# 高尔基体 -> 蛋白质加工
# 溶酶体 ->消化分解
# 细胞核
# 核糖体（附着 游离）
# 液泡
#
# 细胞质
# 细胞质基质
# 水 无机盐 脂质 糖 核苷酸  氨基酸 酶
# 囊泡

    diameter:
    prokaryotic cell:1-10um
    eukaryotic cell:3-30um
    '''
    def __init__(self,radius):
        self.wall=None
        self.chloroplasts=None #叶绿体
        self.vacuole=None #液泡
        self.membrance=Biomembrance(radius)
        self.cytoplasm=None #细胞质基质
        self.mitochondrion=None #线粒体
        self.endoplasmic_reticulum=None #内质网
        self.golgi_body=None #高尔基体
        self.ribosomes=None # 核糖体
        self.centrosome=None #中心体
        self.microtubule=None #微丝及微管
        self.nucleus=None #细胞核

    def cytoplasm(self):
        pass

    def move(self):
        pass

    def differentiation(self):
        'cell differentiation(细胞分化)'
        pass

    def grow(self):
        pass

    def apoptosis(self):
        'apoptosis(程序性凋亡)'
        pass

    def death(self):
        pass

    def proliferation(self):
        '''
        cell proliferation(细胞增殖)
        mitosis,amitosis and meissis(有丝分类，无丝分裂和减数分裂)
        fertilization(受精)
        '''
        pass

    def metabolism(self):
        '''
        cellular metabolism(细胞代谢)
        anabolism and catabolism(合成代谢和分解代谢)
        '''
        pass

    def plasmodesmata(self):
        '胞间连丝'
        pass



class ProkaryoticCell(Cell):
    '''原核细胞
    '''
    pass

class EukaryoticCell(Cell):
    '''真核细胞
    '''
    pass


class Prokaryotes:
    '''原核生物
    1-10µm
    '''
    pass
class Eukaryotes:
    '''真核生物
    '''
    pass


#-------------------living

# ecosystem 生态系统
# self-regulating 自我调节
# 阳光。热能，水。空气，无机盐
# decomposer 分解动植物残骸

# biosphere 生物圈
# food chain 食物链
# food web 食物网

class  Bacteria:
    '''细菌
    '''
    pass

class Fungus:
    '''真菌
    '''
    pass

class Virus:
    '''病毒
    '''
    pass

class Plant:
    pass

class Animal:
    pass

class BlueGreenAlgae(Cell):
    '''蓝藻
    '''
    self.thylakoid=Thylakoid()

# 1e+19 kJ/d
# 1%
# 10-20%能量流动

# internal environment 内环境
# body fluid 体液
# 细胞内液体2/3 外液1/3
# hormoral regulation 激素调节

# phototropsim 向光性
# auxin生长素
#
# population density 种群密度
# birth and death rate
# sex ratio
#
#
