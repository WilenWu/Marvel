#-*- coding:utf8 -*-

import pandas as pd
import numpy as np

from Marvel import Science
from scipy import constants
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
# 密码子
# gene mutation基因突变
# gene recombination基因重组
# chromosome variations染色体变异 结构和数目



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
class Biomembrance(Matter):
    '''
# 分割+控制物质流+信息交流

# biomembrance生物膜系统
# 跨膜运输：选择性透过
# 磷脂双分子层 镶嵌蛋白
# 被动运输passive transport 不耗能
# h20 co2 自由扩散free diffusion
# 葡萄糖 载体蛋白协助扩散facilitated diffusion
# 主动运输ative transport 耗能
# Na K Ca
    '''
    def __init__(self,radius,thickness=10,seed=None):
        '''
        biomembrance thickness:8-10nm
        phospholipid molecules:5e+6/um2

        脂质（磷脂）50%+蛋白质40%+糖2-10%
        transport protein:carrier protein and channel protein
        '''
        self.thickness=thickness
        self.surface_aera=4*PI*radius**2
        self.phospholipid_count=np.ceil(self.surface_aera*5e+6)
        self.carrier_protein=None
        self.channel_protein=None

    def passive_transport(self):
        '''
        facilitated diffusion or free diffusion
        '''
        pass

    def active_transport(self,ATP):
        pass

    def exocytosis(self,ATP,Ca):
        '胞吞和胞吐'
        pass

    def phagocytosis(self,ATP,Ca):
        pass

    def create(self):
        pass

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
        '细胞分化'
        pass

    def grow(self):
        pass

    def apoptosis(self):
        '程序性凋亡'
        pass

    def death(self):
        pass

    def proliferation(self):
        '''
        细胞增殖
        mitosis,amitosis and meissis(有丝分类，无丝分裂和减数分裂)
        '''
        pass

    def metabolism(self):
        '''
        新陈代谢
        anabolism and catabolism(合成代谢和分解代谢)
        '''
        pass

    def plasmodesmata(self):
        '胞间连丝'
        pass

class CellWall(Matter):
    '''
    # 细胞壁
# 纤维素+果胶
    bacterial(细菌):peptidoglycan(肽聚糖)
    plant:cellulose,hemi-cellulose and lignin
    fungus(真菌):chitin(可壳糖)
    '''

    def hold(self):
        '保护细胞，避免变形破裂'
        pass

    def potect(self):
        pass

    def create(self):
        pass

class Chloroplasts(Matter):
    '''
    diameter:5-100um
    chlorophyl a(叶绿素a):C55H72O5N4Mg
    chlorophyl b(叶绿素b):C55H70O6N4Mg
    carotenoids(类胡萝卜素):C40H56O2
    '''

    def __init__(self,radius):
        self.membrance=Biomembrance(radius)
        self.thylakoid_aera=None

    def photosynthesis(self,photo):
        '''
        I: water photolysis(水的光解)
        II:
        III:
        '''
        pass

class Vacuole:
    def __init__(self,radius):
        self.membrance=Biomembrance(radius)
        self.cellsap=None

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
#-------------------living

# ecosystem 生态系统
# self-regulating 自我调节
# 阳光。热能，水。空气，无机盐
# decomposer 分解动植物残骸

# biosphere 生物圈
# food chain 食物链
# food web 食物网

class Cyanobacteria(Cell):
    '''
    blue-green algae(蓝藻)
    '''
    pass

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