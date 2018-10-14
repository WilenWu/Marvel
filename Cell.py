#-*- coding:utf8 -*-
from __future__ import division
import math
import numpy as np
PI=math.pi

class Biomembrance:
    def __init__(self,radius):
        '''
        biomembrance thickness:8-10nm
        phospholipid molecules:5e+6/um2

        transport protein:carrier protein and channel protein
        '''
        self.thickness=10
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

    # 胞吞和胞吐
    def exocytosis(self,ATP,Ca):
        pass

    def phagocytosis(self,ATP,Ca):
        pass

    def create(self):
        pass

class Cell:
    '''
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

    # 细胞分化
    def differentiation(self):
        pass

    def grow(self):
        pass

    # 程序性凋亡
    def apoptosis(self):
        pass

    def death(self):
        pass

    # 细胞增殖
    def proliferation(self):
        'mitosis,amitosis and meissis(有丝分类，无丝分裂和减数分裂)'
        pass

    # 新陈代谢
    def metabolism(self):
        'anabolism and catabolism(合成代谢和分解代谢)'
        pass

    #胞间连丝
    def plasmodesmata(self):
        pass

class cell_wall:
    '''
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

class chloroplasts:
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

class vacuole:
    def __init__(self,radius):
        self.membrance=Biomembrance(radius)
        self.cellsap=None

class mitochondrion:
    '''
    diameter:0.5-1um
    length:1-2um
    count:1000-2000
    '''
    # 呼吸作用
    def respiration(self):
        pass


