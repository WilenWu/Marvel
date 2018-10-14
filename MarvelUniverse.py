#!/usr/bin/env python
#-*- coding:utf8 -*-
'''
Class: Title
method: lower
function: lower
constant: UPPER
'''
import pandas as pd
import numpy as np
import math
import random
from scipy import constants
from itertools import combinations,permutations,product,combinations_with_replacement
import matplotlib.pyplot as plt
import seaborn as sns

#---------------------Universe
class Universe:
    num=0
    def __init__(self,time,space,oblivion,death,survive):
        Universe.num+=1  #平行宇宙数量

        self.TIME=time  #初始时间[init,unit]
        self.SPACE=space['coordinate'] #初始空间
        # unit为每次更新时间和空间位移的最小单位
        self.star=space['star']
        self.planet=space['planet']

    def time_flow(self):
        self.TIME+=self.time['unit'] #时间流逝

    def cosmic_expansion(self,hubble=None):
        pass
        # 宇宙膨胀

    def oblivion(self):
        pass


# -------------------law and constant
class Universe_unit:
    angel = 2 * math.pi/2000  # 单位弧度
    time = 0.5  # 单位day
    mass = 0.5  # 单位kg
    energy = 1  # 单位kJ
    AU = 1.5e8  # 单位km

unit=Universe_unit()

class Universe_law:
    def __init__(self, gravity_constant=9.8):
        self.gravity_constant = gravity_constant

    '''引力'''
    def gravity(self, mass):
        gravity = mass * self.gravity_constant

    '''摩擦力'''
    def friction(self, mass, coefficient):
        friction = mass * coefficient * self.gravity_constant

    '''温度'''
    def temperature(self, photon):
        temperature = photon.density * photon.energy

    '''第二定律'''
    def second_law(self, force,mass):
        a=force/mass

    '''运动方程'''
    def motion(self, x0, y0, t0, t, speed=[0, 0]):
        X = x0 + speed[0] * (t - t0)
        Y = y0 + speed[1] * (t - t0)

    '''做功与能量'''
    def work(self,force,dist):
        energy=force*dist

    '''动量守恒'''
    def momentum(self,mass,speed):
        momentum=mass*speed
        #velocity
    '''熵增加定律'''

newton_law=Universe_law()


class photon:
    wavelength=np.random.rand()  #N(0.5um, 1)

# 引力场和光场
#能量守恒=动能+重力势能+化学能+光能/核能+生物能

#---------------------------star
class Star:
    def __init__(self,power,location=(0,0,unit.AU)):
        self.location=location
        self.power=power

    def solar_wind(self):
        pass

    def solar_flare(self):
        pass
#---------------------------planet
#Atom>>Molecule>>cell

class Planet:
    def __init__(self,R=6371,location=(0,0,0)):
        self.location=location
        self.radius=R  # 单位km
        self.mass=6e24 #kg

        '''球坐标信息'''
        self.surface=pd.DataFrame() #1000*2000
        self.core
        self.mantile
        self.crust #1000*2000*15
        ''''71%sea,29%land'''

        self.structure={'core':self.core,'mantile':self.mantile,'crust':self.crust}

        '''物质，空气，生物层叠加'''
        self.matter
        # 矿藏每次循环不断叠加流逝，达到动态平衡

        self.air
        self.living

    def cloud(self):
        pass
    def rain(self):
        pass

#--------------------atom

#--------------------molecule
#光合作用
#6CO2+6H2O（ 光照、酶、 叶绿体）→C6H12O6（CH2O）+6O2
#所有生物移动能量为mgh的倍数

#不同分子，吸收不同种类能量
#离子活性，氧化还原反应

#化学反应优先级，中毒
#有氧呼吸公式：C6H12O6+6H2O+6O2酶→6CO2+12H2O+38ATP
#无氧呼吸公式：C6H12O6→2C2H5OH+2CO2+2ATP
#ATP~ADP+PI+能量
#在有氧呼吸过程中，葡萄糖彻底氧化分解，1mol的葡萄糖在彻底氧分解以后，共释放出2870kJ的能量，其中有1161kJ的能量储存在ATP中，1709kJ以热能形式散失。利用率为40.45%
#ATP活性和生成速度和温度

#--------------------cell
#占用一个空间单元的cell体积和数量

#--------------matter
class matter:
    #初始矿石成分
    def __init__(self,composition):
        composition=composition
    def density(self):
        pass
    def state(self):
        '''There are three states of matter: solid,liquid and gas'''
        pass #根据分子表匹配形态

#--------------------crust #1000*2000*15
#----构造坐标层，物质层叠加
#lon=[-pi,pi)
#lat=[-pi/2,pi/2)
#min=5,mean=17,max=35
#depth~norm(-15,5)
#height~norm(0,4)
#2倍标准差95%
#n=height-depth

# xy=product(np.arange(2000),np.arange(1000))
# lon,lat=tuple(zip(*xy))

# xy=pd.MultiIndex.from_product([range(2000),range(1000)],names=['lon','lat'])

def polygon(x_min,x_max,y_min,y_max,num_side):
    #region
    radius=[0]
    dist=[]

    for i in range(num_side+1):
        r_add[0]=math.pi*2/(num_side-i)*0.5
        r_add[1]=radius_low=math.pi*2/(num_side-i)*3
        if i>0:
            r_add[0]+=radius[i-1]
            r_add[1]+=radius[i-1]

        radius[i]=random.uniform(r_add[0],r_add[1])

        d_add[0]=dist[i-1]*np.sin(math.pi/12)/np.sin(radius[i])
        d_add[1]
        dist[i]=random.uniform(0,high)



#地形地貌topography
#高原，山地，平原，丘陵，盆地

#山脉长度0-500
#宽度0-100

#地壳中最多的化学元素是氧，它占总重量的48.6%；其次是硅，占26.3%；以下是铝、铁、钙、钠、钾、镁。
#地壳上层为花岗岩层（岩浆岩），主要由硅－铝氧化物构成；下层为玄武岩层（岩浆岩），主要由硅－镁氧化物构成。
#O=1/2
#Si=1/4
def polygon(region):
    pass

#--------------------mantile
#厚度约2865公里
#地幔上层物质具有固态特征，主要由铁、镁的硅酸盐类矿物组成，由上而下，铁、镁的含量逐渐增加。

#--------------------core
#平均厚度约3400公里
#地核。地核还可分为外地核、过渡层和内地核三层，
#外地核厚度约2080公里，物质大致成液态，可流动；
#过渡层的厚度约140公里；
#内地核是一个半径为1250公里的球心，物质大概是固态的，主要由铁、镍等金属元素构成。
#地核的温度和压力都很高，估计温度在5000℃以上，压力达1.32亿千帕以上，密度为每立方厘米13克。

#-------------air
#atmosphere
#shader
#air1000km
#1.205kg/m3
#氮气，占78.1%；氧气占20.9%；氩气占0.93%；还有少量的二氧化碳0.07%
#空气阻力=(1/2)CρSV^2
#系数C=垂直平面体风阻系数大约1.0 球体风阻系数大约0.5
#F浮=ρ液gV排

#最外层电子



#-------------living
#cells体积和数量
#物种(人口)或者个体进化
#animal and plant

class living:
    sex=None

    def genetic(self):
        '''遗传和变异，正态分布参数'''
        pass

    def variation(self):
        '''进化'''
        pass

    def copy(self):
        pass


def dig(volume,matterr):
    '''dig a hole'''
    density=matter.density()
    energy=volume*(density**2)


#--------------------------------运行漫威宇宙

xyz = pd.DataFrame(sun + earth)
'Spherical Coordinates'
[lon, lat, dist]

xyz=1
sun=Star()
earth=Planet()

# 永恒，无限，湮灭，死亡，吞星
Eternity={'init':0,'unit':1}
Infinity={'coordinate':xyz,'star':sun,'planet':earth,'unit':1}
Oblivion={'init':0,'unit':1}
Death={'init':0,'unit':1}
Galactus={'init':0,'unit':1}

Marvel=Universe(Eternity,Infinity,Oblivion,Death,Galactus)



# 云类天气方法
# 气候 climate 天气 weather
# 台风 typhoon 龙卷风 tornado
# windy 刮风的 cloudy 阴天的 sunny 晴朗的 snowy 下雪的 showery 阵雨 foggy 有雾的 storm 暴风雨 sand storm沙尘暴 thunder雷冰雹 hail 雪 snow 雨 rain
# shower 阵雨
# storm 暴风雨
# thunder storm 雷beautiful day 好天气

from Sphere import Sphere

earth=Sphere(radius=2)
lon,lat,h=earth.create_sphere(octaves=10, persistence=0.5,multiplier=2,seed=57)
earth.draw_sphere(lon,lat,h,map='ellipse')

earth.distance([0,1],[0,1],radius=1)

# Creating spherical planetary terrain
# planetary surface
# terrainType
# plains
# hills
# mountains
# badlands
# river valleys



