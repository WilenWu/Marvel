#-*- coding:utf8 -*-

math=r'''
--------2D--------
circle(圆):x^2+y^2=radius^2
aera=π*radius**2
circumference=2*π*radius

ellipse(椭圆):x^2/a^2+y^2/b^2=1
aera=π*ab
circumference=2πb+4(a-b)

rectangle(矩形):
aera=ab
circumference=2(a+b)

square(正方形):
aera=a**2
circumference=4a

rhombus(菱形):
aera=edge*height
circumference=4*edge

--------3D--------
cuboid(长方体):
volume=length*width*height
surface area=2*(length*width+length*height+width*height)

cube(正方体):
volume=edge**3
surface area=6*edge**2

ellipsoid(椭球体):x^2/a^2+y^2/b^2+z^2/c^2=1
volume=4/3*π*abc
surface area=4/3*π*(ab+bc+ca)

globe(球体):x^2+y^2+z^2=radius^2
volume=4/3*π*radius**3
surface area=4*π*radius**2

cone(圆锥体):
volume=1/3*πr**2*h
surface area=πr*(r+sqrt(r**2+h**2))

cylinder(圆柱体):
volume=πr**2*h
surface area=2πr(r+h)

regular tetrahedron(正四面体):
volume=sqrt(2)/12*a**3
surface area=sqrt(3)*a**2
'''

chemistry=r'''
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

physics=r'''
diffusion of rate(扩散率):
热扩散率: α=k/(ρ*Cp)
k:热传导系数, ρ:密度, Cp: 等压比热容
分子扩散率:
扩散系数是由于分子扩散和物质浓度梯度(或扩散驱动力)的摩尔通量之间的比例常数
通常情况下,化合物空气中的扩散系数大约为水中的一万倍左右

'''
