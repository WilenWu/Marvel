# Marvel <img src="man/IronMan.jpg" width="20%" height="20%" align="right" />
*Create Marvel Universe*


## Overview

## Usage

## Theory
### geometry

1D|中文|equation|length
---|---|---|---
archimedean spiral|阿基米德螺线|`r=aθ`|`0.5a(θ*sqrt(1+θ**2)+ln(θ+sqrt(1+θ**2)))`
circular helix|圆柱螺旋线|`x=rcosθ,y=rsinθ,z=bθ`|`sqrt((2θr)**2+(bθ)**2)`


2D|中文|equation|aera|circumference
---|---|---|---|---
equilateral triangle|正三角形||`sqrt(3)/4*edge**2`|`3*edge`
rectangle|矩形||`ab`|`2(a+b)`
square|正方形|`abs(x)+abs(y)=a`<br>`edge=sqrt(2)*a`|`edge**2`|`4*edge`
rhombus|菱形|`abs(x/a)+abs(y/b)=1`<br>`edge=sqrt(a**2+b**2)`|`2ab`|`4*edge`
circle|圆形|`x**2+y**2=r**2`|`πr**2`|`2πr`
ellipse|椭圆|`x**2/a**2+y**2/b**2=1`|`πab`|`2πb+4(a-b)`
ring|圆环|`(x**2+y**2-R**2)(x**2+y**2-r**2)=0`|`π(R**2-r**2)`|`2π(R+r)`

3D|中文|equation|volume|surface area
---|---|---|---|---
regular tetrahedron|正四面体||`sqrt(2)/12*edge**3`|`sqrt(3)*edge**2`
cuboid|长方体||`abc`|`2*(ab+ac+bc)`
cube|正方体|`abs(x)+abs(y)=a (0<=z<=a)`<br>`edge=sqrt(2)*a`|`edge**3`|`6*edge**2`
regular octahedron|正八面体|`abs(x)+abs(y)+abs(z)=a`<br>`edge=sqrt(2)*a`|`(4/3)a**3`|`4*sqrt(3)*a**2`
rhombus octahedron|八面菱体|`abs(x/a)+abs(y/b)+abs(z/c)=1`|`(4/3)abc`|`4abc*sqrt(1/a**2+1/b**2+1/c**2)`
ellipsoid|椭球体|`(x/a)**2+(y/b)**2+(z/c)**2=1`|`(4/3)πabc`|`(4/3)π(ab+bc+ca)`
globe|球体|`x**2+y**2+z**2=r**2`|`(4/3)πr**3`|`4πr**2`
cone|圆锥体|`(x/r)**2+(y/r)**2=(z/h)**2 (0<=z<=h)`|`(1/3)hπr**2`|`πr(r+sqrt(r**2+h**2))`
cylinder|圆柱体|`x**2+y**2=r**2 (0<=z<=h)`|`hπr**2`|`2πr(r+h)`
torus|圆环|`(R-sqrt(x**2+y**2))**2+z**2=r**2`|`(2πR)(2πr)`|`(πr**2)(2πR)`

### chemistry

- 化学反应的本质是旧化学键断裂和新化学键形成的过程
> 物质形态(State):<br>
> s: solid(固体)<br>
> l: liquid(液体)<br>
> g: gas(气体)<br>
> aq: Aqueous solution(溶液)<br>
>
> 化学计量数(stoichiometric number)：化学反应方程式(chemical equation)中,参与反应的物质前的系数,称化学计量数

- 反应是否进行由体系的吉布斯自由能(Gibbs free energy)变化确定<br>
ΔG < 0 : 自发反应<br>
ΔG = 0 : 不能反应<br>
ΔG > 0 : 逆反应自发进行<br>
```
ΔG(Gibbs free energy change) = ΔH - TΔS
ΔS: entropy change
T: Kelvin temperature
```
- 反应热(reaction heat)为吸热或放热反应，完全由体系的反应焓变(reaction enthalpies)决定<br>
ΔH < 0 : exothermic reaction(放热反应)<br>
ΔH > 0 : endothermic reaction(吸热反应)<br>
```
Qp(reaction heat) = ΔU + pΔV
Qp = ΔH
ΔH(reaction enthalpies) = ΣΕ(reactant) — ΣΕ(product) (利用键能计算反应热)
E: chemical bold energy
```
- 键能(Bond Energy)是从能量因素衡量化学键强弱的物理量。
- 化学反应速率：为单位时间内反应物或生成物浓度的变化量
国际单位制建议
```
v = -1/a*d[A]/dt = -1/b*d[B]/dt = 1/c*d[C]/dt = 1/d*d[D]/dt
```
化学反应速率方程(the rate law or rate equation for a chemical reaction)<br>
aA+bB=cC+dD
```
v = k*[A]**m*[B]**n
Arrhenius equation： k = A*exp(-Ea/(RT))
```
> k: 反应速率常数(reaction rate constant)<br>
> [X]: 表示物质X的浓度<br>
> m,n: 反应级数(order of reaction)<br>
> A: 指前因子(Pre-exponential factor)<br>
> Ea: 活化能(Activation energy)<br>
> R: 摩尔气体常数(Molar gas constant)<br>
> T: 开尔文温度(kelvin temperature)<br>

Enzyme catalysis(酶促反应)
```
Michaelis-Menten equation：
v = Vmax[S]/(Km+[S])
Vmax=Kcat[E]
```
> [S]: substrate concentration<br>
> [E]: enzyme concentration<br>
> Kcat: turnover number(转化数)<br>
> Km: Michaelis constant(米氏常数)<br>

```
|   -- reaction process --
|      /          \
|  Ea(forward)     \
|    /           Ea(reverse)
|-------------|      \
| A+B        ΔU       \
|           --|----------
|                    C+D
energy
```
- 电离(ionization)<br>
离子积常数是化学平衡常数的一种形式，多用于纯液体和难溶电解质的电离
```
物质R电离出M^m和N^n-方程式：rR=nM^m++mN^n-
# --------------
离子积(ion product):
K=[M^m+]**n*[N^n-]**m
水的离子积:
Kw=[H+]*[OH-]=1e-14 mol/L
pH=-log[H+]
电离平衡常数(ionization constant):
Ka=[M^m+]**n*[N^n-]**m/[R]**r
pKa=-logKa
```

### physics
- diffusion of rate(扩散率):<br>
热扩散率: α=k/(ρ*Cp)<br>
k:热传导系数, ρ:密度, Cp: 等压比热容<br>
分子扩散率:扩散系数是由于分子扩散和物质浓度梯度(或扩散驱动力)的摩尔通量之间的比例常数。
通常情况下,化合物空气中的扩散系数大约为水中的一万倍左右
