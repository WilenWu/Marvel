import pandas as pd
import numpy as np
from Marvel.Geometry import *

b=ATOMS[20]


n = len(b.electron_shells)
subshells = [i + str(j).join(['^{', '}']) for i, j in b.electron_subshells.items()]
subshells = '$' + ''.join(subshells) + '$'

Z = b.protons
A = b.protons + b.neutrons
inner = '-' + str(Z) if b.anti else '+' + str(Z)
formula = '$' + str(A).join(['^{', '}']) + str(Z).join(['_{', '}']) + str(b.element) + '$'


plt.figure()
plt.text(-3, -0.2, formula, bbox=dict(facecolor='white', edgecolor='white'), fontsize=25)
plt.text(-3, -4, subshells, bbox=dict(facecolor='white', edgecolor='white'), fontsize=10)
plt.text(-0.7, -0.2, inner, bbox=dict(alpha=0), fontsize=25)
draw_circle([0, 0], 1)

for i,j in b.electron_shells.items():
    draw_circle([-1, 0], 2+int(i), -np.pi / 6, np.pi / 6)
    plt.text(0.8+int(i), -0.2, str(j), bbox=dict(facecolor='white', edgecolor='white'), fontsize=25)

plt.xlim(-3, 2 + n)
plt.ylim(-1 - n, 1 + n)
plt.xticks([])
plt.yticks([])
plt.show()


draw_circle([-1, 0], 3, -np.pi / 6, np.pi / 6)
plt.text(1.8, -0.2, '2', bbox=dict(facecolor='white', edgecolor='white'), fontsize=25)
draw_circle([-1, 0], 4, -np.pi / 6, np.pi / 6)
plt.text(2.8, -0.2, '8', bbox=dict(facecolor='white', edgecolor='white'), fontsize=25)
draw_circle([-1, 0], 5, -np.pi / 6, np.pi / 6)
plt.text(3.8, -0.2, '8', bbox=dict(facecolor='white', edgecolor='white'), fontsize=25)




