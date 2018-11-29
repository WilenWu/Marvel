import numpy as np
import matplotlib.pyplot as plt


def draw_circle(o,r,a0=0,a1=2*np.pi):
    a,b=o
    theta=np.arange(a0,a1,0.01)
    x=a+r*np.cos(theta)
    y=b+r*np.sin(theta)
    plt.plot(x,y,linewidth=3)

