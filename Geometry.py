import numpy as np
import matplotlib.pyplot as plt


def draw_circle(center,r,a0=0,a1=2*np.pi,axes=plt,**kwargs):
    a,b=center
    theta=np.arange(a0,a1,0.01)
    x=a+r*np.cos(theta)
    y=b+r*np.sin(theta)
    axes.plot(x,y,**kwargs)

def draw_polygon(center,r,a0=0,a1=2*np.pi,n=6,axes=plt,**kwargs):
    a, b = center
    theta = np.arange(a0, a1+0.01, 2*np.pi/n)
    x=a+r*np.cos(theta)
    y=b+r*np.sin(theta)
    axes.plot(x, y, **kwargs)
    return x,y
