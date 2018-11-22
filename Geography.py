#-*- coding:utf8 -*-

import numpy  as np
import pandas as pd
import random
from  numpy import sin,cos,arccos,arcsin,sqrt
from noise import pnoise2,pnoise3


class cube:
    pass


def timer(func):
    import datetime
    from functools import wraps

    @wraps(func)
    def decorated(*args, **kwargs):
        starttime = datetime.datetime.now()
        res=func(*args, **kwargs)
        endtime = datetime.datetime.now()
        print('time used {} sec'.format((endtime - starttime).seconds))
        return res
    return decorated

class Sphere:

    def __init__(self,radius):
        self.radius=radius

    def coord_trans(self,longitude,latitude,radius=None):
        '''
        coordinate transformation
        Cartesian coordinates
        '''
        lat,lon=latitude,longitude
        if radius is None:
            r=self.radius
        else:
            r=radius

        x = r * cos(lat) * cos(lon)
        y = r * cos(lat) * sin(lon)
        z = r * sin(lat)
        return x,y,z

    def distance(self,lon,lat,angle='degrees',radius=None):
        if angle=='degrees':
            lon,lat=np.radians([lon,lat])

        a1,a2=lat
        b1,b2=lon
        if radius is None:
            r1=r2 =self.radius
        elif isinstance(radius,(int,float)):
            r1=r2 = radius
        else:
            r1,r2,*_=radius

        'Cartesian distance'
        tmp = cos(a1) * cos(a2) * cos(b1 - b2) + sin(a1) * sin(a2)
        L = sqrt(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * tmp)

        if r1==r2:
            'spherical distance'
            S = r1 * arccos(tmp)
            return L,S
        else:
            return L,None

    @timer
    def create_sphere(self,unit=0.1,multiplier=1.0,stretch=1.0,
                      seed=None,*args,**kwargs):
        '''
        Spherical coordinate system (lon,lat,r)
        0 <= r < math.inf
        0 <= lon <= PI * 2
        -PI / 2 <= lat <= PI / 2
        '''

        lon=np.arange(-180,180,unit)
        lat=np.arange(-90,90+unit,unit)

        'coordinate transformation'
        # coord = pd.MultiIndex.from_product([lon, lat], names=['lon', 'lat'])
        # coord = pd.DataFrame(index=coord).reset_index()
        lon1,lat1=np.meshgrid(lon,lat)
        x, y, z = self.coord_trans(np.radians(lon1), np.radians(lat1), radius=stretch)

        'default arguments'
        if seed is None:
            seed=random.randint(0,256)
        else:
            seed=int(seed)

        'Define numpy ufunc(universal function)'
        ufunc_pnoise3=lambda x,y,z:pnoise3(x, y, z,base=seed,*args,**kwargs)
        self.ufunc_pnoise3=np.frompyfunc(ufunc_pnoise3, 3, 1)

        h = self.ufunc_pnoise3(x,y,z)
        # h=pd.pivot(coord.lat, coord.lon, h)

        print('seed={}'.format(seed))
        return lon,lat,h

    def draw_sphere(self,lon,lat,h,map='ellipse'):
        import matplotlib.pyplot as plt
        h1=h.copy()
        h2=h.copy()
        h1[h<0]=None
        h2[h>0]=None

        fig=plt.figure()
        if map=='ellipse':
            lon,lat=np.meshgrid(lon,lat)
            lon=sqrt(1-(lat/90)**2)*lon
            plt.contourf(lon, lat, h1, cmap='Greens_r')
            plt.contourf(lon, lat, h2, cmap='Blues_r')
        elif map=='cosine':
            lon,lat=np.meshgrid(lon,lat)
            lon=cos(np.radians(lat))*lon
            plt.contourf(lon, lat, h1, cmap='Greens_r')
            plt.contourf(lon, lat, h2, cmap='Blues_r')
        else:
            plt.contourf(lon, lat, h1, cmap='Greens_r')
            plt.contourf(lon, lat, h2, cmap='Blues_r')
        plt.show()
