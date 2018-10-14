#-*- coding:utf8 -*-
import pandas as pd

def excel2df(*args,**kwargs):
    xlsx=pd.read_excel(*args,**kwargs)
    data=xlsx.apply(lambda x:"','".join(x),axis=1)
    data="['" + data + "'],"
    data.to_csv('excel2df.csv',index=False)
    return data

excel2df('Marvel\SI.xlsx',sheet_name='units',header=None)