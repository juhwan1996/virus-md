# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 23:46:49 2020

@author: user
"""
from matplotlib import pyplot as plt
f=open("tmp.rdf")
lines=f.readlines()
x=[]
y=[]
tf=False
for line in lines:
    splited=line.split(" ")
    if splited[0]=="15000":
        tf=True
    if tf and len(splited)==4:
        x.append(float(splited[1]))
        y.append(float(splited[2]))

from matplotlib import pyplot as plt
plt.plot(x,y)
plt.savefig("rdf.jpg")
plt.clf()
