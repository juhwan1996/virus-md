# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 02:55:55 2020

@author: user
"""
#x y z increasing order
import numpy as np
def atomplace(pointvectors,divtable):
    veclst=[]
    for i in range(len(divtable)):
        zsetvec=pointvectors[0:4]+i/(len(divtable)-1)*(pointvectors[4:8]-pointvectors[0:4])
        for j in range(len(divtable[i])):
            yzsetvec=zsetvec[0:2]+j/(len(divtable[i])-1)*(zsetvec[2:4]-zsetvec[0:2])
            for k in range(divtable[i][j]):
                vec=yzsetvec[0]+k/(divtable[i][j]-1)*(yzsetvec[1]-yzsetvec[0])
                veclst.append(vec)
    veclst=np.array(veclst)
    return veclst

def autodivtable(pointvectors,divtablexy,r):
    place=[]
    for i in range(len(divtablexy)):
        place2=[]
        zsetvec=pointvectors[0:4]+i/(len(divtablexy)-1)*(pointvectors[4:8]-pointvectors[0:4])
        for j in range(divtablexy[i]):
            yzsetvec=zsetvec[0:2]+j/(divtablexy[i]-1)*(zsetvec[2:4]-zsetvec[0:2])
            place2.append(np.int(np.ceil(np.sqrt(np.sum((yzsetvec[0]-yzsetvec[1])**2))/r)+1))
        place.append(place2)
    return place

def partcoord(lenofpart,partnum,height):
    pointvectors=np.zeros(shape=(8,3))
    if partnum==0:
        pointvectors[0][0]=-lenofpart/6.0
        pointvectors[1][0]=lenofpart/2.0
        pointvectors[2][1]=lenofpart/6.0*1.618
        pointvectors[3][0]=lenofpart/3.0
        pointvectors[3][1]=lenofpart/6.0*1.618
    elif partnum==1:
        pointvectors[0][0]=-lenofpart/6.0
        pointvectors[0][1]=lenofpart/3.0*1.618
        pointvectors[1][0]=-lenofpart/2.0
        pointvectors[2][1]=lenofpart/6.0*1.618
        pointvectors[3][0]=-lenofpart/6.0
    elif partnum==2:
        pointvectors[0][0]=lenofpart/3.0
        pointvectors[0][1]=lenofpart/6.0*1.618
        pointvectors[1][1]=lenofpart/2.0*1.618
        pointvectors[2][1]=lenofpart/6.0*1.618
        pointvectors[3][0]=-lenofpart/6.0
        pointvectors[3][1]=lenofpart/3.0*1.618
    k=1-2/1.618/1.618*height
    pointvectors[4:8]=k*pointvectors[0:4]
    for i in range(4,8):
        pointvectors[i][2]=height
    return pointvectors

def univec(planum): #vec[2]=1
    vec=np.zeros(shape=3)
    const=1.618
    if planum==0:
        vec[1]=-1.0
    elif planum==1:
        vec[0]=-const*const
        vec[1]=const
        vec[2]=1
    elif planum==2:
        vec[0]=const*const
        vec[1]=const
        vec[2]=1
    elif planum==3:
        vec[0]=-const*const*3
        vec[1]=const*3
        vec[2]=1
    elif planum==4:
        vec[0]=const*const*3
        vec[1]=const*3
        vec[2]=1
    elif planum==5:
        vec[1]=const*3
        vec[2]=1
    norm = np.linalg.norm(vec)
    return vec/norm

def magplace(lenofpart,partnum,height,r):
    pointvectors= partcoord(lenofpart,partnum,height/2)[4:8] #middle
    print(pointvectors)
    attractionpoints=np.zeros(shape=(10,3))
    for i in range(4):
        lst=[1,2,4,5]
        attractionpoints[i]=pointvectors[0]+(pointvectors[1]-pointvectors[0])*(float(lst[i])/6)
        if partnum==0:
            attractionpoints[i]+=r/2*univec(0)
        elif partnum==1:
            attractionpoints[i]+=r/2*univec(1)
        elif partnum==2:
            attractionpoints[i]+=r/2*univec(2)
    for i in range(2):
        attractionpoints[i+4]=pointvectors[0]+(pointvectors[2]-pointvectors[0])*(float(i+1)/3)
        if partnum==0:
            attractionpoints[i+4]+=r/2*univec(3)
        elif partnum==1:
            attractionpoints[i+4]+=r/2*univec(4)
        elif partnum==2:
            attractionpoints[i+4]-=r/2*univec(5)
    for i in range(2):
        attractionpoints[i+6]=pointvectors[2]+(pointvectors[3]-pointvectors[2])*(float(i+1)/3)
        if partnum==0:
            attractionpoints[i+6]+=r/2*univec(5)
        elif partnum==1:
            attractionpoints[i+6]-=r/2*univec(3)
        elif partnum==2:
            attractionpoints[i+6]-=r/2*univec(4)
    for i in range(2):
        attractionpoints[i+8]=pointvectors[3]+(pointvectors[1]-pointvectors[3])*(float(i+1)/3)
        if partnum==0:
            attractionpoints[i+8]+=r/2*univec(2)
        elif partnum==1:
            attractionpoints[i+8]+=r/2*univec(0)
        elif partnum==2:
            attractionpoints[i+8]+=r/2*univec(1)
    return attractionpoints

def atomwrite(place,fileplace):
    f=open(fileplace,'w+')
    f.write("\n%d atoms\n\nCoords\n\n" % len(place))
    for i in range(len(place)):
        f.write("%d %f %f %f\n" % (i,place[i][0],place[i][1],place[i][2]))
    f.write('\n\nTypes\n\n')
    for i in range(len(place)):
        f.write("%d 1\n" % i)
    f.close()
    
partlink=[[9,9,8,7,10,11,15,14,5,6],[3,4,5,6,12,13,11,10,8,7],[4,3,2,1,14,15,13,12,2,1]]
def atomwrite2(place,bondsite,partnum,fileplace):
    f=open(fileplace,'w+')
    f.write("\n%d atoms\n\nCoords\n\n" % (len(place)+len(bondsite)))
    for i in range(len(place)):
        f.write("%d %f %f %f\n" % (i+1,place[i][0],place[i][1],place[i][2]))
    for i in range(len(bondsite)):
        f.write("%d %f %f %f\n" % (i+1+len(place),bondsite[i][0],bondsite[i][1],bondsite[i][2]))
    f.write('\n\nTypes\n\n')
    for i in range(len(place)):
        f.write("%d %d\n" % (i+1,partnum+1))
    for i in range(len(bondsite)):
        f.write("%d %d\n" % (i+1+len(place),(3+partlink[partnum][i])))
    f.close()

place1=atomplace(partcoord(5.0,0,1.0),autodivtable(partcoord(5.0,0,1.0),[3,3,2],1.16))
bondsite=magplace(5.0,0,1.0,1.16)
atomwrite2(place1,bondsite,0,"V1.molc")
place2=atomplace(partcoord(5.0,1,1.0),autodivtable(partcoord(5.0,1,1.0),[3,3,2],1.16))
bondsite=magplace(5.0,1,1.0,1.16)
atomwrite2(place2,bondsite,1,"V2.molc")
place3=atomplace(partcoord(5.0,2,1.0),autodivtable(partcoord(5.0,2,1.0),[3,3,2],1.16))
bondsite=magplace(5.0,2,1.0,1.16)
atomwrite2(place3,bondsite,2,"V3.molc")