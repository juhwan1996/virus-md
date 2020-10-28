import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import scipy

def simulation_init(xnum,ynum,znum,xlen,ylen,zlen,temp,dt):
    sumv=0
    sumv2=0
    n=0
    places=np.zeros(shape=(xnum*ynum*znum,3))
    vs=np.random.random_sample(size=(xnum*ynum*znum,3))-0.5
    meanv=np.mean(vs,axis=0)
    sumv2=np.sum(np.var(vs,axis=0))
    for i in range(xnum):
        for j in range(ynum):
            for k in range(znum):
                places[n][0]=i*xlen
                places[n][1]=j*ylen
                places[n][2]=k*zlen
                n+=1
    fs=np.sqrt(3*temp/sumv2)
    vs-=meanv
    vs*=fs
    print(np.sum(np.mean(vs*vs,axis=0)))
    print(places,vs)
    return places,vs

def calc_force(places,box,rc):
    n=len(places)
    forces=np.zeros(shape=(n,3))
    en=0
    rc2=rc*rc
    ecut=4*1/rc**12-4*1/rc**6
    for i in range(n-1):
        for j in range(i+1,n):
            xr=places[i]-places[j]
            xr-=box*np.rint(xr/box)
            r2=np.sum(xr*xr)
            if r2<=rc2:
                r2i=1/r2
                r6i=r2i**3
                ff=48*r2i*r6i*(r6i-0.5)
                forces[i]+=ff*xr
                forces[j]-=ff*xr
                en+=4*r6i*(r6i-1)-ecut
    meane=en/n
    return forces,meane

def intergrate(places,places_before,vels,forces,meane,delta):
    places_next=2*places-places_before+delta**2*forces
    vels_next=(places_next-places_before)/(2*delta)
    meanv=np.mean(vels_next,axis=0)
    meanv2=np.sum(np.mean(vels_next*vels_next,axis=0))
    temp=meanv2/3
    etot_mean=meane+0.5*meanv2
    return places_next,places,vels_next,temp,etot_mean

def simulate(xnum,ynum,znum,xlen,ylen,zlen,temp,dt,box,rc,tnum,collectnum):
    pl,vs=simulation_init(xnum,ynum,znum,xlen,ylen,zlen,temp,dt)
    pl_m=pl-dt*vs
    pl_collect=np.zeros(shape=(collectnum,xnum*ynum*znum,3))
    vs_collect=np.zeros(shape=(collectnum,xnum*ynum*znum,3))
    T_collect=np.zeros(shape=collectnum)
    for i in range(tnum):
        forces,meane=calc_force(pl,box,rc)
        pl,pl_m,vs,temp,etot=intergrate(pl,pl_m,vs,forces,meane,dt)
        if i==0:
            print(forces)
            print("T:",temp)
            print("E:",etot)
        if i%10==0:
            print("processing: ",i)
        if i>=tnum-collectnum:
            pl_collect[i-(tnum-collectnum)]=pl
            vs_collect[i-(tnum-collectnum)]=vs
            T_collect[i-(tnum-collectnum)]=temp
    print("T:",temp)
    print("E:",etot)
    return pl,vs,pl_collect,vs_collect,T_collect


'''
def calc_force_multiprocess(id,places,box,rc,start,end):
    n=len(places)
    forces=np.zeros(shape=(n,3))
    en=0
    rc2=rc*rc
    ecut=4*1/rc**12-4*1/rc**6
    npart=0
    for i in range(start,end):
        for j in range(i+1,n):
            xr=places[i]-places[j]
            xr-=box*np.rint(xr/box)
            r2=np.sum(xr*xr)
            npart+=1
            if r2<=rc2:
                r2i=1/r2
                r6i=r2i**3
                ff=48*r2i*r6i*(r6i-0.5)
                forces[i]+=ff*xr
                forces[j]-=ff*xr
                en+=4*r6i*(r6i-1)-ecut
    return forces,en,npart

def calc_force_total_multiprocess(places,box,rc,numofprocess):
    n=len(places)
    totalpair=(n*(n-1))/2
    k=1
    sumval=0
    processrangelst=[0]
    totalforces=np.zeros(shape=(n,3))
    result=Queue()
    totalen=0
    for i in range(n-1):
        sumval+=n-1-i
        if sumval>=totalpair*k/numofprocess:
            processrangelst.append(i)
            k+=1
    if len(processrangelst)!=numofprocess+1:
        print("error!")
        print(len(processrangelst))
    th=[]
    for i in range(numofprocess):
        th.append(Process(target=calc_force_multiprocess, args=(i,places,box,processrangelst[i],processrangelst[i+1])))
    for i in range(numofprocess):
        th[i].start
    for i in range(numofprocess):
        th[i].join
    result.put('STOP')
    while True:
        tmp=result.get()
        if tmp=='STOP':
            break
        else:
            forces,en,npart=tmp
            print(forces)
            totalforces+=forces
            totalen+=en*npart
    return totalforces,totalen/n


def simulate_multiprocess(xnum,ynum,znum,xlen,ylen,zlen,temp,dt,box,rc,tnum,collectnum,numofprocess):
    pl,pl_m,vs=simulation_init(xnum,ynum,znum,xlen,ylen,zlen,temp,dt)
    pl_collect=np.zeros(shape=(collectnum,xnum*ynum*znum,3))
    for i in range(tnum):
        forces,meane=calc_force_total_multiprocess(pl,box,rc,numofprocess)
        forces2,meane2=calc_force(pl,box,rc)
        print(np.mean(forces*forces))
        pl,pl_m,vs,temp,etot=intergrate(pl,pl_m,vs,forces,meane,dt)
        if i==0:
            print(temp)
            print(etot)
        if i%10==0:
            print("processing: ",i)
        if i>=tnum-collectnum:
            pl_collect[i-(tnum-collectnum)]=pl
    print(temp)
    print(etot)
    return pl,vs,pl_collect

'''
def gr(x_collect,nhis,box):
    collectnum=len(x_collect)
    npart=len(x_collect[0])
    delg=box/(2*nhis)
    g=np.zeros(shape=nhis)
    for num in range(collectnum):
        x=x_collect[num]
        for i in range(npart-1):
            for j in range(i+1,npart):
                xr=x[i]-x[j]
                xr-=box*np.rint(xr/box)
                r=np.sqrt(np.sum(xr*xr))
                if r<box/2:
                    ig=int(r/delg)
                    g[ig]+=2
    print(g)
    for i in range(nhis):
        rho=npart/(box**3)# ???
        r=delg*(i+0.5)
        vb=((i+1)**3-i**3)*delg**3
        nid=(4/3)*np.pi*vb*rho
        g[i]/=npart*nid*collectnum
    r=(np.array(range(nhis))+0.5)*delg
    print(r)
    print(g)
    return r,g

def vdistr(vs_collect,nhis,box):
    collectnum=len(vs_collect)
    npart=len(vs_collect[0])
    v_collect=np.sqrt(np.sum(vs_collect*vs_collect,axis=2))
    v_collect=v_collect.reshape((-1,))
    delg=np.max(v_collect)/nhis
    vdist=np.zeros(shape=nhis)
    for v in v_collect:
        ig=int(v/delg)
        if nhis==ig:
            vdist[ig-1]+=1
        else:
            vdist[ig]+=1
    v_metric=(np.array(range(nhis))+0.5)*delg
    return v_metric,vdist/collectnum

def Tequibriump(places,vs,box):
    places_revised=places-box*np.rint(places/box)
    vs1sum=0
    n1=0
    vs2sum=0
    n2=0
    for i in range(len(places)):
        if places_revised[i][0]>=0:
            n1+=1
            vs1sum+=np.sum(vs[i]*vs[i])
        else:
            n2+=1
            vs2sum+=np.sum(vs[i]*vs[i])
    if n1!=0 and n2!=0:
        return vs1sum/n1/3-vs2sum/n2/3
    return None

def Tequibriump_vec(places_collect,vs_collect,box,dt):
    diff1=[]
    xs=[]
    for i in range(len(places_collect)):
        val=Tequibriump(places_collect[i],vs_collect[i],box)
        if val:
            diff1.append(val)
            xs.append(i*dt)
    return xs,diff1



if __name__=="__main__":
    box=6.0
    temp=0.3
    dt=0.01
    nbin=100
    #pl,pl_m,vs=simulation_init(3,3,3,5,5,5,2,0.1)
    #vs=np.zeros(shape=(2,3))
    #print(np.mean(vs*vs))
    x_collect=np.zeros(shape=(20,125,3))
    vs_collect=np.zeros(shape=(20,125,3))
    '''
    for i in range(5):
        pl,vs,x_collect[i],vs_collect[i]=simulate(5,5,5,3,3,3,temp,dt,box,100000000,1000,20)
        '''
    pl,vs,x_collect,vs_collect,T_collect=simulate(5,5,5,1.2,1.2,1.2,temp,dt,box,100000000,800,200)
    #distx=x_collect[:,0,:]-x_collect[:,1,:]
    #print(distx-box*np.rint(distx/box))
    #pl,vs,x_collect=simulate_multiprocess(5,5,5,2,2,2,1,0.01,10,100000000,1000,920,4)
    vdistr(vs_collect,100,10)

    r,g=gr(x_collect,nbin,box)
    print(r,g)
    plt.title("g(r) vs r by MD simulation. (Box size 6.0, lattice size 1.2, N=125,T=6~8s)")
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.plot(r,g)
    plt.savefig("crosstesttest_1.png")
    plt.clf() # clear figure

    plt.title("t vs T by MD simulation. (Box size 6.5, lattice size 1.3, N=125,T=0~15s)")
    plt.xlabel("t")
    plt.ylabel("T")
    plt.plot(np.array(range(200))*dt,T_collect)
    #plt.savefig("crosstesttest_1_T.png")
