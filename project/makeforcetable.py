# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 19:08:27 2020

@author: user
"""

fi=open("ass.force",'w+')
fi.write("\n")
fi.write("ASS_FORCE\n")
fi.write("N 1000 R 0.003 3.0\n\n")
e_const=0.89
ra=3
rh=0.3
for i in range(1000):
    r=(i+1)/1000*5
    e,f=None,None
    if r<=rh:
        e=e_const*(1/ra/ra+r*r/rh/rh/rh/rh-2/rh/rh)
        f=e_const*2*r/rh/rh/rh/rh
    elif r<=ra:
        e=e_const*(1/ra/ra-1/r/r)
        f=e_const*2/r/r/r
    else:
        e=0
        f=0
    data=str(i+1)+" "+str(r)+" "+str(e)+" "+str(f)+"\n"
    fi.write(data)


fi.write('\n\n')

fi.write("\n")
fi.write("CAPSID_REP\n")
fi.write("N 1000 R 0.00112 1.12\n\n")
for i in range(1000):
    r=(i+1)/1000*1.12
    e=4*(r**-12-r**-6)
    f=-48*r**-13+24*r**6
    data=str(i+1)+" "+str(r)+" "+str(e)+" "+str(f)+"\n"
    fi.write(data)
    

fi.close()
    