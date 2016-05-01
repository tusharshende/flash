# -*- coding: utf-8 -*-
"""
Created on Tue Apr 07 12:39:58 2015

@author: vhd
"""

import scipy
import matplotlib.pyplot as plt
from scipy import integrate
def deriv(arrC,t,arrK):
    [CA,CB,CC,CD]=arrC
    [k1,k2]=arrK
    dCA=-k1*CA*CB
    dCB=-k1*CA*CB-k2*CB*CC
    dCC=k1*CA*CB-k2*CB*CC
    dCD=k2*CC*CB
    arrdC=[dCA,dCB,dCC,dCD]
    return arrdC
k1,k2=0.1,0.05
arrC0=[1.0,2.0,0.0,0.0]
tc=1.0/(k1*arrC0[1])
t=scipy.linspace(0.0,1,1000)
C=scipy.integrate.odeint(deriv,arrC0,t,args=([k1,k2],))
fig=plt.figure()
ax=fig.add_subplot(111);fig.show()
ax.plot(t,C[:,0],'r')
ax.plot(t,C[:,1],'b')
ax.plot(t,C[:,2],'m')
ax.plot(t,C[:,3],'g')
fig.canvas.draw()