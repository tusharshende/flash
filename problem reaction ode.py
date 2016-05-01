# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 23:52:00 2016

@author: Tushar
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def rxn(C,t):
    Ca=C[0]
    Cb=C[1]
    k=2.0
    dAdt=-k*Ca
    dBdt= k*Ca
    return (dAdt,dBdt)
print (rxn([0.2,0.5],0.5))

t=np.linspace(0,5,100)
C0=[1,0]
C=odeint(rxn,C0,t)
print (C)

plt.plot(t,C[:,0],'r--',linewidth=2.0)
plt.plot(t,C[:,1],'b-',linewidth=2.0)
plt.xlabel('time (sec)')
plt.ylabel('COncentration')
plt.legend(['Ca','Cb'])