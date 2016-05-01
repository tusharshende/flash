# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 19:58:53 2016

@author: Eeshani Godbole
"""

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt

P=1 #atm #total pressure
kLa_m=1.0367 #s-1 #overall mass transfer coefficient for Methane
kLa_c=0.0253e-2 #m/s #overall mass transfer coefficient for Carbon dioxide
Hm=1.4e-3 #M/atm #Henry's law constant for Methane '''Wilhelm et al [1977]
Hc=3.4e-2 #M/atm #Henry's law constant for Carbon dioxide
L=5 #m #length of absorption column
G=10 #mol/s #total gas flow rate
L=500 #mol/s #total liquid flow rate
A=0.06 #m2 #cross sectional area of the column
kGW=1.367e-3 #m/s
a=100 #m2/m3
Pwsat=0.0418 #atm #saturation pressure of water at 25 deg celcius
#Flow rates have flowrates in kmol/s flowrates

def gm(Gm,z):
    
    
    Gm0=Gm[0]
    Pm=P*(Gm/(Gm+Gc+Gw))
    CmL=Lm/L
    Gm1=