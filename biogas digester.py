# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 23:18:47 2016

@author: Tushar
"""
''' Coding for Anaerobic Digestion Model'''
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
'''Given is the animal waste conversion system  which converts
cow manure into biogas utilizing Anaerobic Digestion (AD).We consider a reactor, where the feed is described by a volumetric feed rate  Vf [L/ d] (control input) with a
given concentration Z_vsf of volatile solids (disturbance).
The liquid level of the reactor is made constant by the use of a weir system, and it is possible
to control the reactor temperature T accurately using electric heating (potential control input). 
We described the operations of the bio-reactor by 4 states, viz. hydrolysis,acidogenesis,acetogenesis and methanogenesis.
Terminologies
bvs= biodegradable volatile solids, vfa=volatile fatty acids,
Xa=Acetogenic bacteria, Xm=Methanogenic bacteria       
'''

Vf=50 #volumetric feed rate (L/d)
V=250 #Reactor Volume (in Litres)
Z_vsf=32.4 #Feed concentration of volatile solids(g/l)
theta_Xbvs=theta_Xvfa=1 #residence time correction for solids
theta_Xa=theta_Xm=2.9 #Correction of residence time for bacteria due to nonideal flow
Y_bvsXa=3.90 #(Inverse) yield: consumption of bvs per growth of bacteria
Y_vfaXa=1.76 #(Inverse) yield: production of vfa per growth of bacteria
Y_vfaXm=31.7 #(Inverse) yield: consumption of vfa per growth of bacteria
Y_CH4Xm=26.3 #(Inverse) yield: production of methane per growth of bacteria
K_bvs=15.5 #Half-velocity constant for bvs substrate
K_vfa=3.0 #Half-velocity constant for vfa substrate
muhat_35=0.326 #Maximal growth rate at T = 35 degrees C,
alpha_muhat=0.013 #Temperature sensitivity of maximal growth rate, valid T--[20,60] degrees C
kd_a=kd_m=0.02 #Death rate constants for acetogenic and methanogenic bacteria
b0=0.25 #Fraction biodegradable volatile solids in volatile solids feed
af=0.69 #Fraction volatile fatty acids in biodegradable volatile solids feed

def muhat(T):
    return muhat_35+alpha_muhat*(T-35)
print muhat(50)

def rxn(Z,t):
    Ra=(muhat(50)/(1+(K_bvs/Z[0])))*Z[2]
    Rm=(muhat(50)/(1+(K_vfa/Z[1])))*Z[3]
    
    R_bvs=-Y_bvsXa*Ra
    R_vfa=Y_vfaXa*Ra-Y_vfaXm*Rm
    R_Xa=Ra-kd_a*Z[2]
    R_Xm=Rm-kd_m*Z[3]
   
    #Feed Concentrations of states
    Z_bvsf=b0*Z_vsf
    Z_vfaf=af*Z_bvsf
    Z_Xaf=0
    Z_Xmf=0
    
    #Differential Equations for each state
    dz_bvsdt=(1/theta_Xbvs)*(Vf/V)*(Z_bvsf-Z[0])+R_bvs
    dz_vfadt=(1/theta_Xvfa)*(Vf/V)*(Z_vfaf-Z[1])+R_vfa
    dz_Xadt=(1/theta_Xa)*(Vf/V)*(Z_Xaf-Z[2])+R_Xa
    dz_Xmdt=(1/theta_Xm)*(Vf/V)*(Z_Xmf-Z[3])+R_Xm
    
    return (dz_bvsdt,dz_vfadt,dz_Xadt,dz_Xmdt)
    
t=np.linspace(0,100,10)
Z0=[5.81,1.13,1.32,0.39] #Initial concentrations of solids and bacterias
Conc=odeint(rxn,Z0,t) 
print (Conc)

#COncentration vs time graph
Z_bvs=Conc[:,0]
Z_vfa=Conc[:,1]
Z_Xa=Conc[:,2]   
Z_Xm=Conc[:,3]  
plt.figure(1)
plt.plot (t,Z_bvs)
plt.plot (t,Z_vfa)
plt.plot (t,Z_Xa)
plt.plot (t,Z_Xm)
plt.xlabel('time (sec)')
plt.ylabel('COncentration')
plt.legend(['Z_bvs','Z_vfa','Z_Xa','Z_Xm'])

#pH vs time graph
pH_Z_bvs=-np.log(Conc[:,0])
pH_Z_vfa=-np.log(Conc[:,1])
pH_Z_Xa=-np.log(Conc[:,2])   
pH_Z_Xm=-np.log(Conc[:,3])
plt.figure(2)
plt.plot (t,pH_Z_bvs) 
plt.plot (t,pH_Z_vfa)
plt.plot (t,pH_Z_Xa)
plt.plot (t,pH_Z_Xm)
plt.xlabel('time (sec)')
plt.ylabel('pH')
plt.legend(['pH_Z_bvs','pH_Z_vfa','pH_Z_Xa','pH_Z_Xm'])    


    
    
    