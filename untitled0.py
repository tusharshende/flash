'''
Name-Tushar A Shende
Date-10/01/2016 
'''
#for water(in gas phase)
import numpy as np
A=50
kla=0.00005

pi=0.2
Hi=0.001
def f(x):
    return -kla*A*((pi/Hi)-x)
#define initial condition
x0=100
dt=0.1
t=np.linspace(0,1,int(1/dt)+1)    
x=np.zeros(len(t))
x[0]=x0
for i in xrange(1,len(t)):
    x[i]=x[i-1]+f(x[i-1])*dt
print x[i]

