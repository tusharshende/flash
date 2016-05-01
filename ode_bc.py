import scipy
import scipy.optimize,scipy.integrate
from matplotlib import pyplot as plt

#step 1: write the function deriv
def deriv(Y,z,X0,Y0,h,A,B):
    n=len(Y)
    sumY=sum(Y)
    Ky=A/sumY #This is a Diagonal Matrix

    X=Y+(X0-Y0)
    sumX=sum(X)
    Kx=B/sumX #This is a Diagonal Matrix

    DrvF=scipy.dot(Ky,Y)-scipy.dot(Kx,X) #column matrix

    dYbydz=-scipy.dot(h,DrvF) #Column matrix n*1
    dYbydz=dYbydz.reshape((n,)) #NOw it is a vector of length n and not a n*1 matrix
    return dYbydz

Y0=scipy.array([0.5,0.5,1.0])
XL=scipy.array([0.0,0.0,10.0])

h=scipy.diag(scipy.array([1.0,1.5,2.0]))
A=scipy.diag(scipy.array([2.0,1.0,1.5]))
B=scipy.diag(scipy.array([1.0,1.0,1.0]))

L=1.0
n=50
z=scipy.linspace(0,L,n)

def error(X0,z,XL,Y0,h,A,B):
    soln=scipy.integrate.odeint(deriv,Y0,z,args=(X0,Y0,h,A,B))
    YL=soln[-1,:]
    XL_sm=YL+(X0-Y0)
    error=XL_sm-XL
    return error

X0_guess=XL+0.0
X0_opt,err=scipy.optimize.leastsq(error,X0_guess,args=(z,XL,Y0,h,A,B))
YY=scipy.integrate.odeint(deriv,Y0,z,args=(X0_opt,Y0,h,A,B))
XX=YY+(X0_opt-Y0)

#Plotting the driving force
Kyy=scipy.array([A/sum(Y) for Y in YY])
Kxx=scipy.array([B/sum(X) for X in XX])
KyyYY=scipy.array([scipy.dot(Kyy[i],YY[i]) for i in xrange(len(YY))])
KxxXX=scipy.array([scipy.dot(Kxx[i],XX[i]) for i in xrange(len(XX))])
                                  
fig=plt.figure()
for i in xrange(len(Y0)):
    ax=fig.add_subplot(1,len(Y0),i+1)
    PotY=KyyYY[:,i]
    PotX=KxxXX[:,i]
    maxPot=max([max(PotY),max(PotX)])
    ax.plot(PotY,z,'r')
    ax.plot(PotX,z,'b')
    ax.title.set_text('Driving Force for Component %d'%(i))
    ax.xaxis.label.set_text('Potential')
    ax.yaxis.label.set_text('z')
    ax.axis([0.0,maxPot,0.0,L])
    ax.legend(['PhaseY','PhaseX'])
fig.canvas.draw()                                  
fig.show()

        
                                  


