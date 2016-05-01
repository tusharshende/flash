import scipy
import scipy.integrate
import matplotlib.pyplot as plt

def derivative(Y,x,a,b):
    [y,z]=Y
    dybydx=z
    dzbydx=a*y+b*scipy.exp(x)
    return [dybydx,dzbydx]

y0,z0=1.0,2.0
Y0=[y0,z0]

x=scipy.linspace(0.0,1.0,30)
a,b=4.0,1.0

solution=scipy.integrate.odeint(derivative,Y0,x,args=(a,b))

y=solution[:,0]
z=solution[:,1]

fig=plt.figure();fig.show()
axy,axz=fig.add_subplot(121),fig.add_subplot(122)
axy.title.set_text('y vs x');axz.title.set_text('z vs x')
axy.xaxis.label.set_text('x');axy.yaxis.label.set_text('y')
axz.xaxis.label.set_text('x');axz.yaxis.label.set_text('z')
axy.plot(x,y,'r')
axz.plot(x,z,'r')

y_theo=lambda x:(5.0/4.0)*scipy.exp(2*x)+(1.0/12.0)*scipy.exp(-2*x)-scipy.exp(x)/3.0
axy.plot(x,y_theo(x),'c',linewidth=5.0,alpha=0.5)
fig.canvas.draw()



