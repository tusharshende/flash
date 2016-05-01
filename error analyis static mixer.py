# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 23:00:33 2016

@author: Tushar
"""

#Error Estimation.I have taken "Hydrodynamics of Static Mixer" experiment to evaluate or to analyse the error within it.
import scipy ,numpy
import scipy.optimize,scipy.stats
import numpy.random
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels
import statsmodels.stats
import statsmodels.stats.stattools as stools

def fitdata(f,XData,YData,ErrData,pguess,ax=False,ax2=False):
    
    def error(p,XData,YData,ErrData):
        Y=f(XData,p)
        residuals=(Y-YData)/ErrData
        return residuals
    res=scipy.optimize.leastsq(error,pguess,args=(XData,YData,ErrData),full_output=1)
    (popt,pcov,infodict,errmsg,ier)=res
    
    perr=scipy.sqrt(scipy.diag(pcov))
    
    M=len(YData)
    N=len(popt)
    #residuals
    Y=f(XData,popt)
    residuals=(Y-YData)/ErrData
    meanY=scipy.mean(YData)
    squares=(Y-meanY)/ErrData
    squaresT=(YData-meanY)/ErrData
    
    SSM=sum(squares**2) #Corrected Sum of Squares
    SSE=sum(residuals**2)#Sum of Squares of Errors
    SST=sum(squaresT**2)#Totaal corrected sum of squares
    
    DFM=N-1 #degrees of freedom for model
    DFE=M-N #degrees of freedom for error
    DFT=M-1 #degrees of freedom total
    
    MSM=SSM/DFM #Mean squares for model (explained variance)
    MSE=SSE/DFE #Mean squares for error
    MST=SST/DFT #Mean squares for total
    
    R2=SSM/SST #proportion of explained varience
    R2_adj=1-(1-R2)*(M-1)/(M-N-1)#Adjusted R2
    
    #t-test to see if parameters are different from zero
    t_stat=popt/perr #t-statistic for popt different from zero
    t_stat=t_stat.real
    p_p=1.0-scipy.stats.t.cdf(t_stat,DFE) #should be low for good fit.
    z=scipy.stats.t(M-N).ppf(0.95)
    p95=perr*z
    
    #chisquared analysis on residuals
    chisquared=sum(residuals**2)
    degfreedom=M-N
    chisquared_red=chisquared/degfreedom
    p_chi2=1.0-scipy.stats.chi2.cdf(chisquared,degfreedom)
    stderr_reg=scipy.sqrt(chisquared_red)
    chisquare=(p_chi2,chisquared_red,degfreedom,R2,R2_adj)
    
    #Analysis of Residuals
    w,p_shapiro=scipy.stats.shapiro(residuals)
    mean_res=scipy.mean(residuals) #mean of all residuals
    stddev_res=scipy.sqrt(scipy.var(residuals)) #standard deviation of all residuals
    t_res=mean_res/stddev_res #t-statistic to test that was mean_res is zero.
    p_res=1.0-scipy.stats.t.cdf(t_res,M-1)
         #if p_res<0.05,null-hypothesis is rejected and mean is non-zero
         #should be high for good fit.
    # F-test on residuals
    F=MSM/MSE #explained/un-explained. Should be large
    p_F=1.0-scipy.stats.f.cdf(F,DFM,DFE)
         #if p_F<0.05,null hypothesis is rejected
         #i.e. R^2>0 and at least one of the fitting parameters >0.
    dw=stools.durbin_watson(residuals)
    resanal=(p_shapiro,w,mean_res,p_res,F,p_F,dw)
    if ax:
        formataxis(ax)
        ax.plot(YData,Y,'ro')
        ax.errorbar(YData,Y,yerr=ErrData,fmt='.')
        Ymin,Ymax=min((min(Y),min(YData))), max((max(Y),max(YData)))
        ax.plot([Ymin,Ymax],[Ymin,Ymax],'b')
        
        ax.xaxis.label.set_text('Data')
        ax.yaxis.label.set_text('Fitted')
        
        sigmaY,avg_stddev_data=get_stderr_fit(f,XData,popt,pcov)
        Yplus=Y+sigmaY
        Yminus=Y-sigmaY
        ax.plot(Y,Yplus,'C',alpha=0.6,linestyle='--',linewidth=0.5)
        ax.plot(Y,Yminus,'c',alpha=0.6,linestyle='--',linewidth=0.5)
        ax.fill_between(Y,Yminus,Yplus,facecolor='cyan',alpha=0.5)
        titletext='Parity plot for fit.\n'
        titletext += r'$r^2$ = %5.2f, $r^2_(adj)$= %5.2f'
        titletext +='$\sigma_{exp}$ = %5.2f, $\chi^2_{\nu}$=%5.2f, $p_{\chi^2}$=%5.2f,'
        titletext +='$\sigma_{ree}^{reg}$ = %5.2f'
        
        ax.title.set_text(titletext%(R2,R2_adj,avg_stddev_data,chisquared_red,p_chi2,stderr_reg))
        ax.figure.canvas.draw()
        
    if ax2: #test for homoscedasticity
        formataxis(ax2)
        ax2.plot(Y,residuals,'ro')
        
        ax2.xaxis.label.set_text('Fitted Data')
        ax2.yaxis.label.set_text('Residuals')
        
        titletext = 'Analysis of residuals\n'
        titletext +=r'mean= %5.2f, $p_{res}$= %5.2f, $p_{shapiro}$=%5.2f , $Durbin-watson$=%2.1f'
        titletext +='\n F=%5.2f,$p_F$ = %3.2e'
        ax2.title.set_text(titletext%(mean_res,p_res,p_shapiro,dw,F,p_F))
        
        ax2.figure.canvas.draw()
        
    return popt,pcov,perr,p95,p_p,chisquare,resanal
    
def get_stderr_fit(f,XData,popt,pcov):
    Y=f(XData,popt)
    listdY=[]
    for i in xrange(len(popt)):
        p=popt[i]
        dp=abs(p)/1e6 + 1e-20
        popt[i] +=dp
        Yi=f(XData,popt)
        dY=(Yi-Y)/dp
        listdY.append(dY)
        popt[i] -=dp
    listdY=scipy.array(listdY)#It is an array with N rows and M columns. N=len(popt), M=len(XData[0])
    left = scipy.dot(listdY.T,pcov)#left is an array with M rows and N columns matrix.. pcov is an array with N rows and N columns
    right=scipy.dot(left,listdY)#right is an M*M matrix.the diagonals of this are what we need.
    sigma2Y=right.diagonal()#sigmaY is the standard error of the fit and is a function of X
    mean_sigma2Y=scipy.mean(right.diagonal())
    M=XData.shape[1]
    N=len(popt)
    avg_stddev_data=scipy.sqrt(M*mean_sigma2Y/N)#If experimental error is constant ,then mean_sigma2Y=N/M*sig_dat**2
    sigmaY=scipy.sqrt(sigma2Y)
    return sigmaY,avg_stddev_data
    
def formataxis(ax):
    ax.xaxis.label.set_fontname('Georgia')
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontname('Georgia')
    ax.yaxis.label.set_fontsize(12)
    ax.title.set_fontname('Georgia')
    ax.title.set_fontsize(12)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
     
def import_data(xlfile,sheetname):
    df=pd.read_excel(xlfile,sheetname=sheetname)
    return df
    
def prepare_data(df,Criterion,Predictors,Error=False):
    Y=scipy.array(df[Criterion])
    if Error:
        ErrData=scipy.array(df[Error])
    else:
        ErrData=scipy.ones(len(Y))
    XData=[]
    for X in Predictors:
        X=list(df[X])
        XData.append(X)
    XData=scipy.array(XData)
    return XData,Y,ErrData
    
if __name__=='__main__':
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.show()
    
    fig2=plt.figure()
    ax2=fig.add_subplot(111)
    fig2.show()
    
    def f(X,p): #arbitrary function of two variables
     
        (x,)=X
        Y=p[0]/x
        return Y
        
    df=import_data('Exp15(Hydrodynamics of static mixer).xlsx','Sheet2') #importing excel file 
    XData,YData,ErrData=prepare_data(df,'YData',('x'),Error='err')
    #print YData
    N=1
    pguess=N*[0.0]
    
    popt,pcov,perr,p95,p_p,chisquare,resanal=fitdata(f,XData,YData,ErrData,pguess,ax=ax1,ax2=ax2)
    print popt
    #print chisquare
    #print resanal
