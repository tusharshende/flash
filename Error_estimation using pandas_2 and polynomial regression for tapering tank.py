# -*- coding: utf-8 -*-
"""
Created on Tue Feb 09 18:52:53 2016

@author: Omkar Mehta
"""

# import formula api as alias smf
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd
import statsmodels.api as sm
import os
import matplotlib.pyplot as plt
pre = os.path.dirname(os.path.realpath(__file__))
fname = 'Tapering tank.xlsx'
path = os.path.join(pre, fname)
df = pd.read_excel(path,sheetname='Data')
df.head()

y = df.h  # response
print y.shape
X = df.t  # predictor
print X.shape
X = sm.add_constant(X)  # Adds a constant term to the predictor
print X.head()

# plot x vs Ydata
plt.figure(figsize=(6 , 6))
plt.scatter(df.t, df.h, s=10, alpha=0.3)
plt.xlabel('t')
plt.ylabel('h')
def exp(t):
	return -np.exp(-0.003639933*t)
# points linearlyd space on lstats
x = pd.DataFrame({'t': np.linspace(df.t.min(), df.t.max(), 100)})
# Exponential
exp1 = smf.ols(formula='h ~  exp(t)', data=df).fit()
plt.plot(x.t, exp1.predict(x), 'r-', label='Poly n=1 $R^2$=%.2f' % exp1.rsquared, 
         alpha=0.9)
a=exp1.summary()
print 'a', a

