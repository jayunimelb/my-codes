"""

Converts loglikelihood to likelihood and fit gaussian to calculate 1 sigma significance!


"""

import numpy as np
import sys,pickle,gzip 
from scipy.optimize import curve_fit
from pylab import *
fname = sys.argv[1]
data = pickle.load(gzip.open(fname))
results = data['results']
x = np.zeros(len(results.keys())-1);y = np.zeros(len(results.keys())-1)#x = np.zeros(len(results.keys()));y = np.zeros(len(results.keys())) 

j = 0
for key in sorted(results.keys()):
	if(float(key) != 5.):
		x[j] = float(key)
		y[j] = results[key]
		j = j+1

y -= max(y); y = np.exp(y)
plot(x,y)#;show()
inds = np.where(y >=0.6)
x = x[inds]
y = y[inds]
n = len(x)
mean = sum(x*y)/n
sigma = sum(y*(x-mean)**2)/n

def fit_gaussian(x,a, x0,sigma):
	return a*np.exp(-(x-x0)**2/(2*sigma**2))


popt, pcov = curve_fit(fit_gaussian,x,y,p0 =[1,mean,sigma])
x_new = np.arange(min(x),max(x),0.1)
#plot(x,fit_gaussian(x,*popt), label = 'fit')
print popt
plot(x_new,fit_gaussian(x_new,*popt),'r', label = 'fit')

plot(x,y,'b')
show()