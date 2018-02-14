import numpy as np
import pickle,gzip
import sys
from pylab import *
fname = sys.argv[1]
data = pickle.load(gzip.open(fname))
results = data['results']
x = np.zeros(len(results.keys()));y = np.zeros(len(results.keys()))

#print fname;quit()
#max_lk = max(results.values())
#min_lk = min(results.values())
#print np.sort(results.keys());quit()
#print len(results.keys())-1;quit()
#def get_likelihood():


for j,key in enumerate(sorted(results.keys())):
	x[j] = float(key)
	y[j] = results[key]
	#inds = np.where(x<8.)[0]
	i =i+1
#print x
#print y;
#plot(x,y);show();quit()
y -= max(y)
y = np.exp(y)
#print x
#print y
#quit()
	#return x,y
#print x
#print y
#quit()
plot(x,y,lw =3);show()
