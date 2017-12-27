import numpy as np
import pickle,gzip
import sys
from pylab import *
fname = sys.argv[1]
data = pickle.load(gzip.open(fname))
results = data['results']
x = np.zeros(len(results.keys())-1);y = np.zeros(len(results.keys())-1)

#print fname;quit()
#max_lk = max(results.values())
#min_lk = min(results.values())
#print np.sort(results.keys());quit()
#print len(results.keys())-1;quit()
#def get_likelihood():
i=0
#print sorted(results.keys());quit()
for j,key in enumerate(sorted(results.keys())):
	#if float(key) == 0.000:
		#continue
	if float(key) == 5.0:
		continue
#print float(key)
	#i =0
	x[i] = float(key)
	
	y[i] = results[key]
	#y[i] = (results[key] - min_lk)/(max_lk -min_lk)
#	x = np.sort(x)
	#print results.keys()
	#for k,z in enumerate(x):
		#y[k] = results['%.3f'%(z)]
#print x;quit()
	inds = np.where(x<8.)[0]
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
