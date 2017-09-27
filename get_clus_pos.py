import pickle,gzip
import scl_cmb as scl_cmb
import sky_local
import numpy as np
from pylab import *
#sims = scl_cmb.SPTpol_map_using_SPTpol_pipeline()

#sims = scl_cmb.simulations()

catalog =  'map_150.pkl.gz_cutouts_150ghz_no_downsamplingredmapper_clusters'
cat_data = pickle.load(gzip.open('%s'%(catalog)))['cutouts']

data_map = pickle.load(gzip.open('map_150.pkl.gz'))[0]
imshow(data_map)
print data_map.shape
for key in cat_data.keys():
	ra, dec = key[0], key[1]
#	ra, dec = np.degrees(ra), np.degrees(dec)
	y, x= sky_local.ang2Pix([ra,dec], [0, -57.5], 0.5, np.asarray(data_map.shape), proj = 0) # these are the x and y co-ordinates we get on imshow space
#	imshow(data_map)
	plot(x,y,'k*')

show()