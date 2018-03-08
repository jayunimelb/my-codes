import numpy as np
import sys,pickle,gzip
from pylab import *

# input image file, which is an nd array
image = pickle.load(gzip.open(sys.argv[1]))

# give image resolution along x and y in arcminutes
dx, dy = float(sys.argv[2]),float(sys.argv[3])
try:
    conv_to_rads = int(sys.argv[4])
except:
    conv_to_rads = 1
# convert the resolution to radians
if conv_to_rads:
    dx,dy = np.radians(dx),np.radians(dy)

image_fft = np.fft.fft2(image)
#frequency along x and y directions
wx = np.fft.fftfreq(image.shape[0])*2*np.pi/dx
wy = np.fft.fftfreq(image.shape[1])*2*np.pi/dy

# normalising ffts (not involving phase information)

image_fft_norm = image_fft* dx*dy *(1/np.sqrt(2*np.pi))# *np.exp(-complex(0,1*wx*))
image_fft_norm = np.abs(image_fft_norm)

imshow(image_fft_norm);colorbar();show()

# 
