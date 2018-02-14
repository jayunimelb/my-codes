"""


"""
import numpy as np
import scipy.optimize as opt
from pylab import *
def twoD_gaussian((x,y), amplitude,x_0,y_0,sigma_x,sigma_y,theta = 0, offset =0):

	"""
	Given x,y as numpy meshgrid, amplitude of the gaussian
	initial guesses for x_0, y_0 as the central values of the gaussian
	sigma_x and sigma_y are the widths along x and y direction
	theta --> defines the whether the gaussian distributions principal axes are X and Y and theta =0 is the intial guess principal axes alinged along X and Y
	offset --> by default is zero

	"""
	x_0 = float(x_0)
	y_0 = float(y_0) 
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = offset + amplitude*np.exp( - (a*((x-x_0)**2) + 2*b*(x-x_0)*(y-y_0) + c*((y-y_0)**2)))
	return g.ravel()



x = np.linspace(0, 200, 201)
y = np.linspace(0, 200, 201)
x, y = np.meshgrid(x, y)

data = twoD_gaussian((x, y), 1, 100, 100, 20, 40)
initial_guess = (3,100,100,20,40,0,10)
data_noisy = data + 0.2*np.random.normal(size=data.shape)
popt, pcov = opt.curve_fit(twoD_gaussian, (x, y), data_noisy, p0=initial_guess)
data_fitted = twoD_gaussian((x, y), *popt)
subplot(211);imshow(data_noisy.reshape(len(x),len(y)));colorbar();title('raw data')
subplot(212);imshow(data_fitted.reshape(len(x),len(y)));colorbar();title('fitted data');show()

imshow(data_fitted.reshape(len(x),len(y))-data.reshape(201, 201));colorbar();show()