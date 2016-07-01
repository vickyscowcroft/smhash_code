#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit

def ccm_optical(wavelength, Rv):
	x = 1.0 / wavelength
	y = x - 1.82
	ax = 1.0 + 0.17699*(y) - 0.50447*(y**2) -0.02427*(y**3) + 0.72085*(y**4) + 0.01979*(y**5) - 0.77530*(y**6) + 0.32999*(y**7)
	bx = 1.41338*y + 2.28305*(y**2) + 1.07233*(y**3) - 5.38434*(y**4) - 0.62251*(y**5) + 5.30260*(y**6) - 2.09002*(y**7)
	AlAv = ax + (bx / Rv)
	return(AlAv)
	
def ccm_nearir(wavelength, Rv):
	x = 1.0 / wavelength
	ax = 0.574*(x**1.61)
	bx = -0.527*(x**1.61)
	AlAv = ax + (bx / Rv)
	return(AlAv)
	
def indebetouw_ir(wavelength): ## For use with 1.5 microns onwards
	logw = np.log10(wavelength)
	logwsq = logw**2
	logAlAk = 0.61 - 2.22*logw + 1.21*logwsq
	AlAk = 10**(logAlAk)
	return(AlAk)

## The fit_reddening module takes all the uncorrected distance moduli
## and returns muo, Av, muo_err, Averr
## which are the dereddened distance modulus, the Av and their uncertainties

	
def fit_reddening(b,v,r,i,j,h,k,i1,i2):
	## The input values are the uncorrected distance moduli
	alam = []
	mus = []
	wavelengths = [0.445, 0.551, 0.658, 0.806, 1.22, 1.63, 2.19, 3.545, 4.442]
	mus = b,v,r,i,j,h,k,i1,i2
	mus = np.array(mus)
	num = mus[mus>0].size	
	if num >= 3:
		Ak = ccm_nearir(2.164,3.1)
		for band in range(0,4):
			if mus[band] > 0:
				alam.append(ccm_optical(wavelengths[band], 3.1))
			else:
				alam.append(0.0)
		for band in range(4,7):
			if mus[band] > 0:
				alam.append(ccm_nearir(wavelengths[band],3.1))
			else:
				alam.append(0.0)
		for band in range(7,9):
			if mus[band] > 0:
				alam.append(indebetouw_ir(wavelengths[band])*Ak)
			else:
				alam.append(0.0)
		alam = np.array(alam)
		#print mus
		#print alam
		def redfit(alam, muo, Av):
			return muo + Av*alam
		popt, pcov = curve_fit(redfit, alam[alam>0], mus[mus>0])
		muo = popt[0]
		Av = popt[1]
		#print muo, Av
		muo_err = np.sqrt(pcov[0, 0])
		Averr = np.sqrt(pcov[1, 1])
		#print muo_err
		#print Averr
		return(muo, Av, muo_err, Averr)
	else:
		return(i1, 0.00, 0.00, 0.00)
		
		

	
