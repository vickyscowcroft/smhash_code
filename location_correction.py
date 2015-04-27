#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import re

## Usage: input data file = sys.argv[1], correction image = sys.argv[2]


photFile = open(sys.argv[1], "r")


corImage = fits.open(sys.argv[2])
corData = corImage[0].data

newmag = []
corr = []


id, xc, yc, mag, err = np.loadtxt(photFile, usecols=(0, 1, 2, 3, 4), unpack='TRUE')
flux = 10**(mag / -2.5)
for star in id:
	star = star - 1.
	xcoord = int(floor(xc[star]) - 1.)
	ycoord = int(floor(yc[star]) - 1.)
	correction = corData[ycoord, xcoord]
	corr.append(correction)
	#print correction
	"""
	
	
	if (np.isnan(correction) == False):
	#	print 'correction is good'
		flux[star] = flux[star] * correction
		newmag.append(-2.5 * log10(flux[star]))
	elif (np.isnan(correction) == True):
	#	print 'correction is bad'
		newmag.append(mag)

	"""
corr = np.array(corr)

flux = flux * corr
newmag = -2.5 * log10(flux)


#newmag = np.array(newmag)


	
output_name = re.sub('.off', '.cal',sys.argv[1])
np.savetxt(output_name, np.column_stack((id[np.isnan(corr) == False], xc[np.isnan(corr) == False], yc[np.isnan(corr) == False], newmag[np.isnan(corr) == False], err[np.isnan(corr) == False], corr[np.isnan(corr) == False])), fmt= "%d %.2f %.2f %.3f %.3f %.3f")
		
print ("Finished")


