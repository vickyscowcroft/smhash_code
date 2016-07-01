#!/usr/bin/env python

import sys
import numpy as np

def apply_apcor(input, apcor, channel):

	stripped = input[:-4]

	print channel

	id, xc, yc, mag, err = np.loadtxt(input, skiprows=3, usecols=(0, 1, 2, 3, 4 ), unpack='TRUE')

## Zero point correction
	if (channel == 1):
		zp = -25. + 18.80
	## zmags changed for mosaicked data
	elif (channel == 2):
		zp = -25. + 18.31
	## zmags changed for mosaicked data
	
	else:
		print ("Invalid channel")
		sys.exit(1)

## Aperture correction to 3, 3, 7

	apc =  mag + apcor + zp

## Aperture correction to 10, 12, 20

## Updated 27/06/2016 for new S19.2 pipeline

	flux = 10**(-apc/2.5)
	if (channel == '1'):
		flux = flux * 1.125
	if (channel == '2'):
		flux = flux * 1.123
	
	apc = -2.5*np.log10(flux)

## Save output to file. Still needs to be corrected for location at this point.

	output = stripped + '.apc'

	np.savetxt(output, np.column_stack((id, xc, yc, apc, err)), fmt= "%d %.2f %.2f %.3f %.3f")



