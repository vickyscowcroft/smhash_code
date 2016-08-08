#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp
from astropy.stats import sigma_clip
import aplpy

## Updated 28/09/2015 to use sigmaclipping to remove bad points rather than relying on a hard wired magnitude range

## Usage: argv[1] = input raw file
## input[2] = sigma clipping limit
## See smhash documentation for full details

def calc_apcor(flux_image, input, sigma, target):

	xc, yc, apc, eapc, alf, ealf = np.loadtxt(input, skiprows=3, usecols=(1, 2, 3, 4, 5, 6), unpack='TRUE')

	mp.close()

	fig = mp.figure(figsize=(10,10))

	#target = target_name[0:11] + '.' + target_name[12:20]
	print target
	
	image_stem = flux_image[0:7]
	fitsfile = flux_image + '.fits'
		
	## Cutting the difference arrays so they only include stars in the good image region:
	#difference  = apc - alf

	apc2 = apc[(ealf < 0.1)]
	alf2 = alf[(ealf < 0.1)]
	ealf2 = ealf[(ealf < 0.1)]
	eapc2 = eapc[(ealf < 0.1)]
	difference = apc2 - alf2

#av_diff = np.average(difference)
#sdev_diff = np.std(difference)

## Using astropy.stats.sigma_clip
## This version returns a masked array rather than a cut array

	alf_sorted = np.argsort(alf2)
	brightest = alf_sorted[0:50]
	print brightest

	clipped2 = sigma_clip(difference[brightest], sig=sigma, iters=5)

	av_diff = np.ma.mean(clipped2)
	sdev_diff = np.ma.std(clipped2)


	total_err = np.sqrt(eapc2**2 + ealf2**2)

	axp1 = mp.subplot(211)

#axp1.errorbar(apc, difference, yerr = total_err, color='grey', ls='none')
#axp1.plot(apc, difference, 'k.', ls='none')
#axp1.axhline(av_diff, color='r', ls='--')



	axp1.errorbar(alf2[brightest], clipped2, yerr = total_err[brightest], color='grey', ls='none')
	axp1.plot(alf2[brightest], clipped2, 'k.', ls='none')
	axp1.axhline(av_diff, color='r', ls='--')
	axp1.axhline(av_diff+2*sdev_diff, color='b', ls='--')
	axp1.axhline(av_diff-2*sdev_diff, color='b', ls='--')

	axp2 = mp.subplot(212)

	axp2.plot(alf2, apc2, 'k.', ls='None')
	axp2.plot(alf2[brightest], apc2[brightest], 'r.', ls='None')

	mp.show()

	return(av_diff, sdev_diff)
	