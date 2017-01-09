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

def calc_apcor(flux_image, input, sigma, target, error_limit):

	xc, yc, apc, eapc, alf, ealf = np.loadtxt(input, skiprows=3, usecols=(1, 2, 3, 4, 5, 6), unpack='TRUE')

	mp.close()

	fig = mp.figure(figsize=(10,5))

	#target = target_name[0:11] + '.' + target_name[12:20]
	print target
	
	#image_stem = flux_image[0:7]
	fitsfile = flux_image + '.fits'
		
	## Cutting the difference arrays so they only include stars in the good image region:
	#difference  = apc - alf

	apc2 = apc[(ealf < error_limit)]
	alf2 = alf[(ealf < error_limit)]
	ealf2 = ealf[(ealf < error_limit)]
	eapc2 = eapc[(ealf < error_limit)]
	difference = apc2 - alf2
	diff2 = apc - alf
	err2 = np.sqrt(eapc**2 + ealf**2)

#av_diff = np.average(difference)
#sdev_diff = np.std(difference)

## Using astropy.stats.sigma_clip
## This version returns a masked array rather than a cut array

	alf_sorted = np.argsort(alf2)
	brightest = alf_sorted[0:50]
	print brightest

	clipped2 = sigma_clip(difference[brightest], sigma=sigma, iters=None)

	av_diff = np.ma.mean(clipped2)
	sdev_diff = np.ma.std(clipped2)


	total_err = np.sqrt(eapc2**2 + ealf2**2)
	
	full_sample = np.zeros(len(alf))
	
	alf2_sample = np.zeros(len(clipped2))
	
	if(len(alf)>50):
		altered_mask = np.append(clipped2.mask,np.zeros(len(alf)-50))
	else:
		altered_mask = clipped2.mask
	#clipped_sample = np.zeros(len(clipped2))
	
	## fix this later for samples where clipped2 array not same length as original alf array
	
	full_sample[altered_mask==1] = 1 ## if full_sample == 0, then good star
	full_sample[(ealf>error_limit) | (eapc > error_limit)] = 2
	full_sample[(diff2 - av_diff) > error_limit] = 3 ## high range
	full_sample[(diff2 - av_diff) < -error_limit] = 4 ## low range

	alf2_sample[clipped2==1] = 1 
	alf2_sample[(ealf[brightest]>error_limit) | (eapc[brightest] > error_limit)] = 2
	alf2_sample[(diff2[brightest] - av_diff) > 2*error_limit] = 3 ## high range
	alf2_sample[(diff2[brightest] - av_diff) < -2*error_limit] = 4 ## low range


	axp1 = mp.subplot(111)


	mp.ylim(av_diff - 2*error_limit, av_diff + 2*error_limit)
	
	axp1.errorbar(alf[full_sample==2], diff2[full_sample==2], yerr = err2[full_sample==2], color='r', ls='none')
	axp1.errorbar(alf[full_sample==1], diff2[full_sample==1], yerr = err2[full_sample==1], color='b', ls='none')
	axp1.errorbar(alf2[alf2_sample==0], clipped2[alf2_sample==0], yerr = total_err[alf2_sample==0], color='grey', ls='none')
	
	axp1.plot(alf[full_sample==2], diff2[full_sample==2], 'ro', ls='none', label='Error cut', alpha=0.5)
	axp1.plot(alf[full_sample==1], diff2[full_sample==1], 'bo', ls='none', label='sigma clipped', alpha=0.5)
	axp1.plot(alf2[alf2_sample==0], clipped2[alf2_sample==0], 'ko', ls='none', label='Good stars')

### Plotting points beyond axis range
	if(len(alf[full_sample==3])>0):
		axp1.errorbar(alf[full_sample==3], np.ones_like(alf[full_sample==3])*(av_diff + error_limit), yerr = error_limit, color='r', ls='none', lolims=True)
	if(len(alf[full_sample==4])>0):
		axp1.errorbar(alf[full_sample==4], np.ones_like(alf[full_sample==4])*(av_diff - error_limit), yerr = error_limit, color='r', ls='none', uplims=True)

	
	axp1.axhline(av_diff, color='k', ls='--')
	axp1.axhline(av_diff+2*sdev_diff, color='b', ls='--')
	axp1.axhline(av_diff-2*sdev_diff, color='b', ls='--')
	mp.title(target + ': av_diff = ' + str(np.round(av_diff, decimals=3)) + ' sdev_diff = ' + str(np.around(sdev_diff, decimals=3)) + ' nstars = ' + str(len(alf[full_sample==0])))
	
	mp.legend(loc='upper left')


	#mp.show()
	
	mp.savefig(target + '_apcor.pdf')

	return(av_diff, sdev_diff)
	