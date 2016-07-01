#!/usr/bin/env python

import re
import shutil
from astropy.io import fits
import numpy as np

def spitzer_flux2dn(image, newname=""):
	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	print newname
	shutil.copy(image, newname)
	hdulist = fits.open(newname, mode='update')
#	hdulist.fileinfo()
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	exptime = prihdr['exptime']
	fluxconv = prihdr['fluxconv']
	scidata *= exptime / fluxconv
	return(0)
	
def spitzer_flux2dn_mos(image, median, new_target_stem, channel):
	newname = new_target_stem + "_" + channel + "_dn.fits"
	shutil.copy(image, newname)
	hdulist = fits.open(newname, mode='update')
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	exptime = prihdr['exptime']
	fluxconv = prihdr['fluxconv']
	scidata *= exptime / fluxconv
	newmosname = new_target_stem + "_" + channel + "_med_dn.fits"
	shutil.copy(median, newmosname)
	hdulist = fits.open(newmosname, mode='update')
	scidata = hdulist[0].data
	scidata *= exptime / fluxconv
	return(0)
	
def flux_to_dn_expmap(image, map, new):
	maphdu = fits.open(map)
	mapdata = maphdu[0].data
	median_exptime = np.median(mapdata[mapdata!=0])
	shutil.copy(image, new)
	newhdu = fits.open(new, mode='update')
	scidata = newhdu[0].data
	prihdr = newhdu[0].header
	fluxconv = prihdr['fluxconv']
	scidata *= median_exptime / fluxconv
	return(0)
	
		
	