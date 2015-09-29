#!/usr/bin/env python

import re
import shutil
from astropy.io import fits

def spitzer_flux2dn(image, newname=""):
	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	print newname
	shutil.copy(image, newname)
	hdulist = fits.open(newname)
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	exptime = prihdr['exptime']
	fluxconv = prihdr['fluxconv']
	scidata *= exptime / fluxconv
	return(0)
	
	
	