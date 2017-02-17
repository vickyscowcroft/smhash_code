#!/usr/bin/env python


import numpy as np
import re
import sys
import pexpect ## Use this to control daophot
import shutil
from astropy.io import fits
import os

def init_aper_phot(fitsfile, channel):

	image = re.sub(".fits", "", fitsfile)

## Copy the daophot opt file
	if (channel == 1):
		shutil.copy("/Users/vs522/Dropbox/Python/smhash_code/daophot-spitzer-i1.opt", "daophot.opt")
	if (channel == 2):
		shutil.copy("/Users/vs522/Dropbox/Python/smhash_code/daophot-spitzer-i2.opt", "daophot.opt")

	shutil.copy("/Users/vs522/Dropbox/Python/smhash_code/photo-spitzer.opt", "photo.opt")

	print "Working on " + image
## Clean up previous runs

	extensions = ['.coo', '.lst', '.psf', '.nei', '.ap', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.als']
	for ext in extensions:
		if (os.path.isfile(image + ext)):
			os.remove(image+ext)

## Running daophot

	daophot = pexpect.spawn("daophot")

	hdulist = fits.open(fitsfile)
	prihdr = hdulist[0].header
	n_bcds = np.round(prihdr['medcov']) ### Changing this from the total number of BCDs to the median coverage - the coverage can be really non uniform so using the total number of BCDs in the stack is pretty useless.
	#n_bcds = prihdr['totalbcd'] ### Getting the correct number of BCDs that went into the stack
	#if (n_bcds > 20): #### Changing the threshold for the long exposure frames. Change it to do this a smarter way eventually
	#	daophot.expect("Command:")
	#	daophot.sendline("opt")
	#	daophot.expect("INPUT\):")
	#	daophot.sendline("")
	#	daophot.expect("OPT>")
	#	daophot.sendline("th=2")
	#	daophot.expect("OPT>")
	#	daophot.sendline("")
# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + image)

# find the stars
	daophot.expect("Command:")
	daophot.sendline("find")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline("1,1")
	#daophot.sendline("1," + str(np.sqrt(n_bcds)))
	daophot.expect("File for positions")
	daophot.sendline("")
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

	print "FIND complete"

## Aperture photometry
	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(image + '.coo')
	daophot.expect("Output file")
	daophot.sendline(image + '.ap')

	print "PHOT complete"

## Exit daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

	print "Initial daophot run complete."


