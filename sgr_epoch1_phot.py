#!/usr/bin/env python

import numpy as np
import re
import sys
import pexpect ## Use this to control daophot
import shutil
from astropy.io import fits
import os
fitsfile = sys.argv[1]

image = re.sub(".fits", "", fitsfile)

## Copy the daophot opt file
shutil.copy("/Users/vs/Dropbox/Python/smhash_code/daophot-spitzer-i1.opt", "daophot.opt")
shutil.copy("/Users/vs/Dropbox/Python/smhash_code/photo-spitzer.opt", "photo.opt")


## Clean up previous runs

extensions = ['.coo', '.lst', '.psf', '.nei', '.ap', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.als']
for ext in extensions:
	if (os.path.isfile(image + ext)):
		os.remove(image+ext)

## Running daophot

daophot = pexpect.spawn("/Users/vs/daophot/daophot")

hdulist = fits.open(fitsfile)
prihdr = hdulist[0].header
n_bcds = prihdr['totalbcd'] ### Getting the correct number of BCDs that went into the stack
# attach the image
daophot.expect("Command:")
daophot.sendline("at " + image)

# find the stars
daophot.expect("Command:")
daophot.sendline("find")
daophot.expect("Number of frames averaged, summed:")
daophot.sendline("1," + str(n_bcds))
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

print "Initial daophot run complete."
print "Please make psf by hand for this epoch."

