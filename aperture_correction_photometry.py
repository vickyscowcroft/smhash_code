#!/usr/bin/env python

import numpy as np
import re
import sys
import pexpect ## Use this to control daophot
import shutil
import os
import glob

image = sys.argv[1]
## image in this case should be the flux image
isdn = re.search('_dn', image)
if (isdn != None):
	print "you need to run this on the flux image, not the counts image"
	exit

## Do the aperture photometry on the flux image in the standard aperture



if(os.path.isfile("apcor_photo.opt")): os.remove("apcor_photo.opt")
shutil.copy("/Users/vs/Dropbox/Python/smhash_code/apcor_photo.opt", "apcor_photo.opt")

if(os.path.isfile(image + '.apc')): os.remove(image + '.apc')

## Aperture phot
daophot = pexpect.spawn("/Users/vs/daophot/daophot")
daophot.logfile = sys.stdout
daophot.expect("Command:")
daophot.sendline("at " + image)

daophot.expect("Command:")
daophot.sendline("phot")
daophot.expect("File with aperture radii")
daophot.sendline("apcor_photo.opt")
daophot.expect("PHO>")
daophot.sendline("")
daophot.expect("Input position file")
daophot.sendline(image + '_dn.lst')
daophot.expect("Output file")
daophot.sendline(image + '.apc')

daophot.expect("Command:")
daophot.sendline("exit")
daophot.close(force=True)

## Match up the aperture phot file with the psf phot file
# faking the initial file as they are the same image

apcor_mch = open('apcor.mch', 'w')

apcor_mch.write("'{0:s}.apc '  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))
apcor_mch.write("'{0:s}_dn.cal '  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))

apcor_mch.close()


## get the raw file from daomaster
if(os.path.isfile('apcor.raw')): os.remove('apcor.raw')


daomaster = pexpect.spawn("/Users/vs/daophot/daomaster")
daomaster.expect("File with list of input files")
daomaster.sendline('apcor.mch')
daomaster.expect("Minimum number, minimum fraction, enough frames")
daomaster.sendline("2, 1, 2")
print "2, 1, 2"
daomaster.expect("Maximum sigma")
daomaster.sendline("99")
## desired degrees of freedom:
daomaster.expect("Your choice")
daomaster.sendline("6")
daomaster.expect("Critical match-up radius")
daomaster.sendline("10")

for radius in range (10,-1, -1):
	daomaster.expect("New match-up radius")
	daomaster.sendline(str(radius))
## Only need the raw here
print "making output files"
daomaster.expect("Assign new star IDs")
daomaster.sendline("n")
daomaster.expect("A file with mean magnitudes and scatter")
daomaster.sendline("n")
daomaster.expect("A file with corrected magnitudes and errors")
daomaster.sendline("n")
daomaster.expect("A file with raw magnitudes and errors")
daomaster.sendline("y")
daomaster.expect("Output file name")
daomaster.sendline("apcor.raw")
daomaster.expect("A file with the new transformations")
daomaster.sendline("e")
daomaster.close(force=True)

new_raw  = glob.glob('apcor.raw*')
print new_raw
shutil.move( str(new_raw[0]), "apcor.raw")




