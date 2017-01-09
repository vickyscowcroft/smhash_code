#!/usr/bin/env python

### Last edit November 8th 2016
### This script runs the initial photometry then stops when you would need to make the PSF
### This means that you can run the initial photometry in batch mode
### Then make all the PSFs
### Then run allframe
### Then run the calibration scripts



import sys
import glob
import shutil
import re
import pexpect
import os


import sgr_setup
import sgr_initial_phot



target_name = sys.argv[1]
channel = sys.argv[2]

new_chan = sgr_setup.setup_dir_structure(target_name, channel)

if new_chan == '3p6um': num_chan = 1
if new_chan == '4p5um': num_chan = 2

## Convert the images to counts and give them sensible names
## Also cleans up old runs
target_stem = sgr_setup.sgr_setup(target_name, new_chan)

#target_stem = target_name

file_list = glob.glob(target_stem + '_e*_' + new_chan + '_dn.fits')

## Do initial aperture photometry on all the science images

for image in file_list:
# 	print image
	sgr_initial_phot.init_aper_phot(image, num_chan)
	
print 'Now go and make some psfs!'





















