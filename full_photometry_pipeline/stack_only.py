#!/usr/bin/env python


import sys
import glob
import sgr_setup
import sgr_initial_phot
import sgr_pre_allframe_matching
import sgr_median_image_phot
import shutil
import re
import pexpect
import aperture_correction_photometry
import calculate_aperture_correction
import apply_aperture_correction
import calibrate_all_epochs
import all_location_corrections
import os


target_name = sys.argv[1]

## Convert the images to counts and give them sensible names
## Also cleans up old runs
target_stem = sgr_setup.setup_single_file(target_name)

search_term = target_stem + '_3p6um_dn.fits'
print search_term

file_list = glob.glob(search_term)

## Do initial aperture photometry on all the science images

for image in file_list:
	print image
	sgr_initial_phot.init_aper_phot(image)
	
## Haven't automated this part yet. 
## You have to go to IRAF, select good PSF stars, make a daophot PSF, check the allstar subtraction and come back to this window
## Keep the name of the .psf file as the default
print "Now make a good PSF for epoch 1 by hand"
print "Select stars visually, choosing them from the .coo file and copying their info to the .lst file"
print "Check the subtracted image for epoch 1"
print "IMPORTANT!! If you've copied a psf from somewhere, you still need to run allstar on epoch 1"
print "VS will fix this `feature' in a later version"
print "When you're done, type continue"

psfdone = raw_input('Type continue when PSF is complete: ')


## Begin calibration process

## Aperture correction

flux_epoch_1 = target_stem + '_3p6um'

aperture_correction_photometry.apcor_photo(flux_epoch_1)

apcor, sdev_apcor, period = calculate_aperture_correction.calc_apcor(flux_epoch_1, 'apcor.raw', 2.0, cat_target, cat_name)

apply_aperture_correction.apply_apcor(target_stem + '_3p6um_dn.alf', apcor, 1)

## location correction
all_location_corrections.single_location_correction(target_stem + '_3p6um.apc')

print "Calibration complete"



















