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
import calculate_aperture_correction_testing
import apply_aperture_correction
import calibrate_all_epochs
import all_location_corrections
import os


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
	
## Haven't automated this part yet. 
## You have to go to IRAF, select good PSF stars, make a daophot PSF, check the allstar subtraction and come back to this window
## Keep the name of the .psf file as the default
print "Now make a good PSF for epoch 1 by hand"
print "Select stars visually, choosing them from the .coo file and copying their info to the .lst file"
print "Check the subtracted image for epoch 1"
print "IMPORTANT!! If you've copied a psf from somewhere, you still need to run allstar on epoch 1"
print "VS will fix this `feature' in a later version"
print "When you're done, type continue"
# 
psfdone = raw_input('Type continue when PSF is complete: ')
# 
if (psfdone == 'continue'):
	print "Continuing..."
# 
else:
	print "Have you made a good PSF?"
	psfdone = raw_input('Type continue when PSF is complete: ')	
	
## Now copy the psf to each epoch
shutil.copy(target_stem + '_e1_' + new_chan + '_dn.psf', 'master_' + new_chan + '.psf')
for file in file_list:
	psfname = re.sub(".fits",".psf", file)
	shutil.copy("master_" + new_chan + ".psf", psfname)
	
## Prepare the files for allframe

sgr_pre_allframe_matching.sgr_pre_allframe_matching(target_stem, num_chan)
print 'completed pre allframe matching?'

## Do the photometry on the median image and calculate it's offset from epoch 1

median_name = target_stem + '_' + new_chan + '_med_dn.fits'
print 'median name is ', median_name

sgr_median_image_phot.sgr_median_image_phot(median_name, True)

## ALLFRAEME

## Remove old files from possible previous allframe runs to prevent restarting crashed run
if os.path.isfile(target_stem + '_' + new_chan + '.bck'): os.remove(target_stem + '_' + new_chan + '.bck')
if os.path.isfile(target_stem + '_' + new_chan + '.tfr'): os.remove(target_stem + '_' + new_chan + '.tfr')
# 
print "running allframe"
# 
allframe = pexpect.spawn("allframe")
allframe.logfile = sys.stdout
allframe.expect("OPT>")
allframe.sendline("")
allframe.expect("File with list of images")
allframe.sendline(target_stem + '_' + new_chan + '.mch')
allframe.expect("File with list of stars")
allframe.sendline(target_stem + '_' + new_chan + '.mag')
allframe.expect("Starting image")
allframe.expect("Good bye", timeout=6000)
allframe.close()
print allframe.isalive()

print "ALLFRAME complete"

## Begin calibration process

## Aperture correction

flux_epoch_1 = target_stem + '_e1_' + new_chan

aperture_correction_photometry.apcor_photo(flux_epoch_1)

cat_name = 'globular'

## The 2.0 on the input here refers to the sigma clipping, not the channel
apcor, sdev_apcor = calculate_aperture_correction_testing.calc_apcor(flux_epoch_1, 'apcor.raw', 2.0, target_stem)

apply_aperture_correction.apply_apcor(target_stem + '_e1_' + new_chan + '_dn.alf', apcor, num_chan)

## Calibrate all other epochs
calibrate_all_epochs.calibrate(target_stem + '_' + new_chan + '.mch')

## location correction
all_location_corrections.location_corr(target_stem + '_' + new_chan + '.mch')

print "Calibration complete"

print "Individual epochs are in " + target_stem + "_epoch_' + new_chan + '_dn.cal files"
print "Matched file is " + target_stem + '_' + new_chan + '.cal'




















