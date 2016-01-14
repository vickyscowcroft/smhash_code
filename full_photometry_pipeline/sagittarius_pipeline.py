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
location = sys.argv[2] ## location is the location of the RRL - can currently be orphan, sgr_stream

if(location=='orphan'):
	cat_target = target_name
	cat_name = '/Users/vs/Dropbox/SMHASH/orphan_rrlyrae_catalogue'
if(location=='sgr_stream'):
	cat_target = target_name[0:11] + '.' + target_name[12:20]
	cat_name = '/Users/vs/Dropbox/SMHASH/catalina_sagittarius_rrlyrae_catalogue'

## Convert the images to counts and give them sensible names
## Also cleans up old runs
target_stem = sgr_setup.sgr_setup(target_name)

file_list = glob.glob(target_stem + '_e*_3p6um_dn.fits')

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

if (psfdone == 'continue'):
	print "Continuing..."

else:
	print "Have you made a good PSF?"
	psfdone = raw_input('Type continue when PSF is complete: ')	
	
## Now copy the psf to each epoch
shutil.copy(target_stem + '_e1_3p6um_dn.psf', 'master_3p6um.psf')
for file in file_list:
	psfname = re.sub(".fits",".psf", file)
	shutil.copy("master_3p6um.psf", psfname)
	
## Prepare the files for allframe

sgr_pre_allframe_matching.sgr_pre_allframe_matching(target_stem, '1')
print 'completed pre allframe matching?'

## Do the photometry on the median image and calculate it's offset from epoch 1

median_name = target_stem + '_3p6um_med_dn.fits'
print 'median name is ', median_name

sgr_median_image_phot.sgr_median_image_phot(median_name)

## ALLFRAEME

## Remove old files from possible previous allframe runs to prevent restarting crashed run
if os.path.isfile(target_stem + '_3p6m.bck'): os.remove(target_stem + '_3p6um.bck')
if os.path.isfile(target_stem + '_3p6um.tfr'): os.remove(target_stem + '_3p6um.tfr')

print "running allframe"

allframe = pexpect.spawn("/Users/vs/daophot/allframe")
allframe.logfile = sys.stdout
allframe.expect("OPT>")
allframe.sendline("")
allframe.expect("File with list of images")
allframe.sendline(target_stem + '_3p6um.mch')
allframe.expect("File with list of stars")
allframe.sendline(target_stem + '_3p6um.mag')
allframe.expect("Starting image")
allframe.expect("Good bye", timeout=6000)
allframe.close()
print allframe.isalive()

print "ALLFRAME complete"

## Begin calibration process

## Aperture correction

flux_epoch_1 = target_stem + '_e1_3p6um'

aperture_correction_photometry.apcor_photo(flux_epoch_1)

apcor, sdev_apcor, period = calculate_aperture_correction.calc_apcor(flux_epoch_1, 'apcor.raw', 2.0, cat_target, cat_name)

apply_aperture_correction.apply_apcor(target_stem + '_e1_3p6um_dn.alf', apcor, 1)

## Calibrate all other epochs
calibrate_all_epochs.calibrate(target_stem + '_3p6um.mch')

## location correction
all_location_corrections.location_corr(target_stem + '_3p6um.mch')

print "Calibration complete"

print "Individual epochs are in " + target_stem + "_epoch_3p6um_dn.cal files"
print "Matched file is " + target_stem + '_3p6um.cal'




















