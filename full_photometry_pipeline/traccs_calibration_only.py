#!/usr/bin/env python


import sys
import glob
import sgr_setup
import shutil
import aperture_correction_photometry
import calculate_aperture_correction_testing
import apply_aperture_correction
import calibrate_all_epochs
import all_location_corrections
import os


target_name = sys.argv[1]
channel = sys.argv[2]

if (channel == 1 or channel == '3p6um' or channel == '1'): new_chan = '3p6um'
elif (channel == 2 or channel == '4p5um' or channel == '2'): new_chan = '4p5um'
else: 
	print 'invalid channel'
	exit(1)

if new_chan == '3p6um': num_chan = 1
if new_chan == '4p5um': num_chan = 2
target_stem = sgr_setup.get_target_stem(target_name)


if (len(glob.glob(target_stem + '*' + channel +'*.alf'))==0):
	print 'You need to run ALLFRAME before you can calibrate the photometry'
	print 'This is the CALIBRATION ONLY script'
	print 'Run ``globular_pipeline_anychannel.py`` to do the combined photometry and calibration script'
	exit(1)

## Convert the images to counts and give them sensible names
## Also cleans up old runs

## for omegaCen only - doesn't apply for Sgr stream, orphan stream etc:

## Aperture correction

flux_epoch_1 = target_stem + '_e1_' + new_chan

aperture_correction_photometry.apcor_photo(flux_epoch_1)

cat_name = 'globular'

## The 2.0 on the input here refers to the sigma clipping, not the channel
apcor, sdev_apcor = calculate_aperture_correction_testing.calc_apcor(flux_epoch_1, 'apcor.raw', 2.5, target_stem)

apply_aperture_correction.apply_apcor(target_stem + '_e1_' + new_chan + '_dn.alf', apcor, num_chan)

## Calibrate all other epochs
calibrate_all_epochs.calibrate(target_stem + '_' + new_chan + '.mch')


### August 18 - 3:12pm - pipeline works to here. Need to fix location correction now some naming conventions have been changed.

## location correction
all_location_corrections.location_corr(target_stem, new_chan)

print "Calibration complete"

print "Individual epochs are in " + target_stem + "_epoch_' + new_chan + '_dn.cal files"
print "Matched file is " + target_stem + '_' + new_chan + '.cal'




















