#!/usr/bin/env python

### Last edit November 8 2016
## Runs pre allframe matching and allframe only
## Run this only after you have made the master PSF



import sys
import shutil
import pexpect
import os
import glob
import re

import sgr_setup
import sgr_pre_allframe_matching
import traccs_master_image_phot


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

print target_stem

### Checking that you've made a good psf for epoch 1

if((os.path.isfile(target_stem + '_e1_' + new_chan + '_dn.psf') == False) and  os.path.isfile('master_' + new_chan + '.psf')==False):
	print 'You need to make the master psf first!!'
	exit(1)
	
	
## Now copy the psf to each epoch
if (os.path.isfile('master_' + new_chan + '.psf')==False):
	shutil.copy(target_stem + '_e1_' + new_chan + '_dn.psf', 'master_' + new_chan + '.psf')
	
file_list = glob.glob(target_stem + '_e*_' + new_chan + '_dn.fits')
	
	
for file in file_list:
	psfname = re.sub(".fits",".psf", file)
	shutil.copy("master_" + new_chan + ".psf", psfname)
	
## Prepare the files for allframe

print target_stem

sgr_pre_allframe_matching.sgr_pre_allframe_matching(target_stem, num_chan)
print 'completed pre allframe matching?'

## Do the photometry on the master image and calculate it's offset from epoch 1

master_name = target_stem + '_' + new_chan + '_dn.fits'
print 'master name is ', master_name

traccs_master_image_phot.traccs_master_image_phot(master_name, True)

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

