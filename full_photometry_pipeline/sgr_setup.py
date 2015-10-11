#!/usr/bin/env python

import numpy as np
import re
import sys
import glob
import shutil
import os
from  convert_spitzer_flux_to_dn  import *

## Target name  = sys.argv[1]
## This would be the long CSS_J name without any epoch specification

def sgr_setup(target_name):

	new_target_stem = re.sub("CSS_","", target_name)

	new_target_stem = new_target_stem[0:7]
	
	## Clean up old versions
	old_files = glob.glob(new_target_stem + '*')
	for ofn in old_files:
		os.remove(ofn)
		

## Target stem is correct
#print new_target_stem

## Now need to sort out the epochs
## Also need to change the names of the correction images
## Otherwise the extra epochs fuck up the correction stage

	regex = re.compile(new_target_stem) ## creating a regex to find only the images corresponding to this target
	files = glob.glob('CSS*_e*3p6um.fits') ## These are the original images
	files = filter(regex.search, files)	
	#print files

	ch1 = 12
	ch2 = 12

	for filename in files:
	## Just change file names of correction images
		corr = re.search('correction', filename)
		if (corr != None):
			continue
		if (corr == None):
			extra = re.search('xc', filename)
			if (extra != None):
				ch1 += 1
				newname = new_target_stem + "_e" + str(ch1) + "_3p6um_dn.fits"
				corr_file = re.sub("_3p6um", "_correction_3p6um", filename)
				new_correction_name = new_target_stem + "_e" + str(ch1) + "_correction_3p6um.fits"
			else:
				epoch = re.search("_e", filename)
				fit = re.search(".fits", filename)
				is1 = re.search('3p6um', filename)
				if (is1 != None):
					newname = new_target_stem + filename[epoch.start():fit.start()] + "_dn.fits" 
					corr_file = re.sub("_3p6um", "_correction_3p6um", filename)

					new_correction_name = re.sub("_3p6um_dn.fits", "_correction_3p6um.fits", newname)
		print filename, newname
		spitzer_flux2dn(filename, newname)
		print "copying the correction file"
		shutil.copy(corr_file, new_correction_name)
			

	
## Now also convert the mosaic
## target_name is the target without epoch specification. Can use that
## Want to copy the median and the science stack

#	channel_list = ['3p6um', '4p5um']
	channel_list = ['3p6um']
	for channel in channel_list:
		image_name = target_name + "_" + channel + '.fits'
		median_name = target_name + "_" + channel + '_median.fits'
	#print image_name
		if (os.path.isfile(image_name) and os.path.isfile(median_name)):
			spitzer_flux2dn_mos(image_name, median_name, new_target_stem, channel)
		shutil.copy(target_name + '__e1_' + channel + '.fits', new_target_stem + '_e1_' + channel +'.fits')
		
		
	return(new_target_stem)

