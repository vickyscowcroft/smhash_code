#!/usr/bin/env python

import numpy as np
import re
import sys
import glob
import shutil
from astropy.io import fits
import os
from  convert_spitzer_flux_to_dn  import *

## Target name  = sys.argv[1]
## This would be the long CSS_J name without any epoch specification


target_name = sys.argv[1]

new_target_stem = re.sub("CSS_","", sys.argv[1])

new_target_stem = new_target_stem[0:7]

## Target stem is correct
#print new_target_stem

## Now need to sort out the epochs

regex = re.compile(new_target_stem) ## creating a regex to find only the images corresponding to this target
files = glob.glob('*_e*um.fits') ## These are the original images
files = filter(regex.search, files)
print files


ch1 = 12
ch2 = 12

for filename in files:
	## Do not do to correction images!!
	corr = re.search('correction', filename)
	if (corr != None):
		continue
	extra = re.search('xc', filename)
	if (extra != None):
		is1 = re.search('3p6um', filename)
		if (is1 != None):
			ch1 += 1
			newname = new_target_stem + "_e" + str(ch1) + "_3p6um_dn.fits"
		else:
			ch2 += 1
			newname = new_target_stem + "_e" + str(ch2) + "_4p5um_dn.fits"
	else:
		epoch = re.search("_e", filename)
		fit = re.search(".fits", filename)
		is1 = re.search('3p6um', filename)
		if (is1 != None):
			newname = new_target_stem + filename[epoch.start():fit.start()] + "_dn.fits" 
		else:
			newname = new_target_stem + filename[epoch.start():fit.start()] + "_dn.fits" 
	print newname
	spitzer_flux2dn(filename, newname)

	
## Now also convert the mosaic
## target_name is the target without epoch specification. Can use that
## Want to copy the median and the science stack

channel_list = ['3p6um', '4p5um']

for channel in channel_list:
	image_name = target_name + "_" + channel + '.fits'
	median_name = target_name + "_" + channel + '_median.fits'
	#print image_name
	if (os.path.isfile(image_name) and os.path.isfile(median_name)):
		spitzer_flux2dn_mos(image_name, median_name, new_target_stem, channel)
		


