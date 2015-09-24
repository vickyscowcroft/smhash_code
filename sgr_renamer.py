#!/usr/bin/env python

import numpy as np
import re
import sys
import glob
import os

## Target name  = sys.argv[1]
## This would be the long CSS_J name without any epoch specification

target_name = sys.argv[1]

new_target_stem = re.sub("CSS_","", sys.argv[1])

new_target_stem = new_target_stem[0:7]

## Target stem is correct
#print new_target_stem

## Now need to sort out the epochs

files = glob.glob('*_dn.fits') ## These are the new images created by the converter

no_files = len(files) / 2. ## Divide by 2 for 2 channels

#print no_files
## number of files is correct

ch1 = 12
ch2 = 12

for file in files:
	extra = re.search('xc', file)
	if (extra != None):
		is1 = re.search('3p6um', file)
		if (is1 != None):
			ch1 += 1
			newname = new_target_stem + "_e" + str(ch1) + "_3p6um_dn.fits"
		else:
			ch2 += 1
			newname = new_target_stem + "_e" + str(ch2) + "_4p5um_dn.fits"
	else:
		epoch = re.search("_e", file)
		fits = re.search(".fits", file)
		#print epoch.start(), fits.end()
		newname = new_target_stem + file[epoch.start():fits.end()]
	#print newname
	os.rename(file, newname)

		
			


	



