#!/usr/bin/env python

import numpy as np
import re
import sys
import glob
import shutil
import os
import errno
from astropy.io import fits
from  convert_spitzer_flux_to_dn  import *

## Target name  = sys.argv[1]
## This would be the long CSS_J name without any epoch specification

def sgr_setup(target_name, channel):

### Alter the name for targets with long names
### i.e. CSS, linear etc

	css = re.search('css', str.lower(target_name))
	can_split = re.search('_', target_name)
    
	if (css != None):
		#new_target_stem = re.sub("CSS_","", target_name)
		new_target_stem = target_name[0:11]
	elif (can_split != None):
		new_target_stem = str(target_name.split('_', 1)[0][0]) + '_' + str(target_name.split('_', 1)[1])
	else:
		new_target_stem = target_name
		
	## Clean up old versions

	old_files = [ '.alf', '.apc', '.als', '.coo', '.ap', '.raw', '.nmg', '.tfr', '.mch', '.mtr', '.off']
	for ofn in old_files:
		olds = glob.glob('*' + ofn)
		for files in olds:
			os.remove(files)
		

## Target stem is correct
#print new_target_stem

## Now need to sort out the epochs
## Also need to change the names of the correction images
## Otherwise the extra epochs fuck up the correction stage

	regex = re.compile(re.sub('_', '\_', new_target_stem)) ## creating a regex to find only the images corresponding to this target
	regex_target = re.sub('_', '\_', target_name)
	regex_target = re.sub('\+', '\+', regex_target)
	regex_target = re.sub('\-', '\-', regex_target)

	old_regex = re.compile(regex_target)
	files = list(set(glob.glob('*_e*' + channel + '.fits') + glob.glob('*correction*.fits')))## These are the original images and the correction images
	#files = filter(regex.search, files)	
	original_names = filter(old_regex.search, files)
	#print files

	ch1 = 12
	ch2 = 12

	print 'target name = ', target_name
	print 'new_target_stem = ', new_target_stem
	

	for filename in original_names:
	## Just change file names of correction images
		corr = re.search('correction', filename)
		if (corr != None):
			new_corr_name = re.sub('__e', '_e', filename)
			new_corr_name = re.sub(regex_target, new_target_stem, new_corr_name) 
			shutil.move(filename, new_corr_name)
			continue
		if (corr == None):
			extra = re.search('xc', filename)
			if (extra != None):
				ch1 += 1
				newname = new_target_stem + "_e" + str(ch1) + "_" + channel + "_dn.fits"
				corr_file = re.sub("_" + channel, "_correction_" + channel, filename)
				new_correction_name = new_target_stem + "_e" + str(ch1) + "_correction_" + channel + ".fits"
			else:
				epoch = re.search("_e", filename)
				fit = re.search(".fits", filename)
				is1 = re.search('3p6um', filename)
				if (is1 != None):
					newname = new_target_stem + filename[epoch.start():fit.start()] + "_dn.fits" 
					corr_file = re.sub("_3p6um", "_correction_3p6um", filename)
					new_correction_name = re.sub("_3p6um_dn.fits", "_correction_3p6um.fits", newname)
				if (is1 == None):
					newname = new_target_stem + filename[epoch.start():fit.start()] + "_dn.fits" 
					corr_file = re.sub("_4p5um", "_correction_4p5um", filename)
					new_correction_name = re.sub("_4p5um_dn.fits", "_correction_4p5um.fits", newname)

			print filename, newname
			spitzer_flux2dn(filename, newname)	

	
## Now also convert the mosaic
## target_name is the target without epoch specification. Can use that
## Want to copy the median and the science stack

	channel_list = ['3p6um', '4p5um']
#	channel_list = ['3p6um']
	for channel in channel_list:
		image_name = target_name + "_" + channel + '.fits'
		median_name = target_name + "_" + channel + '_median.fits'
	#print image_name
		if (os.path.isfile(image_name) and os.path.isfile(median_name)):
			spitzer_flux2dn_mos(image_name, median_name, new_target_stem, channel)
		
		elif (os.path.isfile(image_name) and (os.path.isfile(median_name)==False)):	
			new_name = new_target_stem + '_' + channel + '_dn.fits'
			spitzer_flux2dn(image_name, new_name)
		shutil.copy(target_name + '__e1_' + channel + '.fits', new_target_stem + '_e1_' + channel +'.fits')
		
	return(new_target_stem)

def setup_single_file(target_name):

	regex = re.compile(target_name) ## creating a regex to find only the images corresponding to this target
	files = glob.glob('*' + channel + '.fits')
	files = filter(regex.search, files)	
	for filename in files:
		corr = re.search('correction', filename)
		if (corr != None):
			continue
		if (corr == None):
			if (channel == '3p6um'):
				new_name = re.sub('3p6um.fits', '3p6um_dn.fits', filename)
				mapname = re.sub('3p6um.fits', '3p6um_exposure.fits', filename)
			if (channel == '4p5um'):
				new_name = re.sub('4p5um.fits', '4p5um_dn.fits', filename)
				mapname = re.sub('4p5um.fits', '4p5um_exposure.fits', filename)
			flux_to_dn_expmap(filename, mapname, new_name)
	new_target_stem = target_name
	return(new_target_stem)
	
def bad_pixel_mask(file_list): ### Creating bad pixel masks for all the science images
	for file in file_list:
		bp_name = re.sub('_dn.fits', '_bp.fits', file)
		shutil.copy(file, bp_name) ## copied original image to the bad pixel mask image
		hdu_list = fits.open(bp_name, mode='update')
		mask = hdu_list[0].data
		mask[mask=='nan'] = 1
		mask[mask>0] = 0
		mask[mask<=0] = 0
		hdu_list.flush()
		hdu_list.close()
		
def setup_dir_structure(target_name, channel):
	if (channel == 1 or channel == '3p6um' or channel == '1'): channel = '3p6um'
	elif (channel == 2 or channel == '4p5um' or channel == '2'): channel = '4p5um'
	else: 
		print 'invalid channel'
		exit(1)
	curr_dir = os.getcwd().split('/')[-1]
	new_dir = target_name + '_' + channel
	is_there = os.path.exists(new_dir)
	if (curr_dir != new_dir):
		if (is_there == False):
			os.mkdir(new_dir)
		os.chdir(new_dir)
	
	regex_target = re.sub('_', '\_', target_name)
	regex_target = re.sub('\+', '\+', regex_target)
	regex_target = re.sub('\-', '\-', regex_target)
	regex = re.compile(regex_target)

	file_list = glob.glob('../*/*.fits')
	file_list = filter(regex.search, file_list)
	epoch_link_list = [x.split('/', 2)[2] for x in file_list]
	if (len(file_list) == 0):
		file_list = glob.glob('../*.fits')
		file_list = filter(regex.search, file_list)
		epoch_link_list = [x.split('/', 1)[1] for x in file_list]

	
	print file_list
	mosaic_list = glob.glob('../MegaMosaics/*/*.fits')
	mosaic_list = filter(regex.search, mosaic_list)
	mosaic_link_list = [x.split('/', 3)[3] for x in mosaic_list]

	if (len(mosaic_list) == 0):
		mosaic_list = glob.glob('../../MegaMosaics/*.fits')
		mosaic_list = filter(regex.search, mosaic_list)
		mosaic_link_list = [x.split('/', 3)[3] for x in mosaic_list]

	for count in range(len(file_list)):
		#print file_list[count], epoch_link_list[count]
		try:
			os.symlink(file_list[count], epoch_link_list[count])
		except OSError, e:
			if e.errno == errno.EEXIST:
				os.remove(epoch_link_list[count])
				os.symlink(file_list[count], epoch_link_list[count])
	for count in range(len(mosaic_list)):
		#print mosaic_list[count], mosaic_link_list[count]
		try:
			os.symlink(mosaic_list[count], mosaic_link_list[count])
		except OSError, e:
			if e.errno == errno.EEXIST:
				os.remove(mosaic_link_list[count])
				os.symlink(mosaic_list[count], mosaic_link_list[count])	
	print 'Directory structure set up'	
	return(channel)
		
def get_target_stem(target_name):
	
	css = re.search('css', str.lower(target_name))
	can_split = re.search('_', target_name)
    
	if (css != None):
		#new_target_stem = re.sub("CSS_","", target_name)
		new_target_stem = target_name[0:11]
	elif (can_split != None):
		new_target_stem = str(target_name.split('_', 1)[0][0]) + '_' + str(target_name.split('_', 1)[1])
	else:
		new_target_stem = target_name
	return(new_target_stem)

		
		
		
		
		
	
	