#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp
import glob
import os
import re
import pexpect
import shutil
import itertools
from astropy.stats import sigma_clip

def calibrate(input):

	## Replacing the old extensions with the alf extensions
	f1 = open(input, 'r')
	f2 = open('temp', 'w')
	for line in f1:
	    f2.write(line.replace('.als', '.alf'))
	f1.close()
	f2.close()

	f1 = open('temp', 'r')
	f2 = open(input, 'w')
	for line in f1:
	    f2.write(line.replace('.ap', '.alf'))
	f1.close()
	f2.close()

	os.remove('temp')

	file_list = glob.glob('*.alf')

	target = re.sub(".mch","", input)
	if (os.path.isfile(target + '.raw')): os.remove(target + '.raw')
	print "running daomaster"
	num_frames = len(file_list)
	print num_frames

	daomaster = pexpect.spawn("daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline(input)
	daomaster.expect("Minimum number, minimum fraction, enough frames")
	daomaster.sendline(str(num_frames) + ", 1, " + str(num_frames))
	print str(num_frames) + ", 1, " + str(num_frames)
	daomaster.expect("Maximum sigma")
	daomaster.sendline("99")
	## desired degrees of freedom:
	daomaster.expect("Your choice")
	daomaster.sendline("6")
	daomaster.expect("Critical match-up radius")
	daomaster.sendline("10")

	for radius in range (10,-1, -1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(radius))

## Only need the transformations here
	print "making output files"
	daomaster.expect("Assign new star IDs")
	print '1'
	daomaster.sendline("y")
	daomaster.expect("A file with mean magnitudes and scatter")
	print '2'
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors")
	print '3'
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors")
	print '4'
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(target + ".raw")

	daomaster.expect("A file with the new transformations")
	daomaster.sendline("y")
	print "asked about transformations"
	daomaster.expect("Output file name")
	daomaster.sendline(input + "_new")
	print "asked about output file"
	daomaster.expect("A file with the transfer table")
	daomaster.sendline("n")
	print "asked about transfer table"
	daomaster.expect("Individual .COO files")
	daomaster.sendline("n")

	for alf in file_list:
		mtr = re.sub(".alf", ".mtr",alf)
		old_mtr = glob.glob(mtr +'*')
		if len(old_mtr) > 0:
			old_mtr = glob.glob(mtr +'*')[0]
			if (os.path.isfile(old_mtr)): os.remove(old_mtr)



	daomaster.expect("Simply transfer star IDs")
	daomaster.sendline("y")

	daomaster.expect("Good bye")
	daomaster.close(force=True)

## Fixing the bad output file name
	new_mch  = glob.glob('*.mch_new*')
	print new_mch
	shutil.move( str(new_mch[0]), input)

	new_raw = glob.glob(target +'.raw*')
	print new_raw
	shutil.move( str(new_raw[0]), target + '.raw')

	for alf in file_list:
		mtr = re.sub(".alf", ".mtr",alf)
		new_mtr = glob.glob(mtr +'*')[0]
		shutil.move(new_mtr, mtr)
	
	if ((num_frames % 6 ) == 0):
		n_lines_raw = int((num_frames / 6) + 1)
	else:
		n_lines_raw = int(np.ceil(num_frames) / 6 + 1)
	print n_lines_raw


## Example of how to read in a variable length multi line file!!!
	with open(target + '.raw') as raw_file:
		lines = []
		## skip the 3 header lines
		for count in np.arange(0,3): raw_file.readline()
		while True:
			line = list(itertools.islice(raw_file, n_lines_raw))
			if line:
				lines.append(line)
			else:
				break
		#print lines
		
	num_stars = len(lines)

## Read all relavent data into a numpy array

	objects = np.zeros((num_stars, (num_frames*2) + 3))
	for star in np.arange(0,num_stars):
		data = "".join(lines[star])
		data = data.split()
		objects[star][0] = data[0] ## ID
		objects[star][1] = data[1] ## XC
		objects[star][2] = data[2] ## YC
		for epoch in np.arange(0,(num_frames*2), 2):
			objects[star][epoch + 3]  = data[epoch + 3] ## mag
			objects[star][epoch + 4] = data[epoch + 4] ## err
			 
	epoch1 = objects[:,3]
	epoch1_err = objects[:,4]

	offsets = np.zeros(num_frames)
	sdev_offsets = np.zeros(num_frames)

		 
	for epoch in np.arange(2, (num_frames*2), 2): ## starting from 1 because don't match epoch 1 to itself
		difference = epoch1 - objects[ : , epoch + 3]
		ediff = np.sqrt(epoch1_err**2 + objects[ : , epoch + 4]**2)
		clipped = sigma_clip(difference, sig = 4., iters=100)

		av_diff = np.ma.mean(clipped)
		sdev_diff = np.ma.std(clipped)
		offsets[(epoch/2.)] = av_diff
		sdev_offsets[(epoch/2.)] = sdev_diff

	
		print "Epoch " + str((epoch/2.)+1) + " offset " +str(av_diff) + " sdev " + str(sdev_diff)
		mp.close('all')		
		axp1 = mp.subplot(111)
		axp1.errorbar(epoch1, clipped, yerr = ediff, color='grey', ls='none')
		axp1.plot(epoch1, clipped, 'k.', ls='none')	
		axp1.axhline(av_diff, color='r', ls='--')
		axp1.axhline(av_diff+2*sdev_diff, color='b', ls='--')
		axp1.axhline(av_diff-2*sdev_diff, color='b', ls='--')
	
		mp.show()

## Read the file names from the mch file to make sure the order matches up with the order in the raw file

	names = []
	mch_file = open(target + '.mch', 'r')
	for count in np.arange(0, num_frames):
		line = mch_file.readline()
		lines = line.split("\'")[1]
		names.append(lines)
		
	epoch_count = 0
	for count in np.arange(0, len(names)):
		
		names[count] = re.sub(" ", "", names[count])
		mtr_name = re.sub(".alf", ".mtr", str(names[count]))
		off_name = re.sub(".alf", ".off", str(names[count]))
		
		input = open(mtr_name, "r")
		output = open(off_name, "w")
		for count in np.arange(0,3): 
			header = 	input.readline()
			output.write(header)
		input.close()
		output.close()
		id, xc, yc, mag, err = np.loadtxt(mtr_name, skiprows=3, usecols=(0, 1, 2, 3, 4), unpack='TRUE')
		newmag = mag + offsets[epoch_count]
		np.savetxt(off_name, np.column_stack((id, xc, yc, newmag, err)), fmt= "%d %.2f %.2f %.3f %.3f")
		epoch_count = epoch_count + 1
	
	print "Finished offset calibration"



	

		




