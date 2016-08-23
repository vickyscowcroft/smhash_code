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



## Grab the input mch file
## First line should be the calibrated epoch 1 data
## Rest should be the alf files
## First file should have the .apc extension
## All others should have .alf

	input_mch = input
	output_mch = re.sub('.mch', '_cal.mch', input)
	output = open(output_mch, 'w')
	
	for line in open(input_mch, 'r'):
		is_epoch1 = re.search('_e1_', line)
		if (is_epoch1 != None):
			splitline = line.split()
			name = splitline[0]
			new_name = name.split(".")[0]+ '.apcor'
			output.write( "{0:s} ' {1:10s} {2:10s} {3:10s} {4:10s} {5:10s} {6:10s} {7:10s} {8:10s} \n".format(new_name, splitline[2], splitline[3], splitline[4], splitline[5], splitline[6], splitline[7], splitline[8], splitline[9]))
		else:
			splitline = line.split()
			name = splitline[0]
			new_name = name.split(".")[0] + '.alf'
			output.write( "{0:s} ' {1:10s} {2:10s} {3:10s} {4:10s} {5:10s} {6:10s} {7:10s} {8:10s} \n".format(new_name, splitline[2], splitline[3], splitline[4], splitline[5], splitline[6], splitline[7], splitline[8], splitline[9]))
		
	output.close()
	
	file_list = glob.glob('*.alf')

	target = re.sub(".mch","", output_mch)
	if (os.path.isfile(target + '.raw')): os.remove(target + '.raw')
	print "running daomaster"
	num_frames = len(file_list)
	print num_frames

	daomaster = pexpect.spawn("daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline(output_mch)
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
	daomaster.sendline(output_mch + "_new")
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
	shutil.move( str(new_mch[0]), output_mch)

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
	sorted_epoch1 = np.argsort(epoch1)

	object_sample = np.zeros(num_stars)
	epoch1_sample = np.zeros(num_stars)
	difference_sample = np.zeros(num_stars)
	
	epoch1_sample = epoch1[sorted_epoch1]
	object_sample = objects[ : , epoch + 3][sorted_epoch1]

	offsets = np.zeros(num_frames)
	sdev_offsets = np.zeros(num_frames)
	howmany_stars = np.zeros(num_frames)
	offsets_log = re.sub('.mch', '.off_log', output_mch)
	offsets_log_file = open(offsets_log, 'w')

		 
	for epoch in np.arange(2, (num_frames*2), 2, dtype=np.int): ## starting from 1 because don't match epoch 1 to itself
		object_sample = objects[ : , epoch + 3][sorted_epoch1]
		difference_sample = epoch1_sample - object_sample
		mask = np.ma.masked_all(num_stars)
		
		#difference = epoch1 - objects[ : , epoch + 3]
		#ediff = np.sqrt(epoch1_err**2 + objects[ : , epoch + 4]**2)
		#clipped = sigma_clip(difference, sigma = 3., iters=None)

		for mag_bin in np.arange(0,11):
			mask[(mag_bin)*num_stars*0.1:(mag_bin+1)*num_stars*0.1] = difference_sample[(mag_bin)*num_stars*0.1:(mag_bin+1)*num_stars*0.1] 
    ## Sigma clip the sample
			clipped3 = sigma_clip(mask, sigma = 3., iters=None, cenfunc=np.ma.mean)
    ## Copy the sigma clipped sample to the masked array
			mask = clipped3
    
		clipped50 = sigma_clip(mask[:num_stars*0.5], sigma = 3., iters=None, cenfunc=np.ma.mean)
		av_diff50 = np.ma.mean(clipped50)
		sdev_diff50 = np.ma.mean(clipped50)

		av_diff = np.ma.mean(clipped50)
		sdev_diff = np.ma.std(clipped50)
		howmany_stars[(epoch/2.)] = clipped50.count()
		offsets[(epoch/2.)] = av_diff
		sdev_offsets[(epoch/2.)] = sdev_diff

	
		#print "Epoch " + str((epoch/2.)+1) + " offset " +str(av_diff) + " sdev " + str(sdev_diff)
		mp.close('all')		
		axp1 = mp.subplot(111)
		#axp1.errorbar(epoch1, clipped, yerr = ediff, color='grey', ls='none')
		axp1.plot(epoch1_sample[:num_stars*0.5], clipped50, 'k.', ls='none')	
		axp1.axhline(av_diff, color='r', ls='--')
		axp1.axhline(av_diff+2*sdev_diff, color='b', ls='--')
		axp1.axhline(av_diff-2*sdev_diff, color='b', ls='--')
		mp.title("Epoch " + str((epoch/2.)+1) + " offset " +str(av_diff) + " sdev " + str(sdev_diff))
		
		#mp.show()
		mp.savefig(target  + str(int((epoch/2.)+1)) + '.pdf')
## Read the file names from the mch file to make sure the order matches up with the order in the raw file

	names = []
	for line in open(output_mch, 'r'):
		splitline = line.split()
		is_epoch1 = re.search('_e1_', splitline[0])
		## Do not apply an offset for epoch 1 
		## So don't append it to the list of files to work on
		if (is_epoch1==None):
			names.append(splitline[0].split("\'")[1])
		
	epoch_count = 0
	for count in np.arange(0, len(names)):
		mtr_name = re.sub(".alf", ".mtr", str(names[count]))
		off_name = re.sub(".alf", ".off", str(names[count]))
		
		input_mtr = open(mtr_name, "r")
		output_mtr = open(off_name, "w")
		for count in np.arange(0,3): 
			header = input_mtr.readline()
			output_mtr.write(header)
		input_mtr.close()
		output_mtr.close()
		print mtr_name, off_name
		id, xc, yc, mag, err = np.loadtxt(mtr_name, skiprows=3, usecols=(0, 1, 2, 3, 4), unpack='TRUE')
		newmag = mag + offsets[epoch_count]
		np.savetxt(off_name, np.column_stack((id, xc, yc, newmag, err)), fmt= "%d %.2f %.2f %.3f %.3f")				
		print count, names[count], off_name
		epoch_count = epoch_count + 1
	
	print "Finished offset calibration"
	
	std_of_offs = np.ma.std(offsets[1:])
	mean_of_offs = np.ma.mean(offsets[1:])
	print std_of_offs, mean_of_offs
	
	### Now write an updated mch file with the correct file names
	
	final_mch = open('tempmch', 'w')
	
	for frame in np.arange(num_frames):
		offsets_log_file.write("{0:d} {1:8.4f} {2:8.4f} {3:d} \n".format(int(frame), offsets[frame], sdev_offsets[frame], int(howmany_stars[frame])))
		if ((abs(offsets[frame] - mean_of_offs) > 3*std_of_offs) and (frame!= 0)):
			print 'OUTLIER:',  int(frame), offsets[frame], sdev_offsets[frame], int(howmany_stars[frame])
		else: print int(frame), offsets[frame], sdev_offsets[frame], int(howmany_stars[frame])
		
		
	offsets_log_file.close()
	
	for line in open(output_mch, 'r'):
		is_epoch1 = re.search('_e1_', line)
		if (is_epoch1 != None):
			final_mch.write( "{0:s}".format(line))
		else:
			splitline = line.split()
			name = splitline[0]
			new_name = name.split(".")[0] + '.off'
			final_mch.write( " {0:s} ' {1:10s} {2:10s} {3:10s} {4:10s} {5:10s} {6:10s} {7:10s} {8:10s} \n".format(new_name, splitline[2], splitline[3], splitline[4], splitline[5], splitline[6], splitline[7], splitline[8], splitline[9]))
		
	final_mch.close()
	
	shutil.move('tempmch', output_mch)
	
	return(0)




	

		




