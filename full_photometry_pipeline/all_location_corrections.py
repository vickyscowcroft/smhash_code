#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import re
import sys
import glob
import shutil
import pexpect

## Usage: input mch = sys.argv[1]

def location_corr(target_stem, new_chan):

	input = target_stem + '_' + new_chan + '_cal.mch'
	names = []

	for line in open(input, 'r'):
		splitline = line.split()
		names.append(splitline[0].split("\'")[1])

		
	### names contains the list of the calibrated photometry files, before location correction
	
	## construct the correction file names 
	## format is TARGET + EPOCH + correction + CHANNEL.FITS
	num_files = len(names)
	for image in np.arange(0, num_files):
		epoch = re.search('_e([0-9]|[0-9][0-9])_', names[image]).group(0)
	
#		corImage = str(glob.glob("*" + stem + "*" + epoch + "correction_" + chan +".fits")[0])
		corImage = target_stem  + epoch + 'correction_' + new_chan + '.fits'
		print corImage
		corFits = fits.open(corImage, mode="readonly")
		corData = corFits[0].data

		newmag = []
		corr = []

		id, xc, yc, mag, err = np.loadtxt(names[image], usecols=(0, 1, 2, 3, 4), unpack='TRUE')
		flux = 10**(mag / -2.5)
		for star in np.arange(len(id)):
			xcoord = int(np.floor(xc[star]) - 1.)
			ycoord = int(np.floor(yc[star]) - 1.)
			correction = corData[ycoord, xcoord]
			corr.append(correction)

		corr = np.array(corr)

		flux = flux * corr
		newmag = -2.5 * np.log10(flux)

		is_epoch1 = re.search('_e1_', names[image])
		if (is_epoch1 != None):
			output_name = re.sub('.apcor', '.cal',names[image])
		else:
			output_name = re.sub('.off', '.cal',names[image])
		np.savetxt(output_name, np.column_stack((id[np.isnan(corr) == False], xc[np.isnan(corr) == False], yc[np.isnan(corr) == False], newmag[np.isnan(corr) == False], err[np.isnan(corr) == False], corr[np.isnan(corr) == False])), fmt= "%d %.2f %.2f %.3f %.3f %.3f")
		

## Replacing the old extensions with the alf extensions
	final_mch = open('tempmch', 'w')
	
	for line in open(input, 'r'):
		splitline = line.split()
		name = splitline[0]
		new_name = name.split(".")[0] + '.cal'
		final_mch.write( " {0:s} ' {1:10s} {2:10s} {3:10s} {4:10s} {5:10s} {6:10s} {7:10s} {8:10s} \n".format(new_name, splitline[2], splitline[3], splitline[4], splitline[5], splitline[6], splitline[7], splitline[8], splitline[9]))
		
	final_mch.close()
	
	shutil.move('tempmch', input)

	mch_name = input
	cal_file = target_stem + '_' + new_chan + '.temp_cal'
	#re.sub(".mch", ".cal", mch_name)


	daomaster = pexpect.spawn("daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline(mch_name)
	daomaster.expect("Minimum number, minimum fraction, enough frames")
	daomaster.sendline(str(num_files) + ", 1, " + str(num_files))
	print str(num_files) + ", 1, " + str(num_files)
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
	daomaster.sendline("n")
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
	daomaster.sendline(cal_file)

	daomaster.expect("A file with the new transformations")
	daomaster.sendline("y")
	print "asked about transformations"
	daomaster.expect("Output file name")
	daomaster.sendline(mch_name + "_new")
	print "asked about output file"
	daomaster.expect("A file with the transfer table")
	daomaster.sendline("e")
	print "asked about transfer table"
	daomaster.expect("Good bye")
	daomaster.close(force=True)


	new_mch  = glob.glob('*.mch_new*')
	print new_mch
	shutil.move( str(new_mch[0]), mch_name)


	new_cal = glob.glob(cal_file + '*')
	print new_cal
	shutil.move(str(new_cal[0]),cal_file)
	
	cal_file2 = target_stem + '_' + new_chan + '.cal'
	shutil.move(cal_file, cal_file2)

	
	
	return(0)

def single_location_correction(target_stem, new_chan):
	input = target_stem + '_' + new_chan + '_cal.mch'
	
	corImage = target_stem + '_correction_' + new_chan + '.fits'
	corFits = fits.open(corImage, mode="readonly")
	corData = corFits[0].data

	newmag = []
	corr = []

	id, xc, yc, mag, err = np.loadtxt(names[image], usecols=(0, 1, 2, 3, 4), unpack='TRUE')
	flux = 10**(mag / -2.5)
	for star in id:
		count = star - 1.
		xcoord = int(np.floor(xc[count]) - 1.)
		ycoord = int(np.floor(yc[count]) - 1.)
		correction = corData[ycoord, xcoord]
		corr.append(correction)
	ext = input[-4:]
	corr = np.array(corr)

	flux = flux * corr
	newmag = -2.5 * np.log10(flux)


	output_name = re.sub(ext, '.cal',input)
	np.savetxt(output_name, np.column_stack((id[np.isnan(corr) == False], xc[np.isnan(corr) == False], yc[np.isnan(corr) == False], newmag[np.isnan(corr) == False], err[np.isnan(corr) == False], corr[np.isnan(corr) == False])), fmt= "%d %.2f %.2f %.3f %.3f %.3f")


