#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import re
import sys
import glob
import shutil
import pexpect

## Usage: input mch = sys.argv[1]

def location_corr(input):

	stem = re.sub(".mch", "", input)
	mch_file = open(input, "r")
	offs = glob.glob('*_e*_dn.off')
	num_frames = len(offs)
	names = []
	for count in np.arange(0, num_frames):
		line = mch_file.readline()
		lines = line.split("\'")[1]
		lines = re.sub(" ", "",lines)
		names.append(lines)
		print lines

## names currently contains a list of .alf files

	num_files = len(names)
	for image in np.arange(0, num_files):
		off_file = re.sub(".alf", ".off", names[image])
		epoch_s = re.search("_e", off_file)
		#chan = off_file[-12:-7]
		chan_s = re.search("_.p.um", off_file)
		epoch = off_file[epoch_s.start():-12]
		chan = off_file[chan_s.start():-7]
		new_stem = re.sub(chan, "", stem)
	
#		corImage = str(glob.glob("*" + stem + "*" + epoch + "correction_" + chan +".fits")[0])
		corImage = new_stem  + epoch +'correction' + chan + '.fits'
		#corImage = re.sub("\'", "", corImage)
		#corImage
		print corImage
		corFits = fits.open(corImage, mode="readonly")
		corData = corFits[0].data

		newmag = []
		corr = []

		id, xc, yc, mag, err = np.loadtxt(off_file, usecols=(0, 1, 2, 3, 4), unpack='TRUE')
		flux = 10**(mag / -2.5)
		for star in id:
			count = star - 1.
			xcoord = int(np.floor(xc[count]) - 1.)
			ycoord = int(np.floor(yc[count]) - 1.)
			correction = corData[ycoord, xcoord]
			corr.append(correction)

		corr = np.array(corr)

		flux = flux * corr
		newmag = -2.5 * np.log10(flux)


		output_name = re.sub('.off', '.cal',off_file)
		np.savetxt(output_name, np.column_stack((id[np.isnan(corr) == False], xc[np.isnan(corr) == False], yc[np.isnan(corr) == False], newmag[np.isnan(corr) == False], err[np.isnan(corr) == False], corr[np.isnan(corr) == False])), fmt= "%d %.2f %.2f %.3f %.3f %.3f")
		

## Replacing the old extensions with the alf extensions
	mch_file.seek(0)
	f2 = open('temp', 'w')
	for line in mch_file:
	    f2.write(line.replace('.alf', '.cal'))
	f2.close()

	mch_name = input
	cal_file = re.sub(".mch", ".cal", mch_name)

	shutil.move("temp", mch_name)

	daomaster = pexpect.spawn("/Users/vs/daophot/daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline(mch_name)
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
	shutil.move( str(new_cal[0]),cal_file)



