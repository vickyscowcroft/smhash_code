#!/usr/bin/env python

import numpy as np
import re
import sys
import pexpect ## Use this to control daophot
import shutil
import os
import glob
import aplpy

def apcor_photo(image):

	## image in this case should be the flux image
	isdn = re.search('_dn', image)
	if (isdn != None):
		print "you need to run this on the flux image, not the counts image"

## Do the aperture photometry on the flux image in the standard aperture


	if(os.path.isfile("apcor_photo.opt")): os.remove("apcor_photo.opt")
	shutil.copy("/Users/vs/Dropbox/Python/smhash_code/apcor_photo.opt", "apcor_photo.opt")

	if(os.path.isfile(image + '.apc')): os.remove(image + '.apc')
	
	with open(image + '_dn.lst', 'r') as orig:
		header = []
		for count in range(0,2): 
			line = orig.readline()
			header.append(line)
	## Setting lobad data to 0.0 on the flux image so that it can get the photometry done
	photdata = np.loadtxt(image + '_dn.lst', skiprows=3)
	outfile = open(image + '_lst.lst', 'w')
	splitline = header[1].split()
	splitline[3] = 0.0
	outfile.write("{10:s} {0:>2d}{1:>6d}{2:>6d}{3:>8.1f}{4:>8.0f}{5:>8.2f}{6:>8.2f}{7:>8.2f}{8:>8.2f}{9:>8.2f} \n".format(int(splitline[0]), int(splitline[1]), int(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float(splitline[6]), float(splitline[7]), float(splitline[8]), float(splitline[9]), header[0]))
	outfile.write('\n')
	np.savetxt(outfile, photdata, fmt='%d %.2f %.2f %.3f %.3f %.2f')
	outfile.close()


## doing the same as above to the allframe output so i can do substar for a clean aperture correction image

	with open(image + '_dn.alf', 'r') as orig:
		header = []
		for count in range(0,2): 
			line = orig.readline()
			header.append(line)
	## Setting lobad data to 0.0 on the flux image so that it can get the photometry done
	photdata = np.loadtxt(image + '_dn.alf', skiprows=3)
	outfile = open(image + '_alf.alf', 'w')
	splitline = header[1].split()
	splitline[3] = 0.0
	outfile.write("{10:s} {0:>2d}{1:>6d}{2:>6d}{3:>8.1f}{4:>8.0f}{5:>8.2f}{6:>8.2f}{7:>8.2f}{8:>8.2f}{9:>8.2f} \n".format(int(splitline[0]), int(splitline[1]), int(splitline[2]), float(splitline[3]), float(splitline[4]), float(splitline[5]), float(splitline[6]), float(splitline[7]), float(splitline[8]), float(splitline[9]), header[0]))
	outfile.write('\n')
	np.savetxt(outfile, photdata, fmt='  %d %.3f %.3f %.3f %.3f %.2f %0.f %.2f %.3f')
	outfile.close()

	### now match the allframe and lst edited files to make sure that the id numbers match
	
	alf_mch = open('alf.mch', 'w')

	alf_mch.write("'{0:s}_alf.alf'  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))
	alf_mch.write("'{0:s}_lst.lst '  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))

	alf_mch.close()
	for file in glob.iglob('alf_new.mch*'):
		if(os.path.isfile(file)): os.remove(file)
		
	daomaster = pexpect.spawn("daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline('alf.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames")
	daomaster.sendline("2, 1, 2")
	print "2, 1, 2"
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
## Only need the raw here
	print "making output files"
	daomaster.expect("Assign new star IDs")
	daomaster.sendline("n")
	daomaster.expect("A file with mean magnitudes and scatter")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors")
	daomaster.sendline("n")
	daomaster.expect("A file with the new transformations")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline('alf_new.mch')
	daomaster.expect("A file with the transfer table")
	daomaster.sendline("n")

	for file in glob.iglob(image + '*.mtr*'):
		if(os.path.isfile(file)): os.remove(file)
			
	daomaster.expect('Individual .COO files')
	daomaster.sendline('n')
	daomaster.expect("Simply transfer star IDs")
	daomaster.sendline("y")

	daomaster.expect("Good bye")
	daomaster.close(force=True)

	new_mch  = glob.glob('alf_new.mch*')
	shutil.move( str(new_mch[0]), 'alf.mch')

	mtr_file = glob.glob(image + '_alf.mtr*')[0]
	os.remove(mtr_file)
	mtr_file = glob.glob(image + '_lst.mtr*')[0]
	shutil.move(mtr_file, image + '.lst')

## substar
	daophot = pexpect.spawn("daophot")
	daophot.logfile = sys.stdout
	daophot.expect("Command:")
	daophot.sendline("at " + image)
	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("apcor_photo.opt")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(image + '_alf.alf')
	daophot.expect("Output file")
	daophot.sendline(image + '.ap')
	
	if os.path.isfile(image + '.psf'): os.remove(image + '.psf')
	
	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline(image + '.ap')
	daophot.expect("File with PSF stars")
	daophot.sendline(image + '.lst')
	daophot.expect("File for the PSF")
	daophot.sendline(image + '.psf')
	
	
	daophot.expect("Command:")
	daophot.sendline('exit')
	daophot.close(force=True)

	allstar = pexpect.spawn("allstar")
	allstar.logfile = sys.stdout
	allstar.expect("OPT>")
	allstar.sendline("")
	allstar.expect("Input image name")
	allstar.sendline(image)
	allstar.expect("File with the PSF")
	allstar.sendline(image + ".psf")
	allstar.expect("Input file")
	allstar.sendline(image +".ap")
	allstar.expect("File for results")
	allstar.sendline(image +".als")
	allstar.expect("Name for subtracted image")
	allstar.sendline(image +"s.fits")
	allstar.expect("stars")
	allstar.expect("Good bye")
	allstar.close()
	print allstar.isalive()
	print "ALLSTAR complete"
	
	
	daophot = pexpect.spawn("daophot")
	daophot.expect("Command:")
	daophot.sendline('at ' + image)
	daophot.sendline('substar')
 	daophot.expect("File with the PSF")
 	daophot.sendline(image + '.psf')
 	daophot.expect('File with photometry')
 	daophot.sendline(image + '.als')
 	daophot.expect('Do you have stars to leave in?')
 	daophot.sendline('y')
 	daophot.expect('File with star list')
 	daophot.sendline(image + '.lst')
 	daophot.expect('Name for subtracted image')
 	daophot.sendline(image + 's')
 	
 	daophot.expect("Command:")
 	daophot.sendline('at ' + image + 's')

	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("apcor_photo.opt")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(image + '.lst')
	daophot.expect("Output file")
	daophot.sendline(image + '.apc')

	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

## Match up the aperture phot file with the psf phot file
# faking the initial file as they are the same image

	apcor_mch = open('apcor.mch', 'w')

	apcor_mch.write("'{0:s}.apc '  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))
	apcor_mch.write("'{0:s}_dn.alf '  0.0000    0.0000 1.000000000 0.000000000 0.000000000 1.000000000    0.000  0.000\n".format(image))

	apcor_mch.close()

## get the raw file from daomaster
	if(os.path.isfile('apcor.raw')): os.remove('apcor.raw')


	daomaster = pexpect.spawn("daomaster")
	daomaster.expect("File with list of input files")
	daomaster.sendline('apcor.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames")
	daomaster.sendline("2, 1, 2")
	print "2, 1, 2"
	daomaster.expect("Maximum sigma")
	daomaster.sendline("0.2")
## desired degrees of freedom:
	daomaster.expect("Your choice")
	daomaster.sendline("6")
	daomaster.expect("Critical match-up radius")
	daomaster.sendline("10")

	for radius in range (10,-1, -1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(radius))
## Only need the raw here
	print "making output files"
	daomaster.expect("Assign new star IDs")
	daomaster.sendline("n")
	daomaster.expect("A file with mean magnitudes and scatter")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline("apcor.raw")
	daomaster.expect("A file with the new transformations")
	daomaster.sendline("e")
	daomaster.close(force=True)

	new_raw  = glob.glob('apcor.raw*')
	print new_raw
	shutil.move( str(new_raw[0]), "apcor.raw")
	
	## do a bit of clean up
	
	extra_alf = glob.glob('*_alf.alf')[0]
	os.remove(extra_alf)
	
	




