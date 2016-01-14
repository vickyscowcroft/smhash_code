#!/usr/bin/env python

import numpy as np
import re
import sys
import pexpect ## Use this to control daophot
import shutil
from astropy.io import fits
import os
import glob

def sgr_median_image_phot(fitsfile):

    image = re.sub(".fits", "", fitsfile)
    is1 = re.search('3p6um', image)
    if (is1 != None):
        imstem = re.sub("_3p6um_med_dn","", image)
    else:
        imstem = re.sub("_4p5um_med_dn","", image)
    print 'made it past the first hurdle!'

	## Copy the daophot opt file
    shutil.copy("/Users/vs/Dropbox/Python/smhash_code/daophot-spitzer-i1.opt", "daophot.opt")
    shutil.copy("/Users/vs/Dropbox/Python/smhash_code/photo-spitzer.opt", "photo.opt")

    print "Working on " + image
## Clean up previous runs

    extensions = ['.coo', '.lst', '.psf', '.nei', '.ap', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.mag']
    for ext in extensions:
        if (os.path.isfile(image + ext)):
            os.remove(image+ext)

## First daophot run
    daophot = pexpect.spawn("/Users/vs/daophot/daophot")

    hdulist = fits.open(fitsfile)
    prihdr = hdulist[0].header
# attach the image
    daophot.expect("Command:")
    daophot.sendline("at " + image)
    daophot.expect("Command:")
    daophot.sendline("opt")
    daophot.expect("File with parameter values")
    daophot.sendline("")
    daophot.expect("OPT>")
    daophot.sendline("th=10")
    daophot.expect("OPT>")
    daophot.sendline("")

# find the stars
    daophot.expect("Command:")
    daophot.sendline("find")
    daophot.expect("Number of frames averaged, summed:")
    daophot.sendline("1, 10")
    daophot.expect("File for positions")
    daophot.sendline("")
    daophot.expect("Are you happy with this?")
    daophot.sendline("y")

    print "FIND complete"

## Aperture photometry
    daophot.expect("Command:")
    daophot.sendline("phot")
    daophot.expect("File with aperture radii")
    daophot.sendline("")
    daophot.expect("PHO>")
    daophot.sendline("")
    daophot.expect("Input position file")
    daophot.sendline(image + '.coo')
    daophot.expect("Output file")
    daophot.sendline(image + '.ap')

    print "PHOT complete"

## Exit daophot
    daophot.expect("Command:")
    daophot.sendline("exit")
    daophot.close(force=True)

    print "Initial daophot run complete."

    if (is1 != None):
        shutil.copy("master_3p6um.psf", image + ".psf")
    else:
        shutil.copy("master_4p5um.psf", image + ".psf")

## Allstar run 1
    allstar = pexpect.spawn("/Users/vs/daophot/allstar")
    allstar.logfile = sys.stdout
    allstar.expect("OPT>")
    allstar.sendline("")
    allstar.expect("Input image name")
    allstar.sendline(image)
    print '1'
    allstar.expect("File with the PSF")
    allstar.sendline(image + ".psf")
    print '2'

    allstar.expect("Input file")
    allstar.sendline(image +".ap")
    print '3'

    allstar.expect("File for results")
    allstar.sendline(image +".als")
    print '4'

    allstar.expect("Name for subtracted image")
    allstar.sendline(image +"s.fits")
    print '5'

    allstar.expect("stars")
    allstar.expect("Good bye")
    allstar.close()
    print allstar.isalive()
    print "ALLSTAR complete"

## Finding missed stars
    daophot = pexpect.spawn("/Users/vs/daophot/daophot")
    daophot.logfile = sys.stdout

    daophot.expect("Command:")
    daophot.sendline("at " + image +"s")
    daophot.expect("Command:")
    daophot.sendline("opt")
    daophot.expect("File with parameter values")
    daophot.sendline("")
    daophot.expect("OPT>")
    daophot.sendline("th=13")
    daophot.expect("OPT>")
    daophot.sendline("lo=10")
    daophot.expect("OPT>")
    daophot.sendline("")

# find the stars in the subtracted image
    daophot.expect("Command:")
    daophot.sendline("find")
    daophot.expect("Number of frames averaged, summed:")
    daophot.sendline("1, 10")
    daophot.expect("File for positions")
    daophot.sendline("")
    daophot.expect("Are you happy with this?")
    daophot.sendline("y")

## combine and sort the files
    daophot.expect("Command:")
    daophot.sendline("append")
    daophot.expect("First input file")
    daophot.sendline(image +".coo")
    daophot.expect("Second input file")
    daophot.sendline(image + "s.coo")
    daophot.expect("Output file")
    daophot.sendline(image + ".cmb")
    daophot.expect("Command:")
    daophot.sendline("sort")
    daophot.expect("Which do you want") ## Sort by increasing x coord
    daophot.sendline("2")
    daophot.expect("Input file name")
    daophot.sendline(image + ".cmb")
    daophot.expect("Output file name")
    daophot.sendline(image + ".srt")
    daophot.expect("Do you want the stars renumbered")
    daophot.sendline("y")

## aperture photometry on the original image
    daophot.expect("Command:")
    daophot.sendline("at " + image)
    daophot.expect("Command")
    daophot.sendline("phot")
    daophot.expect("File with aperture radii")
    daophot.sendline("")
    daophot.expect("PHO>")
    daophot.sendline("")
    daophot.expect("Profile-fitting photometry")
    daophot.sendline("e")
    daophot.expect("Input position file")
    daophot.sendline(image + '.srt')
    daophot.expect("Output file")
    daophot.sendline(image + 's.ap')

    daophot.expect("Command:")
    daophot.sendline("exit")
    daophot.close(force=True)


## Second allstar run
    allstar = pexpect.spawn("/Users/vs/daophot/allstar")
    allstar.logfile = sys.stdout
    allstar.expect("OPT>")
    allstar.sendline("")
    allstar.expect("Input image name")
    allstar.sendline(image)
    allstar.expect("File with the PSF")
    allstar.sendline(image + ".psf")
    allstar.expect("Input file")
    allstar.sendline(image +"s.ap")
    allstar.expect("File for results")
    allstar.sendline(image +".als")
    allstar.expect("This file already exists")
    allstar.sendline("")
    allstar.expect("Name for subtracted image")
    allstar.sendline(image +"s.fits")

    allstar.expect("stars")
    allstar.expect("Good bye")
    allstar.close()
    print allstar.isalive()
    print "ALLSTAR complete"

## Calculate the offset for the mag file coords

    daomatch = pexpect.spawn("/Users/vs/daophot/daomatch")
    daomatch.logfile = sys.stdout

    if(is1!=None):
        epoch1 = imstem + "_e1_3p6um_dn"
    else:
        epoch1 = imstem + "_e1_4p5um_dn"

    if(os.path.isfile("offset.mch")):
        os.remove("offset.mch")

    daomatch.expect("Master input file")
    daomatch.sendline(epoch1 + ".als")
    daomatch.expect("Output file name")
    daomatch.sendline("offset.mch")
    daomatch.expect("Next input file")
    daomatch.sendline(image +'.als')
    daomatch.expect(" Next input file")
    daomatch.sendline("")

    daomatch.close(force=True)

    daomaster = pexpect.spawn("/Users/vs/daophot/daomaster")
    daomaster.expect("File with list of input files")
    daomaster.sendline("offset.mch")
    daomaster.expect("Minimum number, minimum fraction, enough frames")
    daomaster.sendline("2,1,2")
    daomaster.expect("Maximum sigma")
    daomaster.sendline("99")
## desired degrees of freedom:
    daomaster.expect("Your choice")
    daomaster.sendline("4")
    daomaster.expect("Critical match-up radius")
    daomaster.sendline("10")

    for radius in range (10,-1, -1):
        daomaster.expect("New match-up radius")
        daomaster.sendline(str(radius))

## Only need the transformations here
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
    daomaster.sendline("offset.mch_new")
    daomaster.expect("A file with the transfer table")
    daomaster.sendline("e")

    daomaster.close(force=True)

    new_mch  = glob.glob('offset.mch_new*')
    print new_mch
    shutil.move(str(new_mch[0]), "offset.mch")

    xoff, yoff = np.loadtxt("offset.mch", skiprows=1, usecols=(2, 3), unpack='TRUE')
    print xoff, yoff

## Finally apply the transformation to the med als file and create the mag file

    daophot = pexpect.spawn("/Users/vs/daophot/daophot")
    daophot.logfile = sys.stdout

    stem = re.sub("_e1", "", epoch1)
    stem = re.sub("_dn", "", stem)
    if os.path.isfile(stem + '.mag'): os.remove(stem+'.mag')
	
    daophot.expect("Command:")
    daophot.sendline("offset")
    daophot.expect("Input file name")
    daophot.sendline(image +".als")
    daophot.expect("Additive offsets ID, DX, DY, DMAG")
    daophot.sendline("0, " + str(xoff)+ ", " + str(yoff) + ", 0")
    daophot.expect("Output file name")
    daophot.sendline(stem  +".mag")

    daophot.expect("Command")
    daophot.sendline("exit")
    daophot.close(force=True) 


