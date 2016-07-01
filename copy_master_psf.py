#!/usr/bin/env python

import os
import sys
import re
import glob
import shutil

if (os.path.isfile("master_3p6um.psf")):
	files = glob.glob('*e*_3p6um_dn.fits')
	for file in files:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy("master_3p6um.psf", psfname)

if (os.path.isfile("master_4p5um.psf")):
	files = glob.glob('*e*_4p5um_dn.fits')
	for file in files:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy("master_4p5um.psf", psfname)

