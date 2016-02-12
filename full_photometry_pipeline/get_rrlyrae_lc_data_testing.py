#!/usr/bin/env python

import aplpy
import astropy
import itertools
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as mp
import sys
import glob
import re

## input = target_name from original script, image_stem


target_name = sys.argv[1]
location = sys.argv[2]


if(location=='orphan'):
	cat_target = target_name
	cat_name = '/Users/vs/Dropbox/SMHASH/orphan_rrlyrae_catalogue'
	image_stem = cat_target

if(location=='sgr_stream'):
	cat_target = target_name[0:11] + '.' + target_name[12:20]
	cat_name = '/Users/vs/Dropbox/SMHASH/catalina_sagittarius_rrlyrae_catalogue'
	image_stem = re.sub("CSS_","", cat_target)
	image_stem = new_target_stem[0:7]
	



target = cat_target
fitsfile = image_stem + '_e1_3p6um.fits'
with open(cat_name, 'r') as searchfile:
	for line in searchfile:
		if target in line:
			data = line.split()
			ra = float(data[1])
			dec = float(data[2])
			if(location=='orphan'):
				period = float(data[9])
			if(location=='sgr_stream'):
				period = float(data[4])
print ra, dec

fig = aplpy.FITSFigure(fitsfile)
rr_x, rr_y = fig.world2pixel(ra, dec)
print rr_x, rr_y
	
mch_file = image_stem + '_3p6um.mch'
phot_file = image_stem + '_3p6um.cal'
	
num_frames = len(glob.glob('*_dn.cal'))
	
file_list = []
mch_file = open(mch_file, 'r')
for count in np.arange(0, num_frames):
	line = mch_file.readline()
	lines = line.split("\'")[1]
	file_list.append(lines)
	
if ((num_frames % 6 ) == 0):
	n_lines_raw = int((num_frames / 6) + 1)
else:
	n_lines_raw = int(np.ceil(num_frames) / 6 + 1)
print n_lines_raw

	
with open(phot_file) as input:
	lines = []
	while True:
		line = list(itertools.islice(input, n_lines_raw))
		if line:
			lines.append(line)
		else:
			break
		#print lines
		
num_stars = len(lines)

## Read all relavent data into a numpy array
min_distance = 10000.

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
	distance = np.sqrt((rr_x - objects[star][1])**2 + (rr_y - objects[star][2])**2)
	if distance < min_distance:
		min_distance = distance
		best_match = objects[star][0]
		best_star = star
		
mag = np.zeros(num_frames)
err = np.zeros(num_frames)
for epoch in np.arange(0,num_frames):
	mag[epoch] = objects[best_star][(epoch*2) + 3]  
	err[epoch] = objects[best_star][(epoch*2) + 4] 

times = []
	
for name in np.arange(0, num_frames):
	fits_name = re.sub(".cal", ".fits", file_list[name])
	fits_name = re.sub(" ", "", fits_name)

	hdulist = astropy.io.fits.open(fits_name)
	prihdr = hdulist[0].header
	times.append(prihdr['date_obs'])
	
times = Time(times, format='isot', scale='utc')
mjds = times.mjd		
print times
print mjds

idx = np.argsort(mjds)
mjds = mjds[idx]
mag = mag[idx]
err = err[idx]

phase = (mjds / period) - np.floor(mjds / period)
phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

mag_long  = np.concatenate((mag, mag, mag, mag, mag))

err_long = np.concatenate((err, err, err, err, err))
obs = np.arange(1, num_frames+1, 1)


output = image_stem + '_rrlyrae.data'
output_file = open(output, "w")

for frame in np.arange(0, num_frames):
	output_file.write("{0:.8f} {1:.3f} {2:.3f} \n".format(mjds[frame], mag[frame], err[frame]))
	
output_file.close()
## Sanity check plot

mp.close()
mp.clf()
axp1 = mp.subplot(111)

axp1.axis([0,2.5,(np.average(mag) + 0.4),(np.average(mag) - 0.4)])


axp1.errorbar(phase, mag_long, yerr=err_long, ls='None')
axp1.plot(phase, mag_long, 'ko', ls='None')
mp.xlabel("Phase")
mp.ylabel('[3.6]')
mp.title(target + ', P = ' + str(period) +' d')
#mp.show()

mp.savefig(image_stem +'.pdf')

mp.close()
mp.clf()
axp2 = mp.subplot(111)
axp2.invert_yaxis()

dates = np.concatenate((mjds, mjds, mjds, mjds, mjds))
axp2.axis([0,num_frames+1,(np.average(mag) + 0.4),(np.average(mag) - 0.4)])
axp2.errorbar(obs, mag, yerr=err, ls='None')

axp2.plot(obs, mag, 'ko', ls='None')
mp.xlabel('Observation Number')
mp.ylabel('[3.6]')
mp.title(target + ', P = ' + str(period) +' d')
mp.savefig(image_stem +'_timeseries.pdf')

#mp.close()

			
	
	
			 
	
	