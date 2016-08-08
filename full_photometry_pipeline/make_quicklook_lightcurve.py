#!/usr/bin/env python

import numpy as np
from astropy.time import Time
from astropy.io import fits
import glob
import itertools
import matplotlib.pyplot as mp


input_target = sys.argv[1]
input_channel = sys.argv[2]
catalog = sys.argv[3]

if (input_channel == '1' or input_channel == '3p6um'):
	channel = '3p6um'
if (input_channel =='2' or input_channel == '4p5um'):
	channel = '4p5um'

input_image_list = glob.glob(input_target + '_e*_' + channel + '_dn.fits')
input_lc_data = input_target + '_' + channel + '.cal'

mjds = np.zeros(len(input_image_list))

for count in range(0, len(input_image_list)):
	hdulist = fits.open(input_image_list[count])
	prihdr = hdulist[0].header
	date_obs = Time(prihdr['date_obs'])
	mjds[count] = date_obs.mjd
	
epoch1_ids, kal_ids = np.loadtxt(input_target + '_' + channel + '_rrl' + .tfr', skiprows=14, usecols=(0, 15), unpack=True)

k_ids, periods = np.loadtxt(catalog, usecols=(0,5), unpack=True)

with open(input_lc_data) as input:
	lines = []
	while True:
		line = list(itertools.islice(input, 3))
		if line:
			lines.append(line)
		else:
			break

num_stars = len(lines)

for rrl in range(len(epoch1_ids)):
	kaluzny = kal_ids[rrl]
	period = periods[np.where(kal_ids==kaluzny)]
	num_frames = 12
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
		if objects[star][0] == epoch1_ids[rrl]:
			best_star = star
			
			break
		
	mag = np.zeros(num_frames)
	err = np.zeros(num_frames)
	for epoch in np.arange(0,num_frames):
		mag[epoch] = objects[best_star][(epoch*2) + 3]  
		err[epoch] = objects[best_star][(epoch*2) + 4] 

	idx = np.argsort(mjds)
	mjds = mjds[idx]
	mag = mag[idx]
	err = err[idx]
	phase = (mjds / period) - np.floor(mjds / period)
	phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))
	
	mag_long  = np.concatenate((mag, mag, mag, mag, mag))
	
	err_long = np.concatenate((err, err, err, err, err))
	obs = np.arange(1, num_frames+1, 1)


	output = str(kaluzny) + '_rrlyrae.data'
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
	mp.title(input_target + ', P = ' + str(period) +' d')
	#mp.show()

	mp.savefig(str(kaluzny) +'.pdf')

	mp.close()
	mp.clf()
	axp2 = mp.subplot(111)
	axp2.invert_yaxis()

	dates = np.concatenate((mjds, mjds, mjds, mjds, mjds))
	axp2.axis([0,num_frames+1,(np.average(mag) + 0.4),(np.average(mag) - 0.4)])
	axp2.errorbar(obs, mag, yerr=err, ls='None')

	axp2.plot(obs, mag, 'ko', ls='None')
	mp.xlabel('Observation Number')
	if (channel == '3p6um'): mp.ylabel('[3.6]')
	if (channel == '4p5um'): mp.ylabel('[4.5]')
	mp.title(input_target + ', P = ' + str(period) +' d')
	mp.savefig(str(kaluzny) +'_timeseries.pdf')


	


	

