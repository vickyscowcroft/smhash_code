import numpy as np
import pandas as pd
import matplotlib.pyplot as mp
import aplpy
import montage_wrapper as montage
from astropy.coordinates import SkyCoord
import re
import os
import sgr_setup
import calculate_aperture_correction_testing
import all_location_corrections
import calibrate_all_epochs
import glob
import aperture_correction_photometry
import apply_aperture_correction
import itertools
import astropy
from astropy.time import Time
from astropy.visualization import ZScaleInterval
import astropy.units as u
import gloess_fits as gf




def make_rrl_finder_chart(target_list, counter):
	ra = target_list.ix[counter, 'ra']
	dec = target_list.ix[counter, 'dec']
	target_stem = target_list.ix[counter, 'survey_ident']
	target_stem = re.sub('\.', '_', target_stem)
	target_stem = re.sub('-', '', target_stem)
	target_stem = sgr_setup.get_target_stem(target_stem)
	coords = str(ra) +' ' +str(dec)
	ra= SkyCoord(coords, unit=(u.deg, u.deg)).ra
	dec = SkyCoord(coords, unit=(u.deg, u.deg)).dec
	fitsfile = target_stem + '_e1_3p6um.fits'
	print ra, dec
	inputfile = fitsfile
	fitsdata = astropy.io.fits.open(fitsfile)[0].data
	interval = ZScaleInterval()
	zmin, zmax = interval.get_limits(fitsdata)
	fig = mp.figure(figsize=(10,10))
	mosaic = aplpy.FITSFigure(inputfile, figure = fig)
	mosaic.show_grayscale(vmin=zmin,vmax=zmax, invert='true') ### manually implimenting zscale
	mosaic.tick_labels.set_font(size='small')
	mosaic.tick_labels.set_xformat("hh:mm:ss")
	mosaic.set_theme('publication')
	mosaic.show_markers(ra.deg, dec.deg, edgecolor='magenta', facecolor='magenta', marker='o', s=100, alpha=0.3)
	mosaic.show_markers(ra.deg, dec.deg, edgecolor='magenta', facecolor='magenta', marker='o', s=300, alpha=0.1)
	#mosaic.save(target_stem + '_location.pdf')
	mp.show()
    
def directory_switcher(target_list, counter):
    survey_ident = target_list.ix[counter, 'survey_ident']
    dir_name = re.sub('\.', '_', survey_ident)
    dir_name = re.sub('-', '', dir_name)
    new_chan = '3p6um' ## hardwired to channel 1
    directory = '/Users/vs522/Dropbox/TRACCS/TRACCS_Output/' + dir_name + '/' + dir_name + '_' + new_chan
    os.chdir(directory)

    return(os.getcwd()) ### Print out the current directory at the end so I know I'm in the right place


def find_which_star(target_stem, ra, dec):
	fitsfile = target_stem + '_e1_3p6um.fits'
	fig = mp.figure(figsize=(10,10))
	mosaic = aplpy.FITSFigure(fitsfile)
	rr_x, rr_y = mosaic.world2pixel(ra, dec)

	min_distance = 10000.
	num_frames = len(glob.glob('*_dn.cal'))
	if ((num_frames % 6 ) == 0):
		n_lines_raw = int((num_frames / 6) + 1)
	else:
		n_lines_raw = int(np.ceil(num_frames) / 6 + 1)
	print n_lines_raw
	phot_file = target_stem + '_3p6um.cal'


	with open(phot_file) as input:
		lines = []
		while True:
			line = list(itertools.islice(input, n_lines_raw))
			if line:
				lines.append(line)
			else:
				break
	num_stars = len(lines)



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
        
	return (mag, err)


def grab_mjds(mch_file):
    mch_file = open(mch_file, 'r')
    num_frames = len(glob.glob('*_dn.cal'))

    file_list = []
    for count in np.arange(0, num_frames):
        line = mch_file.readline()
        lines = line.split("\'")[1]
        file_list.append(lines)
    times = []
    num_frames = len(glob.glob('*_dn.cal'))
    for name in np.arange(0, num_frames):
        fits_name = re.sub(".cal", ".fits", file_list[name])
        fits_name = re.sub(" ", "", fits_name)

        hdulist = astropy.io.fits.open(fits_name)
        prihdr = hdulist[0].header
        times.append(prihdr['date_obs'])

    times = Time(times, format='isot', scale='utc')
    mjds = times.mjd
    return(mjds)

def plot_a_lc(target_list, counter, mjds, mag, err):
	ra = targe_list.ix[counter, 'ra']
	dec = target_list.ix[counter, 'dec']
	target_stem = target_list.ix[counter, 'survey_ident']
	target_stem = re.sub('\.', '_', target_stem)
	target_stem = re.sub('-', '', target_stem)
	target_name = target_list.ix[counter, 'target_name']
	period = target_list.ix[counter, 'period']
	idx = np.argsort(mjds)
	mjds = mjds[idx]
	mag = mag[idx]
	err = err[idx]

	phase = (mjds / period) - np.floor(mjds / period)
	phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

	mag_long  = np.concatenate((mag, mag, mag, mag, mag))

	err_long = np.concatenate((err, err, err, err, err))
	obs = np.arange(1, num_frames+1, 1)

	mp.close()
	mp.clf()
	mp.axis([0,2.5,(np.average(mag) + 0.4),(np.average(mag) - 0.4)])
	mp.errorbar(phase, mag_long, yerr=err_long, ls='None')
	mp.plot(phase, mag_long, 'ko', ls='None')
	mp.xlabel("Phase")
	mp.ylabel('[3.6]')
	mp.title(target_name + ', P = ' + str(period) +' d')
	outplot = target_name + '_lc.pdf'
	mp.savefig(outplot)
    
	output = target_name + '_rrlyrae.data'
	output_file = open(output, "w")

	for frame in np.arange(0, num_frames):
		output_file.write("{0:.8f} {1:.3f} {2:.3f} \n".format(mjds[frame], mag[frame], err[frame]))

	output_file.close()

def calibrate_and_plot(row):

	directory_switcher(row, 0)

	survey_ident = row.ix[0,'survey_ident']
	ra = row.ix[0,'ra']
	dec = row.ix[0,'dec']
	period = row['period']
	target_name = row['target_name']
    
	channel = 1 ### hardwired to channel 1 here
	if (channel == 1 or channel == '3p6um' or channel == '1'): new_chan = '3p6um'
	elif (channel == 2 or channel == '4p5um' or channel == '2'): new_chan = '4p5um'
	else: 
		print 'invalid channel'
		exit(1)

	if new_chan == '3p6um': num_chan = 1
	if new_chan == '4p5um': num_chan = 2
	target_stem = sgr_setup.get_target_stem(survey_ident)
	dir_name = re.sub('\.', '_', survey_ident)
	dir_name = re.sub('-', '', dir_name)
	directory = '/Users/vs522/Dropbox/TRACCS/TRACCS_Output/' + dir_name + '/' + dir_name + '_' + new_chan
	os.chdir(directory)
	if (len(glob.glob(target_stem + '*' + new_chan +'*.alf'))==0):
		print 'You need to run ALLFRAME before you can calibrate the photometry'
		print 'This is the CALIBRATION ONLY script'
		print 'Run ``globular_pipeline_anychannel.py`` to do the combined photometry and calibration script'
		exit(1)
	flux_epoch_1 = target_stem + '_e1_' + new_chan
	aperture_correction_photometry.apcor_photo(flux_epoch_1)
	apcor, sdev_apcor = calculate_aperture_correction_testing.calc_apcor(flux_epoch_1, 'apcor.raw', 3, target_stem, 0.3)
	apply_aperture_correction.apply_apcor(target_stem + '_e1_' + new_chan + '_dn.alf', apcor, num_chan)
	calibrate_all_epochs.calibrate(target_stem + '_' + new_chan + '.mch')
	all_location_corrections.location_corr(target_stem, new_chan)
    
    
    #make_rrl_finder_chart(target_stem, ra, dec)

	mch_file = target_stem + '_3p6um_cal.mch'
	phot_file = target_stem + '_3p6um.cal'


	mag, err = find_which_star(target_stem, ra, dec)
	mjds = grab_mjds(mch_file)
    
    #plot_a_lc(target_stem, target_name, period, mjds, mag, err)
    
	
	#mosaic.save(target_stem + '_location.pdf')
	make_rrl_finder_chart(row, 0)

	offsets_df = pd.read_csv(target_stem + '_3p6um_cal.off_log', delim_whitespace=True, header=None, names=('obs', 'offset', 'sdev', 'nstars'))
	
	mp.close()
	mp.clf()
	mp.errorbar(offsets_df.obs, offsets_df.offset, yerr=offsets_df.sdev, color='k', ls='none')
	mp.plot(offsets_df.obs, offsets_df.offset, 'ko', ls='none')
	mp.axhline(offsets_df.offset.mean(), color='r', ls='--')
	mp.axhline(offsets_df.offset.mean()- 2*offsets_df.offset.std(), color='b', ls='--')

	mp.show()
	
	plot_existing_lc(row, 1, 1.5/12.)
	plot_new_lc(row, mjds, mag, err, 1, 1.5/12.)	
	return(0)
	
def plot_existing_lc(row, phased=1, smooth=0.1):

    target_name = row.ix[0, 'target_name']
    rrl_data = target_name + '_rrlyrae.data'
    
    rrl_df = pd.read_csv(rrl_data, delim_whitespace=True, header=None, names=('mjds', 'mag', 'err'))
    
    period = row.ix[0, 'period']
    idx = np.argsort(rrl_df.mjds)
    mjds = rrl_df.mjds[idx]
    mag = rrl_df.mag[idx]
    err = rrl_df.err[idx]

    phase = (mjds / period) - np.floor(mjds / period)
    phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

    mag_long  = np.concatenate((mag, mag, mag, mag, mag))
    
    nir1 = len(mag)

    err_long = np.concatenate((err, err, err, err, err))
    #obs = np.arange(1, num_frames+1, 1)
    ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(mag_long,err_long, phase,len(mag),smooth)
    aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)
    if phased == 1:
        factor = np.sqrt(nir1)
    if phased == 0:
        factor = 1
    sdevir1 = sdevir1/factor

    target_name_display = re.sub('_', ' ', target_name)

    mp.close()
    mp.clf()
    mp.axis([1,3.5,(np.average(mag) + 0.3),(np.average(mag) - 0.3)])
    mp.plot(ir1x,ir11,'k-')
    mp.errorbar(phase, mag_long, yerr=err_long, ls='None', zorder=4)
    mp.plot(phase, mag_long, 'ro',  ls='None', zorder=4)
    mp.axhline(aveir1, color='k',ls='--')
    mp.xlabel("Phase")
    mp.ylabel('[3.6]')
    mp.title('Existing:' + target_name_display + ', P = ' + str(np.around(period, decimals=4)) +' d, [3.6] = ' + str(np.around(aveir1, decimals=3)) + ' $\pm$ ' + str(np.around(sdevir1, decimals=3)) + ' mag')
    mp.show()
    
    return(0)
    
def plot_new_lc(row, mjds, mag, err, phased=1, smooth=0.1):

    target_name = row.ix[0, 'target_name']
    period = row.ix[0, 'period']
    idx = np.argsort(mjds)
    mjds = mjds[idx]
    mag = mag[idx]
    err = err[idx]

    phase = (mjds / period) - np.floor(mjds / period)
    phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))

    mag_long  = np.concatenate((mag, mag, mag, mag, mag))
    
    nir1 = len(mag)

    err_long = np.concatenate((err, err, err, err, err))
    #obs = np.arange(1, num_frames+1, 1)
    ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(mag_long,err_long, phase,len(mag),smooth)
    aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)
    if phased == 1:
        factor = np.sqrt(nir1)
    if phased == 0:
        factor = 1
    sdevir1 = sdevir1/factor

    target_name_display = re.sub('_', ' ', target_name)

    mp.close()
    mp.clf()
    mp.axis([1,3.5,(np.average(mag) + 0.3),(np.average(mag) - 0.3)])
    mp.plot(ir1x,ir11,'k-')
    mp.errorbar(phase, mag_long, yerr=err_long, ls='None', zorder=4)
    mp.plot(phase, mag_long, 'ro',  ls='None', zorder=4)
    mp.axhline(aveir1, color='k',ls='--')
    mp.xlabel("Phase")
    mp.ylabel('[3.6]')
    mp.title('Recal:' + target_name_display + ', P = ' + str(np.around(period, decimals=4)) +' d, [3.6] = ' + str(np.around(aveir1, decimals=3)) + ' $\pm$ ' + str(np.around(sdevir1, decimals=3)) + ' mag')
    mp.show()
    
    return(0)



	

	

 
    
    

    