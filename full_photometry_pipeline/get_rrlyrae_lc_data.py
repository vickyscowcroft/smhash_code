#!/usr/bin/env python

import aplpy
import astropy

## input = target_name from original script, image_stem
def get_lc_data(target_name, image_stem)

	target = target_name[0:11] + '.' + target_name[12:20]
	print target
	
	image_stem = flux_image[0:7]
	fitsfile = flux_image + '.fits'
	with open('/Users/vs/Dropbox/SMHASH/catalina_sagittarius_rrlyrae_catalogue', 'r') as searchfile:
		for line in searchfile:
			if target in line:
				data = line.split()
				ra = float(data[1])
				dec = float(data[2])
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
		distance = sqrt((rr_x - objects[star][1])**2 + (rr_y - objects[star][2])**2)
		if distance < min_distance:
			min_distance = distance
			best_match = objects[star][0]

	times = np.zeros(num_frames)
	
	for name in np.arange(0, num_frames):
		fits_name = re.sub(".cal", ".fits" file_list[name])
		hdulist = fits.open(fits_name)
		prihdr = hdulist[0].header
		n_bcds = prihdr['totalbcd']
		times[name] = prihdr['date_obs']
	
	times = astropy.time.Time(times, format='isot', scale='utc')
	mjds = times.mjd		
	print times
	print mjds
	
		
			
	
	
			 
	
	