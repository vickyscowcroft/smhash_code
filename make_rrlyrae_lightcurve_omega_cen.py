#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp
import glob
import os
import re

input_raw = sys.argv[1]
field = sys.argv[2]

file = open(input_raw, "r")


id = []
rrv = []
rrverr = []
cal1 = []
ecal1 = []
cal2= []
ecal2 = []
cal3= []
ecal3 = []
cal4= []
ecal4 = []
cal5= []
ecal5 = []
cal6= []
ecal6 = []
cal7= []
ecal7 = []
cal8= []
ecal8 = []
cal9= []
ecal9 = []
cal10= []
ecal10 = []
cal11= []
ecal11 = []
cal12= []
ecal12 = []

xc = []
yc = []


firstline = file.readline().strip()

while firstline != "":
	secondline = file.readline().strip()
	thirdline = file.readline().strip()
	
	data1 = firstline.split()
	id.append(int(data1[0]))
	xc.append(float(data1[1]))
	yc.append(float(data1[2]))
	rrv.append(float(data1[3]))
	rrverr.append(float(data1[4]))
	cal1.append(float(data1[5]))
	ecal1.append(float(data1[6]))
	cal2.append(float(data1[7]))
	ecal2.append(float(data1[8]))
	cal3.append(float(data1[9]))
	ecal3.append(float(data1[10]))
	cal4.append(float(data1[11]))
	ecal4.append(float(data1[12]))
	cal5.append(float(data1[13]))
	ecal5.append(float(data1[14]))
	
	data2 = secondline.split()
	cal6.append(float(data2[0]))
	ecal6.append(float(data2[1]))
	cal7.append(float(data2[2]))
	ecal7.append(float(data2[3]))
	cal8.append(float(data2[4]))
	ecal8.append(float(data2[5]))
	cal9.append(float(data2[6]))
	ecal9.append(float(data2[7]))
	cal10.append(float(data2[8]))
	ecal10.append(float(data2[9]))
	cal11.append(float(data2[10]))
	ecal11.append(float(data2[11]))

	data3 = thirdline.split()
	cal12.append(float(data2[0]))
	ecal12.append(float(data2[1]))
		
	firstline = file.readline().strip()
	
mp.close()

id = np.array(id)
cal1 = np.array(cal1)
ecal1 = np.array(ecal1)
cal2= np.array(cal2)
ecal2 = np.array(ecal2)
cal3= np.array(cal3)
ecal3 = np.array(ecal3)
cal4= np.array(cal4)
ecal4 = np.array(ecal4)
cal5= np.array(cal5)
ecal5 = np.array(ecal5)
cal6= np.array(cal6)
ecal6 = np.array(ecal6)
cal7= np.array(cal7)
ecal7 = np.array(ecal7)
cal8= np.array(cal8)
ecal8 = np.array(ecal8)
cal9= np.array(cal9)
ecal9 = np.array(ecal9)
cal10= np.array(cal10)
ecal10 = np.array(ecal10)
cal11= np.array(cal11)
ecal11 = np.array(ecal11)
cal12= np.array(cal12)
ecal12 = np.array(ecal12)
xc = np.array(xc)
yc = np.array(yc)

kaluzny_file = "kaluzny_rrl_field" + field
print kaluzny_file

#k_id, ra, dec, period, type = np.loadtxt(kaluzny_file, usecols=(0, 1, 2, 3, 6 ), unpack='TRUE')

k_id, ra, dec, period, type = np.loadtxt(kaluzny_file, dtype = {'names': ('k_id', 'ra', 'dec', 'period', 'rrl_type'), 'formats': ('S5', 'S11', 'S12', 'f6', 'S3')}, usecols=(0, 1, 2, 3, 6 ), unpack='TRUE')

average_mags = []
average_errs = []
log_p = []

for star in range(size(id)):
	rrl_id = "V" + str(id[star])
	match = np.where(k_id == rrl_id)
	
	flux1 = 10**(cal1[star]/2.5)
	flux2 = 10**(cal2[star]/2.5)
	flux3 = 10**(cal3[star]/2.5)
	flux4 = 10**(cal4[star]/2.5)
	flux5 = 10**(cal5[star]/2.5)
	flux6 = 10**(cal6[star]/2.5)
	flux7 = 10**(cal7[star]/2.5)
	flux8 = 10**(cal8[star]/2.5)
	flux9 = 10**(cal9[star]/2.5)
	flux10 = 10**(cal10[star]/2.5)
	flux11 = 10**(cal11[star]/2.5)
	flux12 = 10**(cal12[star]/2.5)
	
	average_flux = (flux1 + flux2 + flux3 + flux4 + flux5 + flux6 + flux7 + flux8 + flux9 + flux10 + flux11 + flux12) / 12.
	
	average_err = (ecal1[star] + ecal2[star] + ecal3[star] + ecal4[star] + ecal5[star] + ecal6[star] + ecal7[star] + ecal8[star] + ecal9[star] + ecal10[star] + ecal11[star] + ecal12[star]) / sqrt(12.)
	
	average_mag = 2.5*np.log10(average_flux)
	
	average_mags.append(average_mag)
	average_errs.append(average_err)
	log_p.append(np.log10(period[match]))
	
	print k_id[match], period[match], average_mag, average_err, rrl_type[match]
	


	
	
	
	