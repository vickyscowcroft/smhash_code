#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp
import glob
import os
import re

input = sys.argv[1]
channel = sys.argv[2]

file = open(input, "r")

## This needs to be updated to take a specific list or to determine the list from the filename
## Right now it's going to screw up if there's more than one field in the folder.

if channel == '1':
	mtrFiles = sorted(glob.glob('*_3p6um_dn.mtr'))
if channel == '2':
	mtrFiles = sorted(glob.glob('*_4p5um_dn.mtr'))


id = []
cal = []
ecal = []
alf2= []
ealf2 = []
alf3= []
ealf3 = []
alf4= []
ealf4 = []
alf5= []
ealf5 = []
alf6= []
ealf6 = []
alf7= []
ealf7 = []
alf8= []
ealf8 = []
alf9= []
ealf9 = []
alf10= []
ealf10 = []
alf11= []
ealf11 = []
alf12= []
ealf12 = []
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
	cal.append(float(data1[3]))
	ecal.append(float(data1[4]))
	alf2.append(float(data1[5]))
	ealf2.append(float(data1[6]))
	alf3.append(float(data1[7]))
	ealf3.append(float(data1[8]))
	alf4.append(float(data1[9]))
	ealf4.append(float(data1[10]))
	alf5.append(float(data1[11]))
	ealf5.append(float(data1[12]))
	alf6.append(float(data1[13]))
	ealf6.append(float(data1[14]))
	
	data2 = secondline.split()
	alf7.append(float(data2[0]))
	ealf7.append(float(data2[1]))
	alf8.append(float(data2[2]))
	ealf8.append(float(data2[3]))
	alf9.append(float(data2[4]))
	ealf9.append(float(data2[5]))
	alf10.append(float(data2[6]))
	ealf10.append(float(data2[7]))
	alf11.append(float(data2[8]))
	ealf11.append(float(data2[9]))
	alf12.append(float(data2[10]))
	ealf12.append(float(data2[11]))
	

	firstline = file.readline().strip()
	
mp.close()

id = np.array(id)
cal = np.array(cal)
ecal = np.array(ecal)
alf2= np.array(alf2)
ealf2 = np.array(ealf2)
alf3= np.array(alf3)
ealf3 = np.array(ealf3)
alf4= np.array(alf4)
ealf4 = np.array(ealf4)
alf5= np.array(alf5)
ealf5 = np.array(ealf5)
alf6= np.array(alf6)
ealf6 = np.array(ealf6)
alf7= np.array(alf7)
ealf7 = np.array(ealf7)
alf8= np.array(alf8)
ealf8 = np.array(ealf8)
alf9= np.array(alf9)
ealf9 = np.array(ealf9)
alf10= np.array(alf10)
ealf10 = np.array(ealf10)
alf11= np.array(alf11)
ealf11 = np.array(ealf11)
alf12= np.array(alf12)
ealf12 = np.array(ealf12)
xc = np.array(xc)
yc = np.array(yc)


fig = mp.figure(figsize=(10,10))
"""
difference  = apc - alf
av_diff = np.average(difference[(apc > 16.5) & (apc < 20.5)])
sdev_diff = np.std(difference[(apc > 16.5) & (apc < 20.5)])

#av_diff = np.np.average(difference)
#sdev_diff = np.std(difference)


total_err = sqrt(eapc**2 + ealf**2)
"""
axp1 = mp.subplot(6,2,1)
axp2 = mp.subplot(6,2,2)
axp3 = mp.subplot(6,2,3)
axp4 = mp.subplot(6,2,4)
axp5 = mp.subplot(6,2,5)
axp6 = mp.subplot(6,2,6)
axp7 = mp.subplot(6,2,7)
axp8 = mp.subplot(6,2,8)
axp9 = mp.subplot(6,2,9)
axp10 = mp.subplot(6,2,10)
axp11 = mp.subplot(6,2,11)


difference2 = cal - alf2
ediff2 = np.sqrt(ecal**2 + ealf2)
av_diff2 = np.average(difference2[(cal >10) & (cal < 16) & (alf2 < 99)])
sdev2 = np.std(difference2[(cal >10) & (cal < 16) & (alf2 < 99)])

difference3 = cal - alf3
ediff3 = np.sqrt(ecal**2 + ealf3)
av_diff3 = np.average(difference3[(cal >10) & (cal < 16) & (alf3 < 99)])
sdev3 = np.std(difference3[(cal >10) & (cal < 16) & (alf3 < 99)])

difference4 = cal - alf4
ediff4 = np.sqrt(ecal**2 + ealf4)
av_diff4 = np.average(difference4[(cal >10) & (cal < 16) & (alf4 < 99)])
sdev4 = np.std(difference4[(cal >10) & (cal < 16) & (alf4 < 99)])

difference5 = cal - alf5
ediff5 = np.sqrt(ecal**2 + ealf5)
av_diff5 = np.average(difference5[(cal >10) & (cal < 16) & (alf5 < 99)])
sdev5 = np.std(difference5[(cal >10) & (cal < 16) & (alf5 < 99)])

difference6 = cal - alf6
ediff6 = np.sqrt(ecal**2 + ealf6)
av_diff6 = np.average(difference6[(cal >10) & (cal < 16) & (alf6 < 99)])
sdev6 = np.std(difference6[(cal >10) & (cal < 16) & (alf6 < 99)])

difference7 = cal - alf7
ediff7 = np.sqrt(ecal**2 + ealf7)
av_diff7 = np.average(difference7[(cal >10) & (cal < 16) & (alf7 < 99)])
sdev7 = np.std(difference7[(cal >10) & (cal < 16) & (alf7 < 99)])

difference8 = cal - alf8
ediff8 = np.sqrt(ecal**2 + ealf8)
av_diff8 = np.average(difference8[(cal >10) & (cal < 16) & (alf8 < 99)])
sdev8 = np.std(difference8[(cal >10) & (cal < 16) & (alf8 < 99)])

difference9 = cal - alf9
ediff9 = np.sqrt(ecal**2 + ealf9)
av_diff9 = np.average(difference9[(cal >10) & (cal < 16) & (alf9 < 99)])
sdev9 = np.std(difference9[(cal >10) & (cal < 16) & (alf9 < 99)])

difference10 = cal - alf10
ediff10 = np.sqrt(ecal**2 + ealf10)
av_diff10 = np.average(difference10[(cal >10) & (cal < 16) & (alf10 < 99)])
sdev10 = np.std(difference10[(cal >10) & (cal < 16) & (alf10 < 99)])

difference11 = cal - alf11
ediff11 = np.sqrt(ecal**2 + ealf11)
av_diff11 = np.average(difference11[(cal >10) & (cal < 16) & (alf11 < 99)])
sdev11 = np.std(difference11[(cal >10) & (cal < 16) & (alf11 < 99)])

difference12 = cal - alf12
ediff12 = np.sqrt(ecal**2 + ealf12)
av_diff12 = np.average(difference12[(cal >10) & (cal < 16) & (alf12 < 99)])
sdev12 = np.std(difference12[(cal >10) & (cal < 16) & (alf12 < 99)])

axp1.errorbar(cal, difference2, yerr = ediff2, color='k', ls='none')
axp1.plot(cal, difference2, 'k.', ls='none')
axp1.axhline(av_diff2, color='r', ls='-')
axp1.axhline(av_diff2 + 2*sdev2, color='r', ls='--')
axp1.axhline(av_diff2 - 2*sdev2, color='r', ls='--')
axp1.annotate("2", xy=(10,4), xycoords='data', ha='center')


axp2.errorbar(cal, difference3, yerr = ediff3, color='k', ls='none')
axp2.plot(cal, difference3, 'k.', ls='none')
axp2.axhline(av_diff3, color='r', ls='-')
axp2.axhline(av_diff3 + 2*sdev3, color='r', ls='--')
axp2.axhline(av_diff3 - 2*sdev3, color='r', ls='--')
axp2.annotate("3", xy=(10,4), xycoords='data', ha='center')

axp3.errorbar(cal, difference4, yerr = ediff4, color='k', ls='none')
axp3.plot(cal, difference4, 'k.', ls='none')
axp3.axhline(av_diff4, color='r', ls='-')
axp3.axhline(av_diff4 + 2*sdev4, color='r', ls='--')
axp3.axhline(av_diff4 - 2*sdev4, color='r', ls='--')
axp3.annotate("4", xy=(10,4), xycoords='data', ha='center')

axp4.errorbar(cal, difference5, yerr = ediff5, color='k', ls='none')
axp4.plot(cal, difference5, 'k.', ls='none')
axp4.axhline(av_diff5, color='r', ls='-')
axp4.axhline(av_diff5 + 2*sdev5, color='r', ls='--')
axp4.axhline(av_diff5 - 2*sdev5, color='r', ls='--')
axp4.annotate("5", xy=(10,4), xycoords='data', ha='center')

axp5.errorbar(cal, difference6, yerr = ediff6, color='k', ls='none')
axp5.plot(cal, difference6, 'k.', ls='none')
axp5.axhline(av_diff6, color='r', ls='-')
axp5.axhline(av_diff6 + 2*sdev6, color='r', ls='--')
axp5.axhline(av_diff6 - 2*sdev6, color='r', ls='--')
axp5.annotate("6", xy=(10,4), xycoords='data', ha='center')

axp6.errorbar(cal, difference7, yerr = ediff7, color='k', ls='none')
axp6.plot(cal, difference7, 'k.', ls='none')
axp6.axhline(av_diff7, color='r', ls='-')
axp6.axhline(av_diff7 + 2*sdev7, color='r', ls='--')
axp6.axhline(av_diff7 - 2*sdev7, color='r', ls='--')
axp6.annotate("7", xy=(10,4), xycoords='data', ha='center')

axp7.errorbar(cal, difference8, yerr = ediff8, color='k', ls='none')
axp7.plot(cal, difference8, 'k.', ls='none')
axp7.axhline(av_diff8, color='r', ls='-')
axp7.axhline(av_diff8 + 2*sdev8, color='r', ls='--')
axp7.axhline(av_diff8 - 2*sdev8, color='r', ls='--')
axp7.annotate("8", xy=(10,4), xycoords='data', ha='center')

axp8.errorbar(cal, difference9, yerr = ediff9, color='k', ls='none')
axp8.plot(cal, difference9, 'k.', ls='none')
axp8.axhline(av_diff9, color='r', ls='-')
axp8.axhline(av_diff9 + 2*sdev9, color='r', ls='--')
axp8.axhline(av_diff9 - 2*sdev9, color='r', ls='--')
axp8.annotate("9", xy=(10,4), xycoords='data', ha='center')

axp9.errorbar(cal, difference10, yerr = ediff10, color='k', ls='none')
axp9.plot(cal, difference10, 'k.', ls='none')
axp9.axhline(av_diff10, color='r', ls='-')
axp9.axhline(av_diff10 + 2*sdev10, color='r', ls='--')
axp9.axhline(av_diff10 - 2*sdev10, color='r', ls='--')
axp9.annotate("10", xy=(10,4), xycoords='data', ha='center')

axp10.errorbar(cal, difference11, yerr = ediff11, color='k', ls='none')
axp10.plot(cal, difference11, 'k.', ls='none')
axp10.axhline(av_diff11, color='r', ls='-')
axp10.axhline(av_diff11 + 2*sdev11, color='r', ls='--')
axp10.axhline(av_diff11 - 2*sdev11, color='r', ls='--')
axp10.annotate("11", xy=(10,4), xycoords='data', ha='center')

axp11.errorbar(cal, difference12, yerr = ediff12, color='k', ls='none')
axp11.plot(cal, difference12, 'k.', ls='none')
axp11.axhline(av_diff12, color='r', ls='-')
axp11.axhline(av_diff12 + 2*sdev12, color='r', ls='--')
axp11.axhline(av_diff12 - 2*sdev12, color='r', ls='--')
axp11.annotate("12", xy=(10,4), xycoords='data', ha='center')

mp.show()

print "epoch 2: ", av_diff2, sdev2
print "epoch 3: ", av_diff3, sdev3
print "epoch 4: ", av_diff4, sdev4
print "epoch 5: ", av_diff5, sdev5
print "epoch 6: ", av_diff6, sdev6
print "epoch 7: ", av_diff7, sdev7
print "epoch 8: ", av_diff8, sdev8
print "epoch 9: ", av_diff9, sdev9
print "epoch 10: ", av_diff10, sdev10
print "epoch 11: ", av_diff11, sdev11
print "epoch 12: ", av_diff12, sdev12

cal2 = alf2 + av_diff2
cal3 = alf3 + av_diff3
cal4 = alf4 + av_diff4
cal5 = alf5 + av_diff5
cal6 = alf6 + av_diff6
cal7 = alf7 + av_diff7
cal8 = alf8 + av_diff8
cal9 = alf9 + av_diff9
cal10 = alf10 + av_diff10
cal11 = alf11 + av_diff11
cal12 = alf12 + av_diff12

av_diff1 = 0

av_diffs = [av_diff10, av_diff11, av_diff12, av_diff1, av_diff2, av_diff3, av_diff4, av_diff5, av_diff6, av_diff7, av_diff8, av_diff9]

count = 0
for name in mtrFiles:
	id, xc, yc, mag, err = np.loadtxt(name, skiprows=3, usecols=(0, 1, 2, 3, 4), unpack='TRUE')
	newmag = mag + av_diffs[count]
	output_name = re.sub('.mtr', '.off',name)
	np.savetxt(output_name, np.column_stack((id, xc, yc, newmag, err)), fmt= "%d %.2f %.2f %.3f %.3f")
	count = count + 1

print "Finished"


#np.savetxt(output, np.column_stack((id, xc, yc, cal, ecal, cal2, ealf2, cal3, ealf3, cal4, ealf4, cal5, ealf5, cal6, ealf6, cal7, ealf7, cal8, ealf8, cal9, ealf9, cal10, ealf10, cal11, ealf11, cal12, ealf12, )), fmt= "%d %.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f ")
