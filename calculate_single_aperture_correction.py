#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp

input = sys.argv[1]

apc, eapc, alf, ealf = np.loadtxt(input, skiprows=3, usecols=(3, 4, 5, 6), unpack='TRUE')

mp.close()

fig = mp.figure(figsize=(10,10))

difference  = apc - alf
av_diff = np.average(difference[(apc > 16.5) & (apc < 20.5)])
sdev_diff = np.std(difference[(apc > 16.5) & (apc < 20.5)])

#av_diff = np.average(difference)
#sdev_diff = np.std(difference)


total_err = np.sqrt(eapc**2 + ealf**2)

axp1 = mp.subplot(111)

axp1.errorbar(apc, difference, yerr = total_err, color='k', ls='none')
axp1.plot(apc, difference, 'k.', ls='none')
axp1.axhline(av_diff, color='r', ls='--')

mp.show()

print av_diff, sdev_diff