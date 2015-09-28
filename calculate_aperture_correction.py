#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as mp
from astropy.stats import sigma_clip

## Updated 28/09/2015 to use sigmaclipping to remove bad points rather than relying on a hard wired magnitude range

input = sys.argv[1]

apc, eapc, alf, ealf = np.loadtxt(input, skiprows=3, usecols=(3, 4, 5, 6), unpack='TRUE')

mp.close()

fig = mp.figure(figsize=(10,10))

difference  = apc - alf
#av_diff = np.average(difference)
#sdev_diff = np.std(difference)

## Using astropy.stats.sigma_clip
## This version returns a masked array rather than a cut array
clipped = sigma_clip(difference, sig = 2.5, iters=100)

av_diff = np.average(clipped)
sdev_diff = np.std(clipped)


total_err = np.sqrt(eapc**2 + ealf**2)

axp1 = mp.subplot(111)

#axp1.errorbar(apc, difference, yerr = total_err, color='grey', ls='none')
#axp1.plot(apc, difference, 'k.', ls='none')
#axp1.axhline(av_diff, color='r', ls='--')

axp1.plot(apc, clipped, 'k.', ls='none')
axp1.axhline(av_diff, color='r', ls='--')
axp1.axhline(av_diff+2*sdev_diff, color='b', ls='--')
axp1.axhline(av_diff-2*sdev_diff, color='b', ls='--')


mp.show()

print av_diff, sdev_diff