#!/usr/bin/env/python

## Want to read in a set of files to grab the period, mag and error
## for each cepheid. It's possible that these three variables will be different for each
## band in a given star

import numpy as np
import matplotlib.pyplot as mp
import glob
import re
import os
from scipy.optimize import curve_fit
import reddening_laws as red


mp.close('all')
os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

## Read in the data

## Will have 90+ files ( in generic cases this will be N files)
## One file for each Cepheid
## Number of lines changes depending on how many bands there are available
## Format of each line is the same unless there was only a single observation

files = glob.glob('*.glo_avs')


band = []
period = []
mag = []
err = []
cepid = []
cepname = []
vmag = []
imag = []
ident = []


## SMC SPECIFIC
## Want to neglect and VI data in the avs files if OGLE data is available.
## OGLE did not go through gloess -- is in a separate file

## Read in the OGLE data first. If it has m=99.999 then can overwrite that mag with 
## the gloess one if it's available.

## Need to add something to fix it if no V or I mag available so the 99s don't
## mess up the averages.

for line in open("ogle_mags","r"):
	data = line.split()
	cepname.append(data[0])
	imag.append(float(data[1]))
	vmag.append(float(data[2]))


for name in files:
	cepid.append(re.sub('.glo_avs', '', name))
	#print cepid
	for line in open(name,'r'):
		data = line.split()
		if data[3] != '=' and data[3] != 'std':
		## Need to remove the extra characters with regex -- re module
			ident.append(re.sub('.glo_avs', '', name))
			band.append(re.sub('[\[\]<>]', '', data[3]))
			period.append(float(data[2]))
			mag.append(float(data[5]))
			if data[7] != 'single':
				err.append(float(data[9]))
			else:
				err.append(0.1)

ident = np.array(ident)
cepname = np.array(cepname)
cepid = np.array(cepid)
imag = np.array(imag)
vmag = np.array(vmag)				
period = np.array(period)
mag = np.array(mag)
err = np.array(err)
band = np.array(band)

tident = []
tband = []
tperiod = []
tmag = []
terr = []

## Want to check to see if the data exists, overwrite if it does, append if not

for count in range(0,cepname.size):
	if imag[count] < 99.:
		if mag[(ident==cepname[count]) & (band=='I')].size == 1:
			mag[(ident==cepname[count]) & (band=='I')] = imag[count]
			err[(ident==cepname[count]) & (band=='I')] = 0.1
		if mag[(ident==cepname[count]) & (band=='I')].size == 0:
			tident.append(cepname[count])
			tband.append("I")
			tperiod.append(period[(ident==cepname[count]) & (band=='3.6')])
			tmag.append(imag[count])
			terr.append(0.1)
	if vmag[count] < 99.:
		if mag[(ident==cepname[count]) & (band=='V')].size == 1:
			mag[(ident==cepname[count]) & (band=='V')] = vmag[count]
			err[(ident==cepname[count]) & (band=='V')] = 0.1
		if mag[(ident==cepname[count]) & (band=='V')].size == 0:
			tident.append(cepname[count])
			tband.append("V")
			tperiod.append(period[(ident==cepname[count]) & (band=='3.6')])
			tmag.append(vmag[count])
			terr.append(0.1)
			
tident = np.array(tident)
tband = np.array(tband)
tperiod = np.array(tperiod)
tmag = np.array(tmag)
terr = np.array(terr)

newident = np.append(ident, tident)
newband = np.append(band, tband)
newperiod = np.append(period, tperiod)
newmag = np.append(mag, tmag)
newerr = np.append(err, terr)



ident = newident
band = newband
period = newperiod
mag = newmag
err = newerr

## At this point should have replaced old data with OGLE data, added new data where it
## wasn't available

## IT WORKS! Using OGLE mags for V and I where they exist :)

logp = log10(period)

## Fit the PL relations
## B - I Fouque 2007
## JHK Persson
## IRAC CHP 

def plu(lp, slopeu, zpu):
	return slopeu*(lp - 1.0) + zpu
def plb(lp, zpb):
	return -2.393 * (lp - 1.0) + zpb
def plv(lp, zpv):
	return -2.734 * (lp - 1.0) + zpv
def plr(lp, zpr):
	return -2.742 * (lp - 1.0) + zpr
def pli(lp, zpi):
	return -2.957 * (lp - 1.0) + zpi
def plj(lp, zpj):
	return -3.153 * (lp - 1.0) + zpj
def plh(lp, zph):
	return -3.234 * (lp - 1.0) + zph
def plk(lp, zpk):
	return -3.281 * (lp - 1.0) + zpk
def pl1(lp, zp1):
	return -3.306 * (lp - 1.0) + zp1
def pl2(lp, zp2):
	return -3.207 * (lp - 1.0) + zp2
	
#lpb = logp[band=='B']
#magb = mag[band=='B']
## U band fitting both coefficients
popt, pcov = curve_fit(plu,logp[(band=='U') & (logp<2.0) & (logp>0.77815) ],mag[(band=='U') & (logp< 2.0) & (logp>0.77815)])
slopeu = popt[0]
zpu = popt[1]
eslopeu = sqrt(float(pcov[0][0]))
ezpu = sqrt(float(pcov[1][1]))

popt, pcov = curve_fit(plb,logp[(band=='B') & (logp<2.0) & (logp>0.77815) ],mag[(band=='B') & (logp< 2.0) & (logp>0.77815)])
zpb = popt[0]
ezpb = sqrt(float(pcov[0]))
popt, pcov = curve_fit(plv,logp[(band=='V') & (logp<2.0) & (logp>0.77815)],mag[(band=='V') & (logp< 2.0) & (logp>0.77815)])
zpv = popt[0]
ezpv = sqrt(float(pcov[0]))
popt, pcov = curve_fit(plr,logp[(band=='R') & (logp<2.0) & (logp>0.77815)],mag[(band=='R') & (logp< 2.0) & (logp>0.77815)])
zpr = popt[0]
ezpr = sqrt(float(pcov[0]))
popt, pcov = curve_fit(pli,logp[(band=='I') & (logp<2.0) & (logp>0.77815)],mag[(band=='I') & (logp< 2.0) & (logp>0.77815)])
zpi = popt[0]
ezpi = sqrt(float(pcov[0]))
popt, pcov = curve_fit(plj,logp[(band=='J') & (logp<2.0) & (logp>0.77815)],mag[(band=='J') & (logp< 2.0) & (logp>0.77815)])
zpj = popt[0]
ezpj = sqrt(float(pcov[0]))
popt, pcov = curve_fit(plh,logp[(band=='H') & (logp<2.0) & (logp>0.77815)],mag[(band=='H') & (logp< 2.0) & (logp>0.77815)])
zph = popt[0]
ezph = sqrt(float(pcov[0]))
popt, pcov = curve_fit(plk,logp[(band=='K') & (logp<2.0) & (logp>0.77815)],mag[(band=='K') & (logp< 2.0) & (logp>0.77815)])
zpk = popt[0]
ezpk = sqrt(float(pcov[0]))
popt, pcov = curve_fit(pl1,logp[(band=='3.6') & (logp<2.0) & (logp>0.77815)],mag[(band=='3.6') & (logp< 2.0) & (logp>0.77815)])
zp1 = popt[0]
ezp1 = sqrt(float(pcov[0]))
popt, pcov = curve_fit(pl2,logp[(band=='4.5') & (logp<2.0) & (logp>0.77815)],mag[(band=='4.5') & (logp< 2.0) & (logp>0.77815)])
zp2 = popt[0]
ezp2 = sqrt(float(pcov[0]))

## Setting up fits for plotting

lp1 = np.arange(0,3.0,0.1)

plufit = slopeu*(lp1 - 1.0) + zpu
plbfit = -2.393*(lp1 - 1.0) + zpb
plvfit = -2.734 * (lp1 - 1.0) + zpv
plrfit = -2.742 * (lp1 - 1.0) + zpr
plifit = -2.957 * (lp1 - 1.0) + zpi
pljfit = -3.153 * (lp1 - 1.0) + zpj
plhfit = -3.234 * (lp1 - 1.0) + zph
plkfit = -3.281 * (lp1 - 1.0) + zpk
pl1fit = -3.306 * (lp1 - 1.0) + zp1
pl2fit = -3.207 * (lp1 - 1.0) + zp2

## Calculate the standard deviations

deviationsu = plb(logp[(band=='U') & (logp<2.0) & (logp>0.77815)],zpb) - mag[(band=='U') & (logp< 2.0) & (logp>0.77815)]
devusq = deviationsu**2
sdevu = sqrt((1./float(size(mag[(band=='U') & (logp< 2.0) & (logp>0.77815)])))*sum(devusq))

deviationsb = plb(logp[(band=='B') & (logp<2.0) & (logp>0.77815)],zpb) - mag[(band=='B') & (logp< 2.0) & (logp>0.77815)]
devbsq = deviationsb**2
sdevb = sqrt((1./float(size(mag[(band=='B') & (logp< 2.0) & (logp>0.77815)])))*sum(devbsq))
deviationsv = plv(logp[(band=='V') & (logp<2.0) & (logp>0.77815)],zpv) - mag[(band=='V') & (logp< 2.0) & (logp>0.77815)]
devvsq = deviationsv**2
sdevv = sqrt((1./float(size(mag[(band=='V') & (logp< 2.0) & (logp>0.77815)])))*sum(devvsq))
deviationsr = plr(logp[(band=='R') & (logp<2.0) & (logp>0.77815)],zpr) - mag[(band=='R') & (logp< 2.0) & (logp>0.77815)]
devrsq = deviationsr**2
sdevr = sqrt((1./float(size(mag[(band=='R') & (logp< 2.0) & (logp>0.77815)])))*sum(devrsq))
deviationsi = pli(logp[(band=='I') & (logp<2.0) & (logp>0.77815)],zpi) - mag[(band=='I') & (logp< 2.0) & (logp>0.77815)]
devisq = deviationsi**2
sdevi = sqrt((1./float(size(mag[(band=='I') & (logp< 2.0) & (logp>0.77815)])))*sum(devisq))
deviationsj = plj(logp[(band=='J') & (logp<2.0) & (logp>0.77815)],zpj) - mag[(band=='J') & (logp< 2.0) & (logp>0.77815)]
devjsq = deviationsj**2
sdevj = sqrt((1./float(size(mag[(band=='J') & (logp< 2.0) & (logp>0.77815)])))*sum(devjsq))
deviationsh = plh(logp[(band=='H') & (logp<2.0) & (logp>0.77815)],zph) - mag[(band=='H') & (logp< 2.0) & (logp>0.77815)]
devhsq = deviationsh**2
sdevh = sqrt((1./float(size(mag[(band=='H') & (logp< 2.0) & (logp>0.77815)])))*sum(devhsq))
deviationsk = plk(logp[(band=='K') & (logp<2.0) & (logp>0.77815)],zpk) - mag[(band=='K') & (logp< 2.0) & (logp>0.77815)]
devksq = deviationsk**2
sdevk = sqrt((1./float(size(mag[(band=='K') & (logp< 2.0) & (logp>0.77815)])))*sum(devksq))
deviations1 = pl1(logp[(band=='3.6') & (logp<2.0) & (logp>0.77815)],zp1) - mag[(band=='3.6') & (logp< 2.0) & (logp>0.77815)]
dev1sq = deviations1**2
sdev1 = sqrt((1./float(size(mag[(band=='3.6') & (logp< 2.0) & (logp>0.77815)])))*sum(dev1sq))
deviations2 = pl2(logp[(band=='4.5') & (logp<2.0) & (logp>0.77815)],zp2) - mag[(band=='4.5') & (logp< 2.0) & (logp>0.77815)]
dev2sq = deviations2**2
sdev2 = sqrt((1./float(size(mag[(band=='4.5') & (logp< 2.0) & (logp>0.77815)])))*sum(dev2sq))

mub = zpb - 17.356 + 18.48 + 2.393
muv = zpv - 17.052 + 18.48 + 2.734 
mur = zpr - 16.697 + 18.48 + 2.742 
mui = zpi - 16.589 + 18.48 + 2.957
muj = zpj - 16.336 + 18.48 + 3.153
muh = zph - 16.079 + 18.48 + 3.234
muk = zpk - 16.051 + 18.48 + 3.281
mu1 = zp1 -  -5.80
mu2 = zp2 - -5.77

emub = sqrt((ezpb**2) + (0.01**2))
emuv = sqrt((ezpv**2) + (0.007**2))
emur = sqrt((ezpr**2) + (0.020**2))
emui = sqrt((ezpi**2) + (0.005**2))
emuj = sqrt((ezpj**2) + (0.015**2))
emuh = sqrt((ezph**2) + (0.012**2))
emuk = sqrt((ezpk**2) + (0.011**2))
emu1 = sqrt((ezp1**2) + (0.03**2))
emu2 = sqrt((ezp2**2) + (0.03**2))

## Plot the relations

mp.clf() 	
mp.figure(figsize=(10, 8))
		

uniq_ids = np.unique(ident)

wavelengths = float(0.445), float(0.551), float(0.658), float(0.806), float(1.22), float(1.63), float(2.19), float(3.545), float(4.442)
wavelengths = np.array(wavelengths)
invw = 1./wavelengths
tempz = np.zeros(1)
invw = np.concatenate((invw, tempz))

alam = []
Ak = red.ccm_nearir(2.164,3.1)
for wlen in range(0,4):
	alam.append(red.ccm_optical(wavelengths[wlen], 3.1))
for wlen in range(4,7):
	alam.append(red.ccm_nearir(wavelengths[wlen],3.1))
for wlen in range(7,9):
	alam.append(red.indebetouw_ir(wavelengths[wlen])*Ak)
alam.append(0.0)
## Added on the zero to the end to deal with infinite wavelength
alam = np.array(alam)


## Here I am defining what the change in mu would be expected from no redding, 
## just from the dispersions in the PL relations at different bands.

## Fitted a quadratic function to the LMC dispersions (the canonical fits)
## to see how the 1 sigma dispersion changes as a function of wavelength.

plus1sig = 0.3861*invw**2 -0.0166*invw + 0.1105
minus1sig = -plus1sig

fitted = 0, 1, 0, 1, 0, 0, 0, 1, 0
fitted = np.array(fitted)

plotwhite = 1, 0, 1, 0, 0, 0, 0, 0, 1
plotwhite = np.array(plotwhite)

## This is what we expect from no reddening, offset to mean distance modulus of galaxy
invw_t = np.arange(0,2.51,0.01)

plus1sig = 0.03861*invw_t**2 -0.0166*invw_t + 2*0.1105 + 18.96
minus1sig = -0.03861*invw_t**2 +0.0166*invw_t - 2*0.1105 + 18.96


## This section is for a single plot containing all stars. Comment out if you want
## a plot for each star


mp.clf()
ax2 = subplot(111)

	
myaxis1 = [0.0,2.25,18.45,20.25]
mp.axis(myaxis1)
	
#title_text = "SMC Cepheids Reddenings"
#mp.title(title_text)

mp.xlabel('1 / $\lambda (\mu m^{-1})$')
mp.ylabel('$\mu$')

## End of single plot definitions

for cepheid in uniq_ids:
	if mag[(ident==cepheid) & (band=='B')].size == 1:
		vmub =  -17.356 + 18.48 + 2.393*logp[(band=='B') & (ident==cepheid)] + mag[(band=='B') & (ident==cepheid)]
		errmub = sqrt(err[(ident==cepheid) & (band=='B')]**2 + (ezpb)**2)
	else:
		vmub = 0.00
		errmub = 0.00
	if mag[(ident==cepheid) & (band=='V')].size == 1:	
		vmuv =  - 17.052 + 18.48 + 2.734*logp[(band=='V') & (ident==cepheid)] + mag[(band=='V') & (ident==cepheid)]
		errmuv = sqrt(err[(ident==cepheid) & (band=='V')]**2 + (ezpv)**2)
	else:
		vmuv = 0.00
		errmuv = 0.00

	if mag[(ident==cepheid) & (band=='R')].size == 1:	
		vmur =  - 16.697 + 18.48 + 2.742*logp[(band=='R') & (ident==cepheid)] + mag[(band=='R') & (ident==cepheid)]
		errmur = sqrt(err[(ident==cepheid) & (band=='R')]**2 + (ezpr)**2)

	else:
		vmur = 0.00
		errmur = 0.00

	if mag[(ident==cepheid) & (band=='I')].size == 1:	
		vmui =  - 16.589 + 18.48 + 2.957*logp[(band=='I') & (ident==cepheid)] + mag[(band=='I') & (ident==cepheid)]	
		errmui = sqrt(err[(ident==cepheid) & (band=='I')]**2 + (ezpi)**2)

	else:
		vmui = 0.00
		errmui = 0.00
	if mag[(ident==cepheid) & (band=='J')].size == 1:	
		vmuj =  - 16.336 + 18.48 + 3.153*logp[(band=='J') & (ident==cepheid)] + mag[(band=='J') & (ident==cepheid)]
		errmuj = sqrt(err[(ident==cepheid) & (band=='J')]**2 + (ezpj)**2)

	else:
		vmuj = 0.00
		errmuj = 0.00
	if mag[(ident==cepheid) & (band=='H')].size == 1:	
		vmuh =  - 16.079 + 18.48 + 3.234*logp[(band=='H') & (ident==cepheid)] + mag[(band=='H') & (ident==cepheid)]
		errmuh = sqrt(err[(ident==cepheid) & (band=='H')]**2 + (ezph)**2)

	else:
		vmuh = 0.00
		errmuh = 0.00
	if mag[(ident==cepheid) & (band=='K')].size == 1:	
		vmuk =  - 16.051 + 18.48 + 3.281*logp[(band=='K') & (ident==cepheid)] + mag[(band=='K') & (ident==cepheid)]
		errmuk = sqrt(err[(ident==cepheid) & (band=='K')]**2 + (ezpk)**2)

	else:
		vmuk = 0.00
		errmuk = 0.00
	if mag[(ident==cepheid) & (band=='3.6')].size == 1:	
		vmu1 =  5.80 + 3.306*(logp[(band=='3.6') & (ident==cepheid)] - 1.) + mag[(band=='3.6') & (ident==cepheid)]
		errmu1 = sqrt(err[(ident==cepheid) & (band=='3.6')]**2 + (ezp1)**2)
	else:
		vmu1 = 0.00
		errmu1 = 0.00
	# removing r and 4.5 from the fit
	if mag[(ident==cepheid) & (band=='4.5')].size == 1:	
		vmu2 =  5.77 + 3.207*(logp[(band=='4.5') & (ident==cepheid)]  - 1.0)+ mag[(band=='4.5') & (ident==cepheid)]
		errmu2 = sqrt(err[(ident==cepheid) & (band=='4.5')]**2 + (ezp2)**2)

	else:
		vmu2 = 0.00
		errmu2 = 0.00
	# v denotes value
	print cepheid, float(vmub), ",", float(vmuv), ",", float(vmur), ",", float(vmui), ",", float(vmuj), ",", float(vmuh), ",", float(vmuk), ",", float(vmu1), ",", float(vmu2)
	#muo, Av,muoerr, Averr = red.fit_reddening(float(vmub), float(vmuv), float(vmur), float(vmui), float(vmuj), float(vmuh), float(vmuk), float(vmu1), float(vmu2))
	
	##Only fitting V,I,[3.6] but including the others in the plot
	
	## Here I'm making 1 plot per Cepheid
	

	 
	if period[(ident==cepheid) & (band=='3.6') ] > 6. and period[(ident==cepheid) & (band=='3.6')] < 60.:
		muo, Av,muoerr, Averr = red.fit_reddening(0.00, float(vmuv), 0.00, float(vmui), 0.00, 0.00, 0.00, float(vmu1), 0.00)
	
	
		mus = float(vmub), float(vmuv), float(vmur), float(vmui), float(vmuj), float(vmuh), float(vmuk), float(vmu1), float(vmu2), muo
		errs = float(errmub), float(errmuv), float(errmur), float(errmui), float(errmuj), float(errmuh), float(errmuk), float(errmu1), float(errmu2), 
	
		mucor = mus - (Av*alam)
		
		mus = np.array(mus)
		errs = np.array(errs)
		
	
	
	
	## Here I am defining what the change in mu would be expected from no redding, 
	## just from the dispersions in the PL relations at different bands.

	## Fitted a quadratic function to the LMC dispersions (the canonical fits)
	## to see how the 1 sigma dispersion changes as a function of wavelength.


		mus = np.array(mus)
		mugood = []
	

	
		print invw
		Atrue = (alam * Av) + muo
		mugood = mus - (alam*Av)
		print Atrue
		Atrue = np.array(Atrue)
		ax2.plot(invw,Atrue,color='DimGray',ls='-')






ax3 = ax2.twiny()
ax3.axis(myaxis1)

## make the ticks at the band passes
ax3.set_xticks([1.0/0.445, 1.0/0.551, 1.0/0.658, 1.0/0.806, 1.0/1.22, 1.0/1.63, 1.0/2.19, 1.0/3.545, 1.0/4.442])
plt.setp(ax3.get_xticklabels(), visible=False)
## make the labels as the band passes
ax4 = ax2.twiny()
ax4.axis(myaxis1)

ax4.set_xticks([1.0/0.445, 1.0/0.551, 1.0/0.658, 1.0/0.806, 1.0/1.22, 1.0/1.63, 1.0/2.19, 1.0/3.35, 1.0/4.65])
ax4.tick_params(axis='x', which='both', top='off')
ax4.set_xticklabels(['$B$', '$V$', '$R$', '$I$', '$J$', '$H$', '$K$', '[3.6]', '[4.5]'])
ax2.fill_between(invw_t, minus1sig, plus1sig, facecolor='Turquoise', alpha=0.5, label='Instability Strip')
#ax2.plot(invw_t, plus1sig, 'r-',lw=2)
#ax2.plot(invw_t, minus1sig, 'DodgerBlue',lw=2)
ax2.plot(invw_t, plus1sig, color='Turquoise',lw=4, label='Instability Strip')
ax2.plot(invw_t, minus1sig, color='Turquoise',lw=4)

## Add the mean reddening line to the plot

## mean reddening is Av = 0.2338 

Atrue = (alam * 0.2201) + 18.96
mugood = mus - (alam*2201)
Atrue = np.array(Atrue)
ax2.plot(invw,Atrue,color='Crimson',ls='-', lw=4, label='$E(B-V) = 0.071$')

handles, labels = ax2.get_legend_handles_labels() 
ax2.legend(handles[::-1],labels[::-1],loc=2, numpoints=1,prop={'size':14}, fancybox=True)
mp.show()

output_r = 'all_reddenings.pdf'
mp.savefig(output_r,transparent=True)
#output_midir.close()
output_zps.close()
#output_reds.close()



