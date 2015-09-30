#!/usr/bin/env python

import glob
import pexpect
import shutil
import os
import sys
import re
import subprocess

target = sys.argv[1]

channel = sys.argv[2]

print channel

targ_regex = re.compile(target)
if (channel == '1'):
	chan_regex = re.compile('3p6um')
	chan = '3p6um'
if(channel =='2'):
	chan_regex = re.compile('4.5um')
	chan = '4p5um'
if(channel != '1' and channel != '2'):
	print "Incorrect channel specified"
#		return (99)
	exit


file_list = glob.glob('*_e*.ap*')
file_list = filter(targ_regex.search, file_list) ## filters out files that don't correspond to this target
file_list = filter(chan_regex.search, file_list) ## filters out channels that are the wrong channel

## Need to also remove epoch 1 from the file list as it's special 

epoch1 = re.compile('_e1_')
epoch1_file = filter(epoch1.search, file_list)
epoch1_file = re.sub(".ap", ".als", epoch1_file[0])

old_epoch1_file = target + '_e1_' + chan + '_dn.ap'
if old_epoch1_file in file_list: file_list.remove(old_epoch1_file)

#file_list = append(epoch1_file, file_list)
#file_list
# Remove existing mch file

mch_file = target + '_' + chan + '.mch'

if (os.path.isfile(mch_file)):
		os.remove(mch_file)

## run daomatch
print "Running daomatch"

daomatch = pexpect.spawn("/Users/vs/daophot/daomatch")
daomatch.expect('Master input file:')
daomatch.sendline(epoch1_file)
daomatch.expect('Output file name')
daomatch.sendline(mch_file)
print "made the mch file"
## Loop over the file_list

for ap_file in file_list:
	daomatch.expect('Next input file')
	daomatch.sendline(ap_file)
	print "finished " + ap_file

## Exiting line
daomatch.expect('Next input file')
daomatch.sendline('')
daomatch.close(force=True)
#print daomatch.isalive()

## running daomaster
print "running daomaster"
num_frames = len(file_list) + 1
print num_frames

daomaster = pexpect.spawn("/Users/vs/daophot/daomaster")
daomaster.expect("File with list of input files")
daomaster.sendline(mch_file)
print mch_file
daomaster.expect("Minimum number, minimum fraction, enough frames")
daomaster.sendline(str(num_frames) + ", 1, " + str(num_frames))
print str(num_frames) + ", 1, " + str(num_frames)
daomaster.expect("Maximum sigma")
daomaster.sendline("99")
## desired degrees of freedom:
daomaster.expect("Your choice")
daomaster.sendline("6")
daomaster.expect("Critical match-up radius")
daomaster.sendline("10")

for radius in range (10,-1, -1):
	daomaster.expect("New match-up radius")
	daomaster.sendline(str(radius))

## Only need the transformations here
print "making output files"
daomaster.expect("Assign new star IDs")
print '1'
daomaster.sendline("n")
daomaster.expect("A file with mean magnitudes and scatter")
print '2'
daomaster.sendline("n")
daomaster.expect("A file with corrected magnitudes and errors")
print '3'
daomaster.sendline("n")
daomaster.expect("A file with raw magnitudes and errors")
print '4'
daomaster.sendline("n")
daomaster.expect("A file with the new transformations")
daomaster.sendline("y")
print "asked about transformations"
daomaster.expect("Output file name")
daomaster.sendline(mch_file + "_new")
print "asked about output file"
daomaster.expect("A file with the transfer table")
daomaster.sendline("e")
print "asked about transfer table"

daomaster.close(force=True)

## Fixing the bad output file name
new_mch  = glob.glob('*.mch_new*')
print new_mch
shutil.move( str(new_mch[0]), mch_file)









