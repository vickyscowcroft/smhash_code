#!/usr/bin/env python

import sys
import numpy as np
import re
import os

input = sys.argv[1]
outputName = "temp"

output = open(outputName,'w')
count = 0
for line in open(input,'r'):
 	if count < 3:
 		output.write("{0:s}".format(line))
 		count = count + 1
  		continue	

	## Check for bad stars that have *** in the lines
	hasStars = re.search("\*",line)
	if hasStars:
		count = count + 1
		continue
 	## Check for records with the wrong number of columns
 	data = line.split()
 	#print len(data)
 	if len(data)== 9:
 		output.write("{0:s}".format(line))
 		count = count + 1
 		continue
 	count = count + 1	
 
file.close(output)

os.rename("temp", sys.argv[1])
 
 		
 	
	