#!/usr/bin/env python

import pexpect

print "Running daomatch"

epoch1_file = 'LINEAR_14385018_e1_3p6um_dn.als'

mch_file = 'LINEAR_14385018_3p6um.mch'

file_list = ['LINEAR_14385018_e10_3p6um_dn.ap', 'LINEAR_14385018_e11_3p6um_dn.ap', 'LINEAR_14385018_e12_3p6um_dn.ap', 'LINEAR_14385018_e2_3p6um_dn.ap', 'LINEAR_14385018_e3_3p6um_dn.ap', 'LINEAR_14385018_e4_3p6um_dn.ap', 'LINEAR_14385018_e5_3p6um_dn.ap', 'LINEAR_14385018_e6_3p6um_dn.ap', 'LINEAR_14385018_e7_3p6um_dn.ap', 'LINEAR_14385018_e8_3p6um_dn.ap', 'LINEAR_14385018_e9_3p6um_dn.ap']

daomatch = pexpect.spawn("daomatch")
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
