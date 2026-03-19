#!/usr/bin/env python

import sys,re

treline=open(sys.argv[1]).readline()

chdtre=re.sub('(\d+\.\d+)_(\d+\.\d+)',r'[&date_range={\1,\2}]',treline)
print(chdtre)
with open(sys.argv[1]+'.nex.tre','w') as out:
	out.write("#NEXUS\n\nbegin trees;\n\ttree tree_1 = [&R]\n{0}\n\nend;".format(chdtre))

#Replace : (\d+\.\d+)_(\d+\.\d+)
#By: [&date_range={$1,$2}]