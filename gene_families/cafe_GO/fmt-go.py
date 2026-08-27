#!/usr/bin/env python

import sys
from collections import defaultdict

i=0
gid2go=defaultdict(set)

for l in open(sys.argv[1]):
	if l.startswith("#"): continue
	r=l.strip().split('\t')
	#print(r[0],r[9])
	gid=r[0].split('.')[0]
	if r[9]=='-': continue
	for go in r[9].split(','):
		gid2go[gid].add(go)
	i+=1
	#if i==8: break

with open(sys.argv[1].split('.')[0]+'.go.tsv','w') as out:
	out.write("term\tgene\n")
	for g in gid2go:
		for go in list(gid2go[g]):
			out.write("{}\t{}\n".format(go,g))

#print(gid2go)

