#!/usr/bin/env python
import sys
import re
from glob import glob
import os

indir=sys.argv[1];
fout=open("names", "w")

for fastafile in sorted(glob('%s/*[0-9]' %(indir))):
	fin=open(fastafile, "r")
	p=fin.read().split("\n")
	for i in p:
		s=re.search("[a-z]{2}\|[^\|]*\|[^ \t\n\r\f\v]+",i)
		if s!=None:
			ss=s.group()
			fout.write(str(ss)+"\n")

os.popen("cat %s/*.taxid | sort | uniq > taxidlist" %(indir))
		

