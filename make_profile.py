#!/usr/bin/env python
from __future__ import division
import sys
import os
from glob import glob
import re
import numpy
from bitarray import bitarray
from getopt import getopt

names=[]
taxidlist={}


numOfIter=3
params=sys.argv
bottom=0.6
top=1.4

def usage():
	print "usage: ./make_profile.py [options] fasta_dir\n"   
	print "-t, --taxcount_limit=float from [0;1]\n\tsimilarity coefficients are counted only if a proportion of taxons list lengths is between taxcount_limit and 2-taxcount_limit"

try:
	opts, args = getopt(sys.argv[1:], "ht:", ["help", "taxcount_limit="])
except getopt.GetoptError:          
	usage()                        
	sys.exit(2)

for opt, arg in opts:   
	if opt in ("-h", "--help"):     
        	usage()               
        	sys.exit()              
	if opt in ("-t", "--taxcount_limit"):     
        	bottom=float(arg)      



if args==[]:
	usage()
	sys.exit()

fresults=open("results", "w")
indir=args[0]
inputpattern="ls %s/*[0-9]|wc -l" %(indir)
seqcount=os.popen(inputpattern)
seqcount=int(seqcount.readlines()[0].rstrip())
seqcount-=1 #matrix index

try:
	m=open("taxidlist")
	file=open("names", 'r')
except:
	print "there is no taxidlist or names file, please generate it with 'get_names_and_taxidlist.py'"
	sys.exit()

#reading names file
n = file.readline()
while (n !=''):
	n=n.rstrip()
	names.append(n)		
	n = file.readline()
file.close()

#reading taxidlist file and preparing dictionary 
#which will be used to translate bit positions (in protein bitmaps) to taxon numbers for each protein
t = m.readline()
i=0
while (t !=''):
	t=int(t.rstrip())
	taxidlist[t]=i
	i+=1	
	t = m.readline()
m.close()
numoftaxid=i #bitmaps length 


taxidmap=[] #list with all protein bitmaps
counter=[] #list with numbers of taxons in which each protein was found

#preparing protein bitmaps
#protstmp contains for each protein a table with information if it was found in the subsequent taxons from the 'taxons' array
for taxidfile in sorted(glob('%s/*.taxid' %(indir))):
	taxtmp=bitarray(numoftaxid) #single protein bitmap
	taxtmp.setall(False)
	t=open("%s" %(taxidfile))
	tmp = t.readline()
	count=0
	while (tmp !=''):
		tmp=int(tmp.rstrip())
		if taxidlist.has_key(tmp):
			taxtmp[taxidlist[tmp]]=True
		tmp = t.readline()
		count+=1
	t.close()
	taxidmap.append(taxtmp)
	counter.append(count)

#score calculating
for i in xrange(seqcount):
	line1=counter[i]
	name1=names[i]
	taxidmap1=taxidmap[i]
	if line1<1: continue
	for j in xrange(i+1,seqcount+1):
		discore=(taxidmap1&taxidmap[j]).count()
		line2=counter[j]
		if discore==0: discore=1
		if line2<1: continue

		if line1<=line2:
			if line1>=line2*bottom:
				fresults.write("%s %s %15.13f\n" %(name1, names[j], (discore*discore)/(line1*line2)))
		else:
			if line1<=line2*top:
				fresults.write("%s %s %15.13f\n" %(name1, names[j], (discore*discore)/(line1*line2)))



