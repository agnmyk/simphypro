#!/usr/bin/env python
from __future__ import division
import re
import os
import sys
from glob import glob
from getopt import getopt

blastdb=""
gi2taxid_file=""
do_psiblast=True
htaxid={}
taxhash={}
numOfIter=3
out=False
rm_psiblast=False
psi_path=""


try:
	opts, args = getopt(sys.argv[1:], "hd:g:p:i:orl:", ["help", "database=", "gi2taxid_file=", "psiblast=", "num_of_psiblast_iteations=", "output", "remove_psiblast_files", "psiblast_path="])
except getopt.GetoptError:          
	usage()                        
	sys.exit(2)

for opt, arg in opts:               
	if opt in ("-h", "--help"):     
        	usage()                     
        	sys.exit()  
	elif opt in ("-d", "--database"):
		blastdb = arg
	elif opt in ("-g", "--gi2taxid_file"):
		gi2taxid_file = arg	                
        elif opt in ("-p", "--psiblast"):
		if arg=="False" or arg=="false":	do_psiblast = False	
	elif opt in ("-i", "--num_of_psiblast_iteations"):
		numOfIter = int(arg)
	elif opt in ("-o", "--output"):
		out=True
	elif opt in ("-r", "--remove_psiblast_files"):
		rm_psiblast=True
	elif opt in ("-l", "--psiblast_path"):
		psi_path=arg


def usage():
	print "\nusage: ./psi_parse.py [options] fasta_dir\n"
	print "-d, --database=path\n\tpath to protein database"
	print "-g, --gi2taxid_file=filename\n\tname of gi_to_taxid file including path"
	print "-p, --psiblast=boolean\n\tcontrols whether the psiblast will be launched.\n\tIf the psiblast files are already prepared it should be set to False (default is True)"
	print "-i, --pnum_of_psiblast_iteations=int\n\tnumber of psiblast iterations (default set to 3)"
	print "-o, --output\n\tif present files names and taxidlist will be produced\n"
	print "-r, --remove_psiblast_files\n\tif present results of psiblast will be removed\n"
	print "-l, --psiblast_path\n\tif the global path for the psiblast in the system is not set, please specify it in this option\n"
	sys.exit()
	

def parse(fastafile, htaxid, names, do_psiblast, out, rm_psiblast):
	if do_psiblast:
		print >> sys.stderr, "PSI_BLAST for %s...\n" %(fastafile)
		os.system("%spsiblast -num_threads 9 -num_iterations %d -db %s -outfmt 5 -evalue 0.001 -query %s -out %s.psiblast" %(psi_path, numOfIter,  blastdb, fastafile, fastafile))
	file=open("%s.psiblast" %(fastafile))
	ftaxid=open("%s.taxid" %(fastafile), 'w')
	t=file.read()
	t=t.split('\n')

	for i in xrange(len(t)-1,-1,-1):
		if re.search("<Iteration_iter-num>", t[i])!=None:
			iterstart=i
			break
	query=t[iterstart+2]
	query=re.search("[a-z]{2}\|[^\|]*\|[^ \t\n\r\f\v]+",query)
	if query!=None:			
		query=query.group() 
		print >> sys.stderr,"Query:", query
		if out:
			names.write(query+"\n")
		taxhash={}
		querylength=0
		for i in xrange(iterstart, len(t)):
			if re.search("<Hit_num>", t[i])!=None:
				hitname=re.search("<Hit_id>([^<]*)</Hit_id>",t[i+1]).group(1) 
				for j in xrange(i, len(t)):
					hitlength=re.search("<Hsp_align-len>",t[j])
					if hitlength!=None:
						hitlength=int(re.search("<Hsp_align-len>([^<]*)</Hsp_align-len>",t[j]).group(1)) 
						i=j
						break
				if querylength==0:	querylength=hitlength
				if querylength/hitlength < 0.8 or querylength/hitlength > 1.2: continue
				gi=int(re.search("gi\|(\d+)\|", hitname).group(1))
				if not htaxid.has_key(gi): htaxid[gi]=1 
				htaxidvar=htaxid[gi]
				if not taxhash.has_key(htaxidvar):
					taxhash[htaxidvar]=1
					ftaxid.write(str(htaxidvar)+"\n")
	file.close()
	ftaxid.close()
	if rm_psiblast:
		os.popen("rm %s.psiblast" %(fastafile))

indir=args
filelist=sorted(glob('%s/*[0-9]' %(indir[0])))
if len(filelist)==0:
	filelist=sorted(glob('%s/*.fasta' %(indir)))
	if len(filelist)==0:
		print >> sys.stderr, "no fasta file in this directory\n"
		usage()
		sys.exit()
	for i in xrange(len(filelist)):
		fr=re.search("(%s/*).fasta", filelist[i])
		if fr==None:
			print >> sys.stderr, "wrong file name: ", filelist[i]	
			sys.exit()
		filelist[i]=fr.group(1)	

if out:	names=open("names", "w")
else:	names=""


print >> sys.stderr, "importing gi_taxid dump file\n\n"

try:
	file=open(gi2taxid_file, 'r')
except IOError:
	print gi2taxid_file, "not found.\n\n"
	sys.exit()

content = file.readline()
while (content !=''):
	content=content.rstrip()
	gin, taxid=content.split()
	htaxid[int(gin)]=int(taxid)		
	content = file.readline()
file.close()


#run blast
print >> sys.stderr, "running blast\n"

for fastafile in filelist:
	parse(fastafile, htaxid, names, do_psiblast, out, rm_psiblast)
if out:
	os.popen("cat %s/*.taxid | sort | uniq > taxidlist" %(indir))
	names.close()
	
