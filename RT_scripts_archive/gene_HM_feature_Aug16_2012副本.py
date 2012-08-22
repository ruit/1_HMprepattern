#!/usr/bin/python
#Tian R., Sep 6, 2011
#divide genes into bins,count how many short reads fall onto each bin.
#The final vector for each gene should be weighted by gene length and
#seqencing depth. RPKM, reads per kb per 10 million total reads
#----------------------------------------------------

#The HM prepattern  project (Sep,2011-Oct., 2012)
#Core code for HM feature generation
#@Time stamp: June 27, 2012
#@Aug 13, 2012
#@Aug 16, 2012
#-----------------------------------------------------
#Requirements for successful running of this programme
#BEDtools pre-installed and make sure all the commands from BEDtools are accessible from the current directory
#This has something to do with your PATH environment settings.
#-------------Only work in Linux/Mac systems-----------

############################################################################
import sys, re, os
def Makebins3(ref_file, cover_len=2000):
	'''
	Take care, gene body region might need to be re-defined. June 27, 2012 \ 
	#input: ref_file(format: refseq_geneID \t chrname \t strand \t TSS \t TTS)  \ 
	#name	chrom	strand	txStart	txEnd \ 
	#NM_032291	chr1	+	66999824	67210768 \ 
	#Divide gene regions into bins: 5' -2000bp, -1000bp, TSS to TTS divided into 10 regions of equal lengths, TTS +1000bp, TTS +2000bp \ 
	#Be cautious with genes on the "-" strand!
	#@Aug 13, 2012
	'''

	list_scale=[]#range
	dict_points={}
	gene_len=0
	fh_ref=open(ref_file, "r")
	bed_list=[]
	eachline=""
	for line in fh_ref:
		if re.search("NM_", line):
			line=line.rstrip()
			gene_ID, chr_name, strand, start, end=line.split("\t")
			
			if strand=="+":##### how to fetch 5' upstream and 3' downstream depends on the strand (+?-?)
				list_scale.append((int(start)-cover_len))# up and down 2kb (TSS)
				list_scale.append((int(start)+cover_len))
				list_scale.append(int(start))
				list_scale.append(int(end))
				dict_points[gene_ID]=list_scale
				for x in range(0,3,2):
					eachline=chr_name+"\t"+str(list_scale[x])+"\t"+str(list_scale[x+1])+"\t"+gene_ID+"_bin"+str(x+1)+"\t"+strand					
					bed_list.append(eachline)
				list_scale=[]
			
			if strand=="-":
				list_scale.append((int(end)+cover_len))# upstream ########!!!!!!!!!!!!!!!!##range, descendingbin
				list_scale.append((int(end)-cover_len))#downstream, following the same schema as postive strand
				list_scale.append(int(end))
				list_scale.append(int(start))
				dict_points[gene_ID]=list_scale
				for x in range(0,3,2):
					eachline=chr_name+"\t"+str(list_scale[x+1])+"\t"+str(list_scale[x])+"\t"+gene_ID+"_bin"+str(x+1)+"\t"+strand					
					bed_list.append(eachline)
				list_scale=[]
				
	return bed_list



#Read intersectBed output (-c) files (.cns) into NM_*****	x1	x2	...vector(14 elements)
def Read_intersectBed2(infile):
	dict_tag_counts={}
	array=[]
	list_counts=[]
	total=0
	frag_len=0.0
	nm, shuzi, part, vector="", "", "", ""
	for line in open(infile, "r"):
		line=line.strip()
		array=line.split("\t")
		gene_part=array[3]
		frag_len=float(int(array[2])-int(array[1]))/1000.0 ### fragment length is considered
		nm, shuzi, part=gene_part.split("_")
		gene_id=nm+"_"+shuzi####NM gene id
		
		#Give a pseudo count 1, to avoid log (0)
		#March 8, 2012
		if array[5]=="0":
			array[5]="1"
		else:
			pass
		total=str(float(array[5])/frag_len)###tag count in this part
		list_counts.append(total)
		if len(list_counts)==2:
			for i in range(1):
				vector+=list_counts[i]+"\t"
			vector+=list_counts[1]
			dict_tag_counts[gene_id]=vector
			vector=""
			list_counts=[]
	return dict_tag_counts ### attention! when printed, sort by key!



def main():
	import os
	pwd=""
	file_list=""
	dict={}
	bashscript=""
	bed_list=Makebins3(sys.argv[1])
	
	fh_out=open(sys.argv[1]+".bins", "w")
	fh_sh=open("temp.sh", "w")
	
	for i in range(len(bed_list)):
		
		fh_out.write(bed_list[i]+"\n")
	fh_out.close()
	
	#write bash script
	#bashscript="#"+"!/usr/bin/bash"+"\n"
	bashscript="for bedfile in `ls *bed`"+"\n"
	bashscript+="do"+"\n"
	bashscript+="\t"+"intersectBed -a " + "\""+sys.argv[1]+".bins" +"\"" + " -b " + "$bedfile" +" -c > " + "$bedfile" + "\".cns\""+"\n"
	bashscript+="done"+"\n"
	fh_sh.write(bashscript)
	fh_sh.close()
	os.system("chmod 755 temp.sh")
	os.system("./temp.sh")	
	pwd=os.getcwd()
	file_list=os.listdir(pwd)
        #fh_seqdepth=open("seqence_depth.txt", "w")
	for i in range(len(file_list)):
		if re.search(".cns$", file_list[i]):
			#print file_list[i]
			fh_cns=open(file_list[i]+".dst", "w")#Aug 16, 2012
			dict=Read_intersectBed2(file_list[i])
			#print dict
			for gene, vector in sorted(dict.iteritems()):
				fh_cns.write(gene+"\t"+vector+"\n")
			fh_cns.close()
	fh_seqdepth.close()

######

if __name__=="__main__":
	import sys
	from optparse import OptionParser
	usage='python .py refseq.file'
	parser=OptionParser(usage)
	#parser.add_option("-r", dest="file", help="The RefSeq annotation file")
	(options, args)=parser.parse_args()
	
	if len(sys.argv) < 2 :
	#	parser.print_help()
		sys.exit(usage)

	main()

