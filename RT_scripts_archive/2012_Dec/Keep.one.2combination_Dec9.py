#!/usr/bin/python
#Tian R., Dec.9, 2011
#"H2AZ" "H4K20" two orders are essentially the same

import os, re

##for two histone marks
files=os.listdir(os.getcwd())

for i in range(len(files)):
	if re.search("report", files[i]):
		HMa,HMb,x,y=files[i].split(".")
		list=[HMa, HMb]
		list.sort()
		###double check this!
		newfile=list[0]+"."+list[1]+".PERF"
    		os.system("cat "+files[i]+" > "+newfile)

##################################
#files=os.listdir(os.getcwd())
#for tripple histone marks

#for i in range(len(files)):
#	if re.search("report", files[i]):
#		HMa,HMb, HMc, x,y=files[i].split(".")
#		list=[HMa, HMb, HMc]
#		list.sort()
		###double check this!
#		newfile=list[0]+"."+list[1]+"."+list[2]".PER"
#    		os.system("cat "+files[i]+" > "+newfile)
