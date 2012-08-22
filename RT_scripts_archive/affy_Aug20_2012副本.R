#Tian R.
#Aug. 20, 2012
#re-analyze HSC-ERY expression
#@cdfName is the core!!!
####################################
#install affy packages
#install all needed packages including "AnnotationDbi".
#source("http://bioconductor.org/biocLite.R") 
#biocLite("affy")

library("affy")
data<-ReadAffy(filenames=dir(pattern=".CEL"))
# get the affy cdf
data@cdfName
#
#[1] "HG-U133_Plus_2"

#download the right cdf.refseq file from Microarray lab
#http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/15.1.0/refseq.asp

#download to current directory
#in Linux terminal run
# R CMD INSTALL hgu133plus2hsrefseqcdf_15.0.0.tar.gz

#load and assign the right cdf file
library("hgu133plus2hsrefseqcdf")
data@cdfName<-"hgu133plus2hsrefseqcdf"

#P/M/A calling
pa<-mas5calls(data)

#RMA expression indexing
lev<-rma(data)

HSC_exp<-cbind(exprs(pa), exprs(lev))

write.table(HSC_exp, file="Expression_table_HSC_ERY_Aug20_2012.tab", sep="\t")

###
###select 4 groups of genes.
###poised pairs
silen2silen<-exp[exp[,1]=="A" & exp[,2]=="A" & exp[,3]=="A" & exp[,4]=="A", ]
silen2on<-exp[exp[,1]=="A" & exp[,2]=="A" & exp[,3]=="P" & exp[,4]=="P", ]

###switching off pairs
on2on<-exp[exp[,1]=="P" & exp[,2]=="P" & exp[,3]=="P" & exp[,4]=="P", ]
on2off<-exp[exp[,1]=="P" & exp[,2]=="P" & exp[,3]=="A" & exp[,4]=="A", ]

write.table(silen2silen, file="silen2silen.exprs.tab", sep="\t")
write.table(silen2on, file="silen2on.exprs.tab", sep="\t")
write.table(on2on, file="on2on.exprs.tab", sep="\t")
write.table(on2off, file="on2off.exprs.tab", sep="\t")

