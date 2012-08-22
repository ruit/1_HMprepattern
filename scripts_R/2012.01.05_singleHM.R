#Tian R., Zhang Y., Dec.8, 2011
#Jan 05, 2012
#systematically assess the discrimative abilities of the feature combinations

setwd("/Users/tianr/Desktop/Jan5.2012.Redcell/Y2012_Jan05_DataCleaning/SingleHM_model")

source("SVMmodule_RT_Jan05_2012.R")

posi<-read.table("Posi.HSC.cleaned.TAB", header=T)
nega<-read.table("Nega.HSC.cleaned.TAB", header=T)

###((()))numbers must be equal!!!

##for single HM

report<-c()

each<-c()

for (i in 1:100){
	
	each<-Summa_nega (induced=subset(posi, select=c(gene, H2AZ_Pro)), SS=subset(nega, select=c(gene, H2AZ_Pro)))
	report<-rbind(report, each)
	print (dim(report))
}


write.table(report, file='H2AZ_Pro.model', sep='\t')



#for double HMs

#pf<-Summa_nega (induced=subset(posi, select=c(gene, H2AZ_Body1, H2AZ_Body2)), SS=subset(nega, select=c(gene, H2AZ_Body1, H2AZ_Body2)))

#write.table(pf, file='H2AZ_Body1.H2AZ_Body2.perf.report', sep='\t')


###############################
#for tripple HMs
#pf<-Summa_nega (induced=subset(posi, select=c(gene, H2AZ_BodyA, H2AZ_BodyB, H2AZ_BodyC)), SS=subset(nega, select=c(gene, H2AZ_BodyA, H2AZ_BodyB, H2AZ_BodyC))) 
#write.table(pf, file='H2AZ_BodyA.H2AZ_BodyB.H2AZ_BodyC.perf.report', sep='\t')
