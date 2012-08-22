#Tian R., Dec.8, 2011
#systematically assess the discrimative abilities of the feature combinations

setwd("/Users/tianr/Desktop/Dec7.2011.Redcell/results_12.8.subset_features")

source("/Users/tianr/Desktop/Dec7.2011.Redcell/Rscripts/SVMmodule_RT_Dec6_2011.R")

#dataset<-Name_feature(induced.file="Redcell.Posi", SS.file="Redcell.Nega", fea.order="feature.order.redcell.copy")

#named.fea.posi<-dataset$posi
#named.fea.nega<-dataset$nega

posi<-read.table("named.9fea.posi.redcell.Dec8.2011")
nega<-read.table("named.9fea.nega.redcell.Dec8.2011")


###((()))numbers must be equal!!!

##for single HM
#pf<-Summa_nega (induced=subset(posi, select=c(gene, H2AZ_Body)), SS=subset(nega, select=c(gene, H2AZ_Body)))
#write.table(pf, file='H2AZ_Body.perf.report', sep='\t')



#for double HMs

#pf<-Summa_nega (induced=subset(posi, select=c(gene, H2AZ_Body1, H2AZ_Body2)), SS=subset(nega, select=c(gene, H2AZ_Body1, H2AZ_Body2)))

#write.table(pf, file='H2AZ_Body1.H2AZ_Body2.perf.report', sep='\t')


###############################
#for tripple HMs
pf<-Summa_nega (induced=subset(posi, select=c(gene, H4K20me1_Pro, H3K9me3_Body, H3K9me1_Body)), SS=subset(nega, select=c(gene, H4K20me1_Pro, H3K9me3_Body, H3K9me1_Body))) 
write.table(pf, file='H4K20me1_Pro.H3K9me3_Body.H3K9me1_Body.perf.report', sep='\t')
