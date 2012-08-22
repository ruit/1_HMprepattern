#Tian R.
#Dec. 6, Tue., 2011
#Perform five fold cross validation SVM, taking account into consideration the excessive negative samples

#import core functions
source("SVMmodule_RT_Jan05_2012.R")

#assign column names for the features
#dataset<-Name_feature(induced.file="Redcell.Posi", SS.file="Redcell.Nega",fea.order="feature.order.redcell.copy" )

#repeat 100 times to get an idea about the negative samples

posi<-read.table("Posi.HSC.cleaned.TAB", header=T)
nega<-read.table("Nega.HSC.cleaned.TAB", header=T)


report<-c()

each<-c()

for (i in 1:100){
	each<-Summa_nega(induced=posi, SS=nega)
	report<-rbind(report, each)
	print (dim(report))
}

write.table(report, file="HSCmodel.9fea,reports.txt", sep="\t")

#empty all the variables in the Neicun
rm (list=ls(all=T))

