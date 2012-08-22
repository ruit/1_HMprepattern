#Tian R.,  March 20, 2012
#NaiveBayes Model for prediction of poised genes
#June 26, 2012, updated

hsc<-read.table("HSCdiffHM.2disc_March27_2012.txt", header=T)
#H2AZ_Pro H3K27me1_Body H3K27me3_Pro H3K36me3_Body H3K4me1_Body H3K4me3_Pro
#1        0             1            0             1            0           0
#2        1             0            2             1            0           1
#H3K9me1_Body H3K9me3_Pro H4K20me1_Body Poising
#1            1           0             0       1
#2            0           0             0       1
hsc.4hm<-subset(hsc, select=c(1,8,6,3, 10))
dim(hsc.4hm)
#	[1] 8107    5
head(hsc.4hm)
#	H2AZ_Pro H3K9me3_Pro H3K4me3_Pro H3K27me3_Pro Poising
#	1        0           0           0            0       1
#	2        1           0           1            2       1
hsc.4hm.f<-data.frame(H2AZ=factor(hsc.4hm$H2AZ_Pro),H3K9me3=factor(hsc.4hm$H3K9me3_Pro), H3K4me3=factor(hsc.4hm$H3K4me3_Pro), H3K27me3=factor(hsc.4hm$H3K27me3_Pro), poising=factor(hsc.4hm$Poising))
summary(hsc.4hm.f)
rm(hsc.4hm)
#	h2az     k9me3    k4me3    k27me3   poising 
#	0:4620   0:6246   0:4903   0:4412   0:7751  
#	1:2166   1:1658   1:2012   1:1248   1: 356  
#	2:1321   2: 203   2:1192   2:2447           

##############################################
library("e1071")

NBayes_h2az_k9me3<-function(hsc.4hm.f, times){
sensi<-c()##!!!!!!
speci<-c()
MCC<-c()

perf<-0.0
opt.model<-c()
test_data<-c()
for (num in 1:times){###March 28, 2012
	
	posi.tr.index<-sample(1:356,200)
	nega.tr.index<-sample(357:8107, 200)

	train<-hsc.4hm.f[c(posi.tr.index, nega.tr.index),]
	test<-hsc.4hm.f[-c(posi.tr.index, nega.tr.index),]

	new.model<-naiveBayes(poising ~ H2AZ + H3K9me3, data=train)
###MCC:Matthews correlation coefficient
###too big number, R can not handle
	TP<-0.01*table(predict(new.model, test[test$poising=="Lev1",][,1:2]))["Lev1"]
	FN<-0.01*table(predict(new.model, test[test$poising=="Lev1",][,1:2]))["Lev0"]
	TN<-0.01*table(predict(new.model, test[test$poising=="Lev0",][,1:2]))["Lev0"]
	FP<-0.01*table(predict(new.model, test[test$poising=="Lev0",][,1:2]))["Lev1"]
	if ((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))>perf){
		perf<-(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
		opt.model<-new.model
		test_data<-test}

#MCC<-c(MCC,(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
	
	sensi<-c(sensi, table(predict(new.model, test[test$poising=="Lev1",][,1:2]))["Lev1"]/nrow(test[test$poising=="Lev1",][,1:2]))
	speci<-c(speci, table(predict(new.model, test[test$poising=="Lev0",][,1:2]))["Lev0"]/nrow(test[test$poising=="Lev0",][,1:2]))
	
	}
#return (list(Sen=sensi, Spe=speci))
#return (MCC)
#return (perf)
return (list(opt.model=opt.model, test_data=test_data))
}


#############################
###select an optimal model, which better reflect the structure of the data, March 29, 2012
### ca 12 min ?
list<-NBayes_h2az_k9me3(hsc.4hm.f, 100)
new.model<-list$opt.model
test<-list$test_data
table(predict(new.model, test[test$poising=="Lev1",][,1:2]))["Lev1"]/nrow(test[test$poising=="Lev1",][,1:2])
table(predict(new.model, test[test$poising=="Lev0",][,1:2]))["Lev0"]/nrow(test[test$poising=="Lev0",][,1:2])

############################March 29, 2012
log_score<-function(di.model=model, df){
	lod<-c()
	for (i in 1:nrow(df)){
		p.posi<-di.model$tables$H2AZ["Lev1",df[i,]$H2AZ]*di.model$tables$H3K9me3["Lev1",df[i,]$H3K9me3]
		p.nega<-di.model$tables$H2AZ["Lev0",df[i,]$H2AZ]*di.model$tables$H3K9me3["Lev0",df[i,]$H3K9me3]
		if (p.nega==0){p.nega<-p.nega+0.000001}
	    lod<-c(lod,round(log2(p.posi/p.nega), 1))
	}
	return (lod)
}

#########################June 26, 2012
###make predictions based on HSC-derived NB classifer in more systems
Assess_Poising<-function (NE.di.file) {#input is a tab delimited file "gene\t H2AZ \t H3K9me3"
    #e.g. ""gm12878.2disc.non.exp.txt"
    gm<-read.table(NE.di.file, header=T)
#table(log_score(di.model=new.model, gm))
    decision.score<-log_score(di.model=new.model, gm)
    gm2<-data.frame(gm, score=decision.score)
	write.table(gm2, file="scored_NE_genes_from_di_HM.tab", sep="\t")
   }
#putative.poised.GM12878.NM<-gm2[gm2$score==0.9,]$gene
#putative.silenced.GM12878.NM<-gm2[gm2$score==-1.5,]$gene
#write.table(putative.poised.GM12878.NM, file="putative.poised.GM12878.NM", sep="\t")
#write.table(putative.silenced.GM12878.NM,file="putative.silenced.GM12878.NM", sep="\t")

############################
###Bernstein's model
NBayes_k4me3_k27me3<-function(hsc.4hm.f){
	sensi<-c()##!!!!!!
	speci<-c()
	for (num in 1:100){###March 28, 2012
		
		posi.tr.index<-sample(1:356,200)
		nega.tr.index<-sample(357:8107, 200)
		
		train<-hsc.4hm.f[c(posi.tr.index, nega.tr.index),]
		test<-hsc.4hm.f[-c(posi.tr.index, nega.tr.index),]
		
		bern.model<-naiveBayes(poising ~ k4me3 + k27me3 , data=train)
		sensi<-c(sensi,table(predict(bern.model, test[test$poising=="Lev1",][,3:4]))["Lev1"]/nrow(test[test$poising=="Lev1",]))#####3:4!!!!!!
		speci<-c(speci,table(predict(bern.model, test[test$poising=="Lev0",][,3:4]))["Lev0"]/nrow(test[test$poising=="Lev0",]))
		
	}
return (list(Sen=sensi, Spe=speci))
}

######### compare my model with B.B.'s bivalent model
perf.new<-NBayes_h2az_k9me3(hsc.4hm.f)
perf.bern<-NBayes_k4me3_k27me3(hsc.4hm.f)
##March 28, 2012
boxplot(list(perf.new$Sen, perf.bern$Sen, perf.new$Spe, perf.bern$Spe), ylim=c(0.45,0.68), col=c("red", "blue", "pink", "lightblue"))


########################################################################################
###K4+K9+K27 model
tri.model<-naiveBayes(poising ~  k9me3 +k4me3 +k27me3, data=train)
#sensi
table(predict(tri.model, test[test$poising=="Lev1",][,2:4]))[["Lev1"]]/nrow(test[test$poising=="Lev1",][,2:4])
#speci
table(predict(tri.model, test[test$poising=="Lev0",][,2:4]))[["Lev0"]]/nrow(test[test$poising=="Lev0",][,2:4])

###Bernstein's model:
bern.model<-naiveBayes(poising ~ k4me3 + k27me3 , data=train)
table(predict(bern.model, test[test$poising=="Lev1",][,3:4]))/nrow(test[test$poising=="Lev1",])
table(predict(bern.model, test[test$poising=="Lev0",][,3:4]))/nrow(test[test$poising=="Lev0",])

###K4+K9 model
k4k9.model<-naiveBayes(poising ~ k9me3 + k4me3 , data=train)
table(predict(k4k9.model, test[test$poising=="Lev1",][,2:3]))/nrow(test[test$poising=="Lev1",])
#86%, 0.8653846
table(predict(k4k9.model, test[test$poising=="Lev0",][,2:3]))/nrow(test[test$poising=="Lev0",])
#23%, 0.2340087
