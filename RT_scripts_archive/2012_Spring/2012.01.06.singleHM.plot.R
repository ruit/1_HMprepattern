#TianR., Dec.9, 2011
#Tian R., Zhang Y., Jan.6, 2012

#Plot the performances of models using all possible single, double and triple combinations of HM features

#for single marks
hm<-c("H2AZ_Pro","H3K27me1_Body","H3K27me3_Pro","H3K36me3_Body",
		"H3K4me1_Body","H3K4me3_Pro","H3K9me1_Body","H3K9me3_Pro","H4K20me1_Body")

df<-c()
Sensi<-c()
Spe<-c()
Accu<-c()

###
#compare full model with single HM models, Jan6, 2012
df<-read.table("/Users/tianr/Desktop/Jan5.2012.Redcell/Y2012_Jan05_DataCleaning/HSCmodel.9fea.reports_Jan5_2012.txt")

Sensi<-cbind(Sensi, df$Sensitivity)
Spe<-cbind(Spe, df$Specificity)
Accu<-cbind(Accu, df$Accuracy)

rm(df)

for (i in 1:9){
#variable<-paste(hm[i], "_sensi", sep="")
	file<-paste(hm[i], ".model", sep="")
	df<-read.table(file)
	Sensi<-cbind(Sensi, df$Sensitivity)
	Spe<-cbind(Spe, df$Specificity)
	Accu<-cbind(Accu, df$Accuracy)
}

colnames(Sensi)<-c("All",hm)
colnames(Spe)<-c("All",hm)

###define the margins"http://rgraphics.limnology.wisc.edu/rmargins_sf.php"
par(mar=c(8,4,5,2) + 0.1) 

boxplot(Sensi, col=rgb(255,0,0,200,maxColorValue=255), 
	ylim=c(20,90), las=3, ylab="Percentage (%)", 
	    cex.axis=0.9,
		main="Sensitivities VS Specificities")

grid<-seq(40,90,by=5)
for (j in 1:length(grid)){
abline(h=grid[j], lty=2, col="grey", lwd=0.5)}

boxplot(Spe, col=rgb(0,255,0,200,maxColorValue=255), 
	ylim=c(20,90), las=3,  add=T,  cex.axis=0.9)

legend(5,40,c("Sensitivity", "Specificity"), col=c(rgb(255,0,0,200,maxColorValue=255), rgb(0,255,0,200,maxColorValue=255)), pch=19, cex=1.2)
