#TianR., Dec.9, 2011
#Plot the performances of models using all possible single, double and triple combinations of HM features

#for single marks
hm<-c("H2AZ_Body","H3K27me1_Body","H3K27me3_Body","H3K36me3_Body",
		"H3K4me1_Body","H3K4me3_Pro","H3K9me1_Body","H3K9me3_Body","H4K20me1_Pro")

df<-c()
Sensi<-c()
Spe<-c()
Accu<-c()
for (i in 1:9){
#variable<-paste(hm[i], "_sensi", sep="")
	file<-paste(hm[i], ".perf.report", sep="")
	df<-read.table(file)
	Sensi<-cbind(Sensi, df$Sensitivity)
	Spe<-cbind(Spe, df$Specificity)
	Accu<-cbind(Accu, df$Accuracy)
}

colnames(Sensi)<-hm
colnames(Spe)<-hm

###define the margins"http://rgraphics.limnology.wisc.edu/rmargins_sf.php"
par(mar=c(8,4,5,2) + 0.1) 

boxplot(Sensi, col="yellow", ylim=c(30,80), las=3, ylab="Percentage (%)", 
	    cex.axis=0.9,
		main="Sensitivities VS Specificities")

grid<-seq(50,80,by=5)
for (j in 1:length(grid)){
abline(h=grid[j], lty=2, col="grey", lwd=0.5)}

legend(7,40,c("Sensitivity", "Specificity"), col=c("yellow", "green"), pch=19, cex=1.2)

boxplot(Spe, col="green", ylim=c(30,80), las=3,  add=T,  cex.axis=0.9)

