#TianR., Dec.13, 2011
#Plot the performances of models using all possible triple combinations of HM features

#for triple marks



#hm<-c("H2AZ_Body","H3K27me1_Body","H3K27me3_Body","H3K36me3_Body",
		"H3K4me1_Body","H3K4me3_Pro","H3K9me1_Body","H3K9me3_Body","H4K20me1_Pro")
Perf_summary<-function(file="triple.hms", sufix=".PER.tab"){
	
triple<-read.table(file)	
hm<-triple[,1]
df<-c()
Sensi<-c()
Spe<-c()
Accu<-c()
for (i in 1:length(hm)){


	file<-paste(as.character(hm[i]), sufix, sep="")
	df<-read.table(file, header=T)###!!!! Caution, this data has a header!!! Dec.13, 2011
#head(df)
	Sensi<-cbind(Sensi, df$Sensitivity)
	Spe<-cbind(Spe, df$Specificity)
	Accu<-cbind(Accu, df$Accuracy)
}

colnames(Sensi)<-hm
colnames(Spe)<-hm
	
return (list(Sensi=Sensi, Spe=Spe))
}


#################################################################


Sensi<-Perf_summary(file="double.hms", sufix=".PERF")$Sensi

Spe<-Perf_summary(file="double.hms", sufix=".PERF")$Spe


###############################################################
###define the margins"http://rgraphics.limnology.wisc.edu/rmargins_sf.php"

Plot_multi_box<-function(Sensi="min'gan",Spe="te'yi" ){
par(mar=c(8,4,5,2) + 0.1) 

boxplot(Sensi, col="yellow", ylim=c(30,80), las=3, ylab="Percentage (%)", 
	    cex.axis=0.9,
		main="Sensitivities VS Specificities")

grid<-seq(50,80,by=5)
for (j in 1:length(grid)){
abline(h=grid[j], lty=2, col="grey", lwd=0.5)}

legend(7,40,c("Sensitivity", "Specificity"), col=c("yellow", "green"), pch=19, cex=1.2)

boxplot(Spe, col="green", ylim=c(30,80), las=3,  add=T,  cex.axis=0.9)

for (i in seq(1,ncol(Spe),by=1)){
	text(i, 75, c(i), cex=0.8, col="purple")
	}
	
}

##################################################################

###select those triple combinations which has better trade off between Sensi and Spe

index<-c(13,18,19,24,40,44,46,55,59,61,65,71,76)

Plot_multi_box(Sensi=Perf_summary(hm[index])$Sensi, Spe=Perf_summary(hm[index])$Spe)


for (i in seq(1,length(index),by=1)){
	text(i, (45+1*i), c(paste(i,(as.character(hm[index[i]])), sep=" ")), cex=0.7, col="red")}

