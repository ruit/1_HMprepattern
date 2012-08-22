#Tian R & Yong Zhang
#Jan. 5, 2012
#Make sure the silenced genes belonging to positive and negative instances have low levels of H3K36me3.
#This is implemented by removing those with high levels of H3K36me3, and the majority of the data are kept.
#####################################################################################################

Name_log<-function(induced.fea="induced.fea", stay_silen.fea="SS.fea", fea_order="feature.order", pseudo=0.0001, take.log=TRUE){
	
	
	induced<-read.table(induced.fea)
	SS<-read.table(stay_silen.fea)
	hm<-read.table(fea_order)
	
	colid<-c("gene")#The first column is gene IDs
	
	for (i in 1:nrow(hm)){
		colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
	}
	
	colnames(SS)<-colid
	colnames(induced)<-colid
	
	
############################	
#pdf("log2.boxplots.pdf")
#	op<-par(mfrow=c(2,4))
	col_n<-(1+2*(nrow(hm)))# As each HM, both promoter and gene body regions are considered.
	
	for (i in 2:col_n){
		
		if (take.log) {
		induced[,i]<-log2(induced[,i]+pseudo)
		SS[,i]<-log2(SS[,i]+pseudo)
					}


#10e-4, as the real data, second min is 10e-2
		###The list structure!!!
#		boxplot(list(Poised=log2(induced[,i]+pseudo), Silenced=log2(SS[,i]+pseudo)), 
#				col=c("purple", "lightblue"), cex.axis=1.4, las=3, 
#				main=colid[i])
	}
	#	dev.off()
##############################
	
	return (list(Poised=induced, Silenced=SS))
}

	


####################################################################
#Dec.31, 2011 discritize the third variable by k-means
###get kmeans cut points,  where k =4 
Get_Kcut<-function(HMdnt.vector, k=4)
{	
    breakpoints<-c()
    brk.mat<-c()
    for (N in 1:100){
	    fit <- kmeans(HMdnt.vector, k)
		aggregate(HMdnt.vector,by=list(fit$cluster),FUN=mean)
		cell.kmeans <- data.frame(HMdnt.vector, fit$cluster)
	    colnames(cell.kmeans)<-c("histM", "Cst")
#table(cell.kmeans$cluster)
#Figure out the order of "1 2 3 4"as assigned by kmeans
        for (i in 1:k){
	        cluster<-cell.kmeans[cell.kmeans$Cst==i, ]
			breakpoints<-c(breakpoints, min(cluster[,1]), max(cluster[,1]))
		}
		brk.s<-sort(breakpoints)
	    breakpoints<-c()
		brk.mat<-rbind(brk.mat,brk.s)
	}
	
    cut<-c()
	for (i in 1:(2*k)){
#   	k=4, breakpoints=8
	    cut<-c(cut, median(brk.mat[,i]))
	}
	return (cut)
}


###functions need to be checked carefully
###O, L, M, H, from the lowest to highest
###1,2,3,4
Get_level<-function(x, vector){#vector is the cutting points
	x.lev<-c()
	for (i in 1:length(x)){
		if (x[i]>=vector[1] && x[i]<=vector[2]) {
			x.lev[i]<-1 
		} else {
			if (x[i]>=vector[3] && x[i]<=vector[4]){
				x.lev[i]<-2 
			}  else {
				if(x[i]>=vector[5] && x[i]<=vector[6]){
					x.lev[i]<-3 
				} else {
					if(x[i]>=vector[7] && x[i]<=vector[8]){x.lev[i]<-4 }
				}
			}
		}    
		
	}
	return (x.lev)
	
}


##############################################################################			
setwd("/Users/tianr/Desktop/Jan5.2012.Redcell/HSC.model.dataset")

dataset<-Name_log(induced.fea="Redcell.Posi", stay_silen.fea="Redcell.Nega", fea_order="feature.order.redcell.copy", pseudo=0.0001, take.log=FALSE)

setwd("/Users/tianr/Desktop/Jan5.2012.Redcell/Y2012_Jan05_DataCleaning")
			

###For each HM, only keep the more informative (in terms of discrinative ability) one

posi<-subset(dataset$Poised, select=c(gene, H2AZ_Pro,H3K27me1_Body,H3K27me3_Pro,
					H3K36me3_Body, H3K4me1_Body, H3K4me3_Pro,  
					H3K9me1_Body, H3K9me3_Pro,H4K20me1_Body))

nega<-subset(dataset$Silenced, select=c(gene, H2AZ_Pro,H3K27me1_Body,H3K27me3_Pro,
					H3K36me3_Body, H3K4me1_Body, H3K4me3_Pro,  
					H3K9me1_Body, H3K9me3_Pro,H4K20me1_Body))




###Add the third dimension to a x-y plot
##discretize the H3K36me3_Body signals
##poised and silenced genes should be discretized together. Jan5, 2011

all.gene.H3K36me3_Body<-c(posi$H3K36me3_Body,nega$H3K36me3_Body)
K36.Kcut<-Get_Kcut(all.gene.H3K36me3_Body,k=4)

posi.H3K36me3_Body.lev<-Get_level(posi$H3K36me3_Body, K36.Kcut)
nega.H3K36me3_Body.lev<-Get_level(nega$H3K36me3_Body, K36.Kcut)

####
posi.disc<-cbind(posi,posi.H3K36me3_Body.lev)
posi.cleaned<-posi.disc[posi.H3K36me3_Body.lev<=2,]
rm(posi.disc)
rm(posi)

nega.disc<-cbind(nega,nega.H3K36me3_Body.lev)
nega.cleaned<-nega.disc[nega.H3K36me3_Body.lev<=2,]
rm(nega.disc)
rm(nega)

####!!!!
posi<-posi.cleaned
nega<-nega.cleaned


####take care!!! Here the orignal order has changed in terms of its correspondence to genes!!!! Jan5, 2012
write.table(subset(nega, select=c(-nega.H3K36me3_Body.lev)), file="Nega.HSC.cleaned", sep="\t")
write.table(subset(posi, select=c(-posi.H3K36me3_Body.lev)), file="Posi.HSC.cleaned", sep="\t")





#######################################################################################################

plot(nega$H3K4me3_Pro ~ nega$H3K27me3_Pro, pch=16,col=rgb(0,0,200,50,maxColorValue=255), ylim=c(0,50),
	cex=nega$nega.H3K36me3_Body.lev/2.5, 
xlab="H3K27me3 signal intensities", ylab="H3K4me3signal intensities",####!!!!take care the x,y labels!!!
)#blue

points(posi$H3K4me3_Pro ~ posi$H3K27me3_Pro, pch=16,col=rgb(255,0,0,50,maxColorValue=255), 
cex=posi$posi.H3K36me3_Body.lev/2.5) #red

lines(lowess(nega$H3K4me3_Pro ~ nega$H3K27me3_Pro), col="blue", lwd=2)
lines(lowess(posi$H3K4me3_Pro ~ posi$H3K27me3_Pro), col="red", lwd=2)

#legend(50,40,
#c("LL", "L", "M", "H"),
#c("H3K36me3_LL","H3K36me3_L", "H3K36me3_M", "H3K36me3_H"), 
#cex=c(1:4)/2, col="black", pch=16)


leg<-c("H3K36me3_LL","H3K36me3_L", "H3K36me3_M", "H3K36me3_H")

for (i in 1:2){
	chang<-48-i*2
	points(40,chang,pch=16,cex=i/2.5)
	text((40+10), chang, leg[i])
}























