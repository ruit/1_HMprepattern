# Tian R. 
# Dec. 5, 2011
# HSC (CD133+), 9-11 days ex vivo induction into erythroid progenitor cells (CD36+).


###----------------------------------------------------------------------------------

Plot_HM<-function(induced.fea="induced.fea", stay_silen.fea="SS.fea", fea_order="feature.order", pseudo=0.0001){
	

	induced<-read.table(induced.fea)
	SS<-read.table(stay_silen.fea)
	hm<-read.table(fea_order)
	
	colid<-c("gene")#The first column is gene IDs
	
	for (i in 1:nrow(hm)){
		colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
		}

	colnames(SS)<-colid
	colnames(induced)<-colid

	

	pdf("log2.boxplots.pdf")
	
op<-par(mfrow=c(2,4))

	col_n<-(1+2*(nrow(hm)))# As each HM, both promoter and gene body regions are considered.
	
	for (i in 2:col_n){
		
#10e-4, as the real data, second min is 10e-2
		
		###The list structure!!!
		boxplot(list(Poised=log2(induced[,i]+pseudo), Silenced=log2(SS[,i]+pseudo)), 
			col=c("purple", "lightblue"), cex.axis=1.4, las=3, 
			main=colid[i])
		}

	dev.off()
	}

