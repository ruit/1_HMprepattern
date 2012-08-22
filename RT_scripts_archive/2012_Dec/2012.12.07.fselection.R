# Tian R. 
# Dec. 7, 2011
# HSC (CD133+), 9-11 days ex vivo induction into erythroid progenitor cells (CD36+).


###----------------------------------------------------------------------------------

Prep_dataset<-function(induced.fea="induced.list.CD4T.fea", stay_silen.fea="SS.list.CD4T.fea", fea_order="feature_CD4T_order2", pseudo=0.0001){

	induced<-read.table(induced.fea)
	SS<-read.table(stay_silen.fea)
	hm<-read.table(fea_order)
	
	colid<-c("gene")
	
	for (i in 1:nrow(hm)){
		colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
		}
	rm (i)
	colnames(SS)<-colid
#print (head(SS))
	colnames(induced)<-colid
	
	col_n<-(1+2*(nrow(hm)))# As each HM, both promoter and gene body regions are considered.
	
	for (k in 2:col_n){
#print (k)
		
#print (head(SS))
		
		SS[,k]<-log2(SS[,k]+pseudo)
		
#print (head(SS))
		
		induced[,k]<-log2(induced[,k]+pseudo)
		}
	rm (k)
	
	return (list(Posi=induced, Nega=SS))

#pdf("log2.boxplots.pdf")
	
#for (i in 2:col_n){
		
#pseudo<-c(0.0001)#10e-4, as the real data, second min is 10e-2
		
		###The list structure!!!
#		boxplot(list(Stay_Silenced=log2(SS[,i]+pseudo), Induced=log2(induced[,i])+pseudo), 
#			col=c("blue", "lightblue"),
#			main=colid[i])
#		}

#dev.off()
	}

