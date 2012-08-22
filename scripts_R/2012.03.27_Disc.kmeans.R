#Tian R., Zhang Y.
#March 27, 2012
#Disretize HM densities data by k-means
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

############################
Get_level<-function(x, vector){#vector is the cutting points
	x.lev<-c()
	for (i in 1:length(x)){
		if (x[i]>=vector[1] && x[i]<=vector[2]) {
			x.lev[i]<-"Lev0"
		} else {
			if (x[i]>=vector[3] && x[i]<=vector[4]){
				x.lev[i]<-"Lev1" 
			}  else {
				if(x[i]>=vector[5] && x[i]<=vector[6]){
					x.lev[i]<-"Lev2" 
				} else {
					if(x[i]>=vector[7] && x[i]<=vector[8]){x.lev[i]<-"Lev3" }
				}
			}
		}    
		
	}
	return (x.lev)
	
}

########################
Disc.main<-function(hm.file){
	
	#first col is gene ids
	vec.disc<-c()
	total<-read.table(hm.file, header=T)
	k.nol<-ncol(total)
	for (m in 2:k.nol){
		vec.disc<-cbind(vec.disc,Get_level(total[,m],Get_Kcut(total[,m])))
		}
	colnames(vec.disc)<-colnames(total[,2:k.nol])

	gm.disc<-data.frame(gene=total[,1], vec.disc)
	
	return (gm.disc)
}
