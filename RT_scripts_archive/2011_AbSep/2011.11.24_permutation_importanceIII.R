#Tian R., Nov 24, 2011
###!!!!! update and keep intact !!!!!!
#Permute by column 'histone mark feature' each 100 times to get an idea the importance of each feature
source("CD4T_induction_SVM_RTNov18_2011rev.R")
# permute all the negative samples by row

#posi_n<-nrow(induced)
#nega_n<-nrow(SS)

#index<-sample(seq(1,nega_n,by=1), size=posi_n)

#SS_permu<-as.data.frame(SS[index,])

#rownames(SS_permu)<-1:length(index)

# the data
data_ori<-Kaotish(Prepare_mat(induced,SS, pseudo=0.0001))
# Permute by each column

#################################################################
# use a fix dataset, Tcell_TCR_data, Nov. 24, 2011


fnum<-ncol(data_ori)

fnames<-as.vector(colnames(data_ori))

rownames(data_ori)<-1:nrow(data_ori)

#make an unchanged copy!!!!Nov 24, 2011
data_copy<-data_ori

for (i in 49:72){
	HMpred<-c()
	###make sure permutation starts with an original data
	data_ori<-data_copy
	for (k in 1:100){
		###repeat 100 times
	
		fea_permu<-Kaotish(as.matrix(subset(data_ori, select=fnames[i])))
		colnames(fea_permu)<-fnames[i]
		data_ori[,i]<-fea_permu
		data_var<-data_ori
		###Caution! make sure the data is permutated as desired!
		HMpred<-rbind(HMpred,FiveFold(data_var))		
		
		}
		write.table(HMpred, file=paste('Pred_Power_',fnames[i], sep=''), sep='\t')
		HMpred<-c()	
	}

