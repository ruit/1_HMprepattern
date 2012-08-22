#Tian R., Nov 23, 2011
#Delete 'histone mark feature' one at each time according to importance of each feature (Delete less important features first)
source("CD4T_induction_SVM_RTNov18_2011.R")
# permute all the negative samples by row

posi_n<-nrow(induced)
nega_n<-nrow(SS)

### This is stochastic
index<-sample(seq(1,nega_n,by=1), size=posi_n)

SS_permu<-as.data.frame(SS[index,])

rownames(SS_permu)<-1:length(index)

# the data
data_ori<-Kaotish(Prepare_mat(induced,SS_permu, pseudo=0.0001))

#performance report
perf<-FiveFold(data_ori)

#this order is different from inherent data_ori feature order
feas<-read.table("/Users/tianr/projects/Prepattern_Histone/Blood_Feature_orders_Oct2011/Exp/CD4T_Activation/Induction_Model/Delete_one_feature_each_time/features.remove.order_N23.2011")

data_copy<-data_ori
colnames(feas)<-c("histone")

for (i in 1:nrow(feas)){
	if (ncol(data_copy) >= 3) {	
		col_n<-which(colnames(data_copy)==as.character(feas$histone[i]))
		data_partial<-data_copy[,-col_n]
		print (ncol(data_partial))
		data_copy<-data_partial###update the data!!!
		perf<-rbind(perf,FiveFold(data_partial))
		}
	}
write.table(perf, file="delete.fea.report")
