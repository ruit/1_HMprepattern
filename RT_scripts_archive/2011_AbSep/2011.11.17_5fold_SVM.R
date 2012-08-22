# Tian R. Select four sets of genes
# Nov. 17, 2011
# process, CD4+T cell activation 18 hrs.
# constitutively silenced genes before and after antigen activation
#SS<-pa[pa$CD4T==0 & pa$Act_CD4T==0,]

#consitutively expressed
#EE<-pa[pa$CD4T==1 & pa$Act_CD4T==1,]

# activated only after activation
#indu<-pa[pa$CD4T==0 & pa$Act_CD4T==1,]

#repressed only after activation
#silen<-pa[pa$CD4T==1 & pa$Act_CD4T==0,]

#write.table(SS, file="SS", sep="\t")
#write.table(EE, file="EE", sep="\t")
#write.table(indu, file="induced", sep="\t")
#write.table(silen, file="silenced", sep="\t")


###--------------------------------------------------
#EE<-read.table("EE.list.CD4T.fea")
SS<-read.table("SS.list.CD4T.fea")
induced<-read.table("induced.list.CD4T.fea")
#silenced<-read.table("silenced.list.CD4T.fea")

hm<-read.table("feature_CD4T_order2")
colid<-c("gene")
for (i in 1:36){
	 colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
	 }

#colnames(EE)<-colid
colnames(SS)<-colid
colnames(induced)<-colid
#colnames(silenced)<-colid

#pdf("72.log2.boxplots.pdf")

#for (i in 2:73){
#	pseudo<-c(0.0001)#10e-4, as the real data, second min is 10e-2
	
#	boxplot(list(Stay_Expressed=log2(EE[,i]+pseudo), Silenced=log2(silenced[,i]+pseudo), Stay_Silenced=log2(SS[,i]+pseudo), Induced=log2(induced[,i])+pseudo), 
#			col=c("blue", "lightblue", "red", "pink"),
#			main=colid[i])
#}

#dev.off()






#################TianR., Nov.10, 2011
# permutation with sample()
Kaotish<-function(mat){
#mat_new : order is randomly changed
#mat is a matrix or data.frame
#row_ord is the index of permutated matrix
	mat_new<-c()
	row_ord<-sample(nrow(mat))#replacement=False
	
	for (i in 1:nrow(mat)){
		mat_new<-rbind(mat_new, mat[row_ord[i],])
		
		}
	return (mat_new)
	}



Prepare_mat<-function(induced,SS, pseudo=0.0001){
		#induced is the positive sample data.frame, 
		#SS is the negative sample data.frame
		#take all the 381 positive samples and for the moment take top 400 negative samples
	posi_sample_n<-dim(induced)[1]
	dataTcell<-rbind(subset(induced, select=c(-gene)), subset(SS[1:posi_sample_n,], select=c(-gene)))
###log tranformation or other kind of scaling method is import for the performance!
### the pseudo number is critical!
#pseudo<-c(0.0001)#10e-4
	dataTcell<-log2(dataTcell+pseudo)
	classesTcell<-c(rep("1",posi_sample_n), rep("0",posi_sample_n))
	data<-data.frame(dataTcell, class=classesTcell)
	return (data)
}
#table(data$class)

#0   1   ### positive samples are labeled with "1" 
#400 381 

#data_rand<-Kaotish(Kaotish(data))
#dim(data_rand)

#table(data_rand[1:100,]$class)

#0  1 
#54 46 
#data_rand_f<-subset(data_rand, select=c(-class))
#dim(data_rand_f)
#[1] 781  72
#classes_rand<-subset(data_rand, select=class)
#dim(classes_rand)
#[1] 781   1

#data_tr<-data_rand_f[1:381,]
#classes_tr<-classes_rand[1:381,]
#data_test<-data_rand_f[382:781,]
#classes_test<-classes_rand[382:781,]

###
#performance evaluation
#library("e1071")
#载入需要的程辑包：class
#model<-svm(data_tr, classes_tr)
#pred<-predict(model, data_test)

#compare<-data.frame(real=classes_test, pred=pred )
#table(compare$real)

#0   1 
#202 198 

####!!!!!
#Nega<-202
#Posi<-198

#TP<-dim(compare[compare$real==1 & compare$pred==1,])
#TN<-dim(compare[compare$real==0 & compare$pred==0,])

#Sensi_TPR<-100*TP[1]/Posi[1]
#Speci_TNR<-100*TN[1]/Nega[1]
#Accuracy<-100*(TP[1]+TN[1])/(Nega+Posi)


##########################
##########################
##########################
#Tian R., Nov. 11, 2011
FiveFold<-function(data){
#perform five fold cross validation for SVM
#data, is all the features and the lables (the last column titled 'class')in a data.frame
#data should be properly scaled for the range of value and also balanced for positive and negative samples.
#this function depends on package 'e1071' and function Kaotish
	library("e1071")
	#permutate the data
	data_rand<-Kaotish(Kaotish(data))
	TP<-c()# True positive, '1' is predicted as '1'
	TN<-c()# True negatvie, '0' is predicted as '0'
	compare<-c()	
	data_rand_f<-subset(data_rand, select=c(-class))
	
	classes_rand<-subset(data_rand, select=class)
	
	n_sample<-nrow(data)
	#print (n_sample)
	group_size<-(n_sample) %/% 5 #### take care %/% the divison by integers
	#print (group_size)
	row_ord_permu<-sample(seq(1:n_sample))# again permutate the order of the samples
	start<-0
	end<-0
	for (i in 1:5){
		#update the training and testing data, make sure each group of samples tested once and only once
		start<-((i-1)*group_size+1)
		if (i<5){ 
			end<-(end+group_size)} else {
							end<-n_sample}              
		index<-row_ord_permu[start:end]

		#take care the last "end"!! you add up to come with a "uplimit", that can be bigger or smaller than the real one
		#make the training data as the rest of the testing data, each one group
		data_tr<-data_rand_f[-index,]
		classes_tr<-classes_rand[-index,]
		#make the testing data
		data_test<-data_rand_f[index,]
		classes_test<-classes_rand[index,]
		
		#build model and make prediction
		model<-svm(data_tr, classes_tr)
	        pred<-predict(model, data_test)
		
		#Assess the performance of the SVM
		compare<-data.frame(real=classes_test, pred=pred )
		#print (head(compare))						
		TP<-c(TP,as.vector(dim(compare[compare$real==1 & compare$pred==1,])[1]))
		TN<-c(TN,as.vector(dim(compare[compare$real==0 & compare$pred==0,])[1]))
		compare<-c()
		#take care the last "end"!! you add up to come with a "uplimit", that can be bigger or smaller than the real one

		}
	
	Posi<-table(data$class)[["1"]]
	Nega<-table(data$class)[["0"]]
	Sensi_TPR<-100*sum(TP)/Posi
	Speci_TNR<-100*sum(TN)/Nega
	Accuracy<-100*(sum(TP)+sum(TN))/(Nega+Posi)
	
	return (data.frame(Sensitivity=Sensi_TPR, Specificity=Speci_TNR, Accuracy=Accuracy))
		
	}


