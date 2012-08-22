# Tian R. Select four sets of genes
# Nov. 25, 2011
# Dec. 6,  2011
# run time over 12 min

####-------------------------------------------------------------------------------------
## the two original files are without column names, class labels are not included either.

Name_feature<-function(induced.file="induced.list.CD4T.fea", SS.file="SS.list.CD4T.fea",fea.order="feature.order" ){
	induced<-read.table(induced.file)
	SS<-read.table(SS.file)
    hm<-read.table(fea.order)
    colid<-c("gene")
    for (i in 1:nrow(hm)){
		colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
		}
	colnames(SS)<-colid
	colnames(induced)<-colid
	
	return(list(posi=induced, nega=SS))
	rm (induced)
	rm (SS)
}



###---------------------------------------------------------------------------------------
Prepare_mat<-function(induced="the positive samples",SS="the negative samples", pseudo=0.0001){
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
	rm (induced)
	rm (SS)
	}

###---------------------------------------------------------------------------------------
# permutation with sample()
Kaotish<-function(mat){
#mat_new : order is randomly changed
#mat is a matrix or data.frame
#row_ord is the index of permutated matrix
	mat_new<-c()
	hang_n<-nrow(mat)
## permute by rows, Nov 25, 2011
	row_ord<-sample(seq(1, hang_n, by=1))#replacement=False
	mat_new<-mat[row_ord,]
# looping is slow!! Nov. 25, 2011 revised
	
	return (mat_new)
	rm (mat)
	}

###---------------------------------------------------------------------------------------
FiveFold<-function(data="complete data.frame including labels"){
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
	
	rm(data)
		
	}


##--------------------------------------------------------------------------------------------
# balance by take partial negative sampel for building classifier; Tian R, Dec 6, 2011, deal with imbanced (excessive) negative samples
Summa_nega<-function(induced="the positive", SS="the excessive negative"){
	#
	partial_SS<-c()	
	result_df<-c()
	
	posi_sample_size<-nrow(induced)
	
	SS_permu<-Kaotish(SS)
	
	nega_sample_n<-nrow(SS)
	
	row_ord_permu<-sample(seq(1:nega_sample_n))
	start<-0
	end<-0
	for (i in seq(1,nega_sample_n, by=posi_sample_size)){
		start<-i
		#print (start)
		end<-(i+posi_sample_size-1)	
		
		if (end > nega_sample_n){
			end <- nega_sample_n}
		
		partial_SS<-SS[start:end,]
		
		result_df<-rbind(result_df,FiveFold(Prepare_mat(induced, SS)))	
		rm (partial_SS)
		
		}
	
	return (result_df)
	rm (SS)
	rm (induced)
	}

##-------------------------------------------------------------------------------------------

