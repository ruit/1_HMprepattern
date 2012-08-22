#Tian R., Dec. 26, 2011
#F score is a classifier independant measure for feature importance.
#Yi-Wei Chen and Chih-Jen Lin, Combining SVMs with Various Feature Selection Strategeies.
#
#
Fscore<-function (posi_data, nega_data){
#input are two data.frames or matrix
#!!!the two datasets should have the same number of columns and the orders should also the same
	positive_df<-as.data.frame(posi_data)
	negative_df<-as.data.frame(nega_data)
	fea_N<-ncol(positive_df)
	fscore<-c()
	for (i in 1:fea_N){
		total_mean<-mean(rbind(positive_df[,i], negative_df[,i]))
		posi_mean<-mean(positive_df[,i])
		nega_mean<-mean(negative_df[,i])
		numerator<-((posi_mean-total_mean)^2+(nega_mean-total_mean)^2)
		demominator<-(var(positive_df[,i])+var(negative_df[,i]))
		fscore<-c(fscore, round(numerator/demominator, 3))
		}

		return (fscore)

	}
