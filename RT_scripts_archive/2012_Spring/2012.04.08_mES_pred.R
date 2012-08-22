###April 8, 2012

plot(sort(log10(silen[sample(1:661, 300),2]+0.00001)), col="blue")
points(sort(log10(poi[sample(1:4346, 300),2]+0.00001)), col="red")
boxplot( list(silen=log2(silen[sample(1:661, 300),2]+0.00001), poi=log2(poi[sample(1:4346, 300),2]+0.00001)), col=c("blue", "red"))


silen<-read.table("silen.mES.eml.expr")
poi<-read.table("poised.mES.eml.expr")
dim(silen)


silen.sum<-c()
poi.sum<-c()

for (i in 1:100){ 
	
	s.ind<-sample(1:661, 300)
	p.ind<-sample(1:4346, 300)
	S<-silen[s.ind, ]
	P<-poi[p.ind, ]
	
	silen.sum<-c(silen.sum, nrow(S[S[,2]>0.001,])/nrow(S))
	
	poi.sum<-c(poi.sum, nrow(P[P[,2]>0.001,])/nrow(P))
	
}

