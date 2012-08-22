#files<-dir(pattern="Pro")
#files<-dir(pattern="Body")

files<-dir(pattern="pred")

sp<-c()
for (i in 1:72){
	pred<-read.table(files[i])
	sp<-cbind(sp, pred$Specificity)
}

boxplot(sp, col=c(1,2))

sen<-c()
for (i in 1:72){
	pred<-read.table(files[i])
	sen<-cbind(sen, pred$Sensitivity)
}
boxplot(sp, col=c("blue", "lightblue"))




##Tian R., Nov 23, 2011

pdf("/Users/tianr/Desktop/permu.pdf")

sen<-c()
for (i in 1:72){
	pred<-read.table(files[i])
	sen<-c(sen, min(pred$Sensitivity))
}

plot(sen, col=c("red", "black"), pch=21, type="b", lty=2, ylab="Predictive performance (%)",
main="Permutation importance of each feature", xlab="Order of the features")



sp<-c()
for (i in 1:72){
	pred<-read.table(files[i])
	sp<-c(sp, min(pred$Specificity))
}
lines(sp, col=c("red", "black"), pch=20, type="b")

legend(1,55, c("Sensitivity", "Specificity"), pch=c(21,20), lty=c(2,1), col="red")


dev.off()




######Tian R., Nov 24, 2001

hm<-read.table("/Users/tianr/projects/Prepattern_Histone/Blood_Feature_orders_Oct2011/Exp/CD4T_Activation/Induction_Model/feature_CD4T_order2")
#colid<-c("gene")
colid<-c()
for (i in 1:36){
	colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
}


pred<-c()
perf<-c()
for (i in 1:length(colid)){
	file<-paste("pred_power_", colid[i], sep="")
	pred<-read.table(file)
	perf<-c(perf, min(pred$Specificity))
}
data.frame(feature=colid, importance=perf)


# delete one feature each step
pdf("/Users/tianr/Desktop/fea.del.pdf")
plot(perf$Sensitivity, type="b", col="blue", ylim=c(50,82), pch=21, main="Delete one feature sequentially", ylab="Predictive performance", xlab="Number of features deleted")
lines(perf$Specificity, type="b", col="red", pch=20)
legend(0, 60, c("Sensitivity", "Specificity"), col=c("blue", "red"), pch=c(21,20))
dev.off()




