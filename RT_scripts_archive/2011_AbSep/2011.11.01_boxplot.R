# Tian R. Select four sets of genes
# Nov. 1, 2011
# process, CD4+T cell activation 18 hrs.
# constitutively silenced genes before and after antigen activation
SS<-pa[pa$CD4T==0 & pa$Act_CD4T==0,]

#consitutively expressed
EE<-pa[pa$CD4T==1 & pa$Act_CD4T==1,]

# activated only after activation
indu<-pa[pa$CD4T==0 & pa$Act_CD4T==1,]

# repressed only after activation
silen<-pa[pa$CD4T==1 & pa$Act_CD4T==0,]

write.table(SS, file="SS", sep="\t")
write.table(EE, file="EE", sep="\t")
write.table(indu, file="induced", sep="\t")
write.table(silen, file="silenced", sep="\t")


###--------------------------------------------------
EE<-read.table("EE.list.CD4T.fea")
SS<-read.table("SS.list.CD4T.fea")
induced<-read.table("induced.list.CD4T.fea")
silenced<-read.table("silenced.list.CD4T.fea")

hm<-read.table("feature_CD4T_order2")
colid<-c("gene")
for (i in 1:36){
	 colid<-c(colid, paste(hm[i,1], "_Pro", sep=""), paste(hm[i,1], "_Body", sep=""))
	 }

colnames(EE)<-colid
colnames(SS)<-colid
colnames(induced)<-colid
colnames(silenced)<-colid

pdf("72.log2.boxplots.pdf")

for (i in 2:73){
	psudo<-c(0.0001)#10e-4, as the real data, second min is 10e-2
	
	boxplot(list(Stay_Expressed=log2(EE[,i]+psudo), Silenced=log2(silenced[,i]+psudo), Stay_Silenced=log2(SS[,i]+psudo), Induced=log2(induced[,i])+psudo), 
			col=c("blue", "lightblue", "red", "pink"),
			main=colid[i])
}

dev.off()

