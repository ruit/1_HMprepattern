
###Tian R., Yong Zhang, Dec.26, 2011
###Redraw the corrplot for redcell model
####



dir()
#[1] "feature.order.redcell.copy" "Redcell.Nega"              
#[3] "Redcell.Posi"              
nega<-read.table("Redcell.Nega")
posi<-read.table("Redcell.Posi")
dim(nega)
#[1] 8007   19
dim(posi)
#[1] 369  19

nega_mat<-nega[,2:19]
posi_mat<-posi[,2:19]

hm.fea<-c("H2AZ_Pro", "H2AZ_Body", "H3K27me1_Pro", "H3K27me1_Body", "H3K27me3_Pro",
 "H3K27me3_Body", "H3K36me3_Pro", "H3K36me3_Body",
 "H3K4me1_Pro", "H3K4me1_Body", "H3K4me3_Pro",  
 "H3K4me3_Body", "H3K9me1_Pro", "H3K9me1_Body", 
 "H3K9me3_Pro", "H3K9me3_Body", "H4K20me1_Pro", 
 "H4K20me1_Body")
 
colnames(posi_mat)<-hm.fea
colnames(nega_mat)<-hm.fea

total_mat<-rbind(posi_mat, nega_mat)

library("corrplot")###From Taiyun

col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
"#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200) ###!!!!!


corrplot(cor(total_mat), order ="hclust", type="upper", col=rev(col), 
	tl.col="black", title="Correlation relationships across histone mark features", 
	mar=c(0,0,3,4), 
	tl.cex=1.2 )       ## your taste
