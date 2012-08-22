#Tian R. & Zhang Y.
#Jan 06, 2012
#Perform PCA, take care scaling of the original variables and keep the order of genes

########################################################################
setwd("/Users/tianr/Desktop/Jan5.2012.Redcell/Y2012_Jan6_PCA")
dir()
#[1] "Nega.HSC.cleaned_Jan5.TAB" "Posi.HSC.cleaned_Jan5.TAB"
nega<-read.table("Nega.HSC.cleaned_Jan5.TAB", header=T)
posi<-read.table("Posi.HSC.cleaned_Jan5.TAB", header=T)

dim(nega)
#[1] 7751   10
dim(posi)
#[1] 356  10

head(posi, 2)
#gene H2AZ_Pro H3K27me1_Body H3K27me3_Pro H3K36me3_Body H3K4me1_Body H3K4me3_Pro
#1 NM_000032     2.75      3.037293         5.25      2.768166     2.537486        0.50
#2 NM_000037     9.50      1.893569        29.00      2.291151     3.396295        6.75
#H3K9me1_Body H3K9me3_Pro H4K20me1_Body
#1     5.344098        1.75      3.844675
#2     3.274999        2.50      3.409772

head(nega, 2)
#gene H2AZ_Pro H3K27me1_Body H3K27me3_Pro H3K36me3_Body H3K4me1_Body H3K4me3_Pro
#1 NM_000014     4.75      5.721940         0.50      7.807865    5.2243804        0.75
#2 NM_000015     2.75      1.145393         4.75      2.433961    0.5011096        0.50
#H3K9me1_Body H3K9me3_Pro H4K20me1_Body
#1     3.559468        1.75      3.272414
#2     0.859045        6.75      1.360155


all<-rbind(posi, nega)

all.9HM<-subset(all, select=c(-gene))

#log transform and scale the variables
all.9HM.log.centered<-scale(log2(all.9HM+0.0001))

#################################################
#plot to see the effect of scaling

op<-par(mfrow=c(2,1))
ori<-all.9HM
colnames(ori)<-NULL
boxplot(ori, main="Original HM ChIP tag Densities")
rm (ori)
#ann=F, 
#axes=F
#box()
new<-all.9HM.log.centered
colnames(new)<-c()
boxplot(new, main="After Scaling")
rm(new)
#########################################################

hm.pca<-princomp(all.9HM.log.centered, cor=T)

summary(hm.pca)

#Importance of components:
#Comp.1    Comp.2    Comp.3    Comp.4    Comp.5    Comp.6    Comp.7
#Standard deviation     1.7325353 1.3424711 1.1551573 0.9746730 0.7922755 0.7064805 0.6227540
#Proportion of Variance 0.3335199 0.2002476 0.1482654 0.1055542 0.0697445 0.0554572 0.0430914
#Cumulative Proportion  0.3335199 0.5337675 0.6820329 0.7875870 0.8573315 0.9127887 0.9558801
#Comp.8    Comp.9
#Standard deviation     0.52229628 0.3525414
#Proportion of Variance 0.03031038 0.0138095
#Cumulative Proportion  0.98619050 1.0000000


loadings(hm.pca)

#Loadings:
#Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7 Comp.8 Comp.9
#H2AZ_Pro      -0.254 -0.250  0.431 -0.499 -0.383  0.340  0.319  0.273       
#H3K27me1_Body -0.327  0.404  0.325 -0.134        -0.340 -0.543  0.437       
#H3K27me3_Pro   0.235 -0.556  0.133                0.332 -0.709              
#H3K36me3_Body -0.248  0.181  0.542  0.317  0.498  0.425        -0.274       
#H3K4me1_Body  -0.498                      -0.414 -0.107 -0.189 -0.624  0.376
#H3K4me3_Pro   -0.106 -0.546  0.209 -0.245  0.461 -0.578  0.124 -0.155       
#H3K9me1_Body  -0.507 -0.157 -0.215  0.205                             -0.788
#H3K9me3_Pro    0.208 -0.123  0.470  0.612 -0.434 -0.348  0.162              
#H4K20me1_Body -0.389 -0.306 -0.292  0.389  0.151  0.105  0.115  0.487  0.483


###############################################
weights<-loadings(hm.pca)
#weights[,1:4]
weights.copy<-weights
#colnames(weights)<-c("")

par(mfrow=c(2,2), las=3)
barplot(weights[,1], col=c(7,7,6,7,7,7,7,6,7),main="1st PC")
barplot(weights[,2], col=c(7,6,7,6,7,7,7,7,7),main="2nd PC")
barplot(weights[,3], col=c(6,6,6,6,7,6,7,6,7),main="3rd PC")
barplot(weights[,4], col=c(7,7,6,6,6,7,6,6,6),main="4th PC")

##############################################
#biplot(hm.pca)

#head(hm.pca$scores)
#Comp.1      Comp.2     Comp.3     Comp.4     Comp.5      Comp.6
#1 -0.484361109  0.45126284 -0.2303485  0.5663207  0.3421223 -0.10296281
#2 -0.003978105 -1.21211643  0.5124933 -0.1315797 -0.1558292  0.14424572

#Comp.7      Comp.8      Comp.9
#1 -0.13208846  0.08223472 -0.50173523
#2 -0.36246430 -0.23849146 -0.03006746

PC.fea<-hm.pca$scores[,1:4]

head(PC.fea)

#Comp.1      Comp.2     Comp.3     Comp.4
#1 -0.484361109  0.45126284 -0.2303485  0.5663207

all.13fea<-cbind(all, PC.fea)

head(all.13fea, 2)
#gene H2AZ_Pro H3K27me1_Body H3K27me3_Pro H3K36me3_Body H3K4me1_Body
#1 NM_000032     2.75      3.037293         5.25      2.768166     2.537486
#2 NM_000037     9.50      1.893569        29.00      2.291151     3.396295
#H3K4me3_Pro H3K9me1_Body H3K9me3_Pro H4K20me1_Body       Comp.1     Comp.2
#1        0.50     5.344098        1.75      3.844675 -0.484361109  0.4512628
#2        6.75     3.274999        2.50      3.409772 -0.003978105 -1.2121164
#Comp.3     Comp.4
#1 -0.2303485  0.5663207
#2  0.5124933 -0.1315797



