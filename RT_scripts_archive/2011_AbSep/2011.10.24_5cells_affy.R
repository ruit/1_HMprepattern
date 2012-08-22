#Tian R.@Compbio.tongji, Zhang Lab
#Oct. 24, 2011
#
library("affy")
raw<-dir(pattern=".CEL") #list all the .CEL files
data<-ReadAffy(filenames=raw)
library()

#check affy platform and assign correct CDF file
data@cdfName
#[1] "HG-U133_Plus_2"
data@cdfName<-"hgu133plus2hsrefseqcdf" 
data@cdfName
#[1] "hgu133plus2hsrefseqcdf"

#Perform P/M/A calls
pa<-mas5calls(data)
#Getting probe level data...
#Computing p-values
#Making P/M/A Calls

write.table(exprs(pa), file="Five_cells_PA.expr", sep="\t")
 
######
Oct. 25, 2011
### for duplicates, call PP(present) as 1, AA(absent) as 0,others NA (all discrepancies, or including any M) 
One_exp<-function(mat){
#mat is a n*2 matrix
	expression <-c()
	for (i in 1:dim(mat)[1]){
		if (mat[i,1]==mat[i,2] & mat[i,1]=="P"){
			expression<-c(expression, 1) # P is 1, expressed
			} else {
				if (mat[i,1]==mat[i,2] & mat[i,1]=="A"){
					expression<-c(expression, 0) # A is 0, silenced
					} else {
						expression<-c(expression,NA) # all others are NA
						}
				}
		}
	return (expression)
}


####make 1 + 10 cols into 1+5 cols
dim(data)
pa<-c()
test<-c()
for (i in seq(2,11, by=2)){
	test<-cbind(as.vector(data[,i]), as.vector(data[,i+1]))
	pa<-cbind(pa,One_exp(test))
}
colnames(pa)=c("hES_H1", "HSC", "ERY", "CD4T", "Act_CD4T")
data_puri<-data.frame(RefSeq=data[,1], pa)


