
rm(list = ls())
load("Biase.RData")
library(matrixStats)
## log-transformed
Mdata<-max(data)  
if(Mdata>10000){
  data<-log(data+1)/log(10)
}else {
  data<-log(data+1)/log(2)
}

## filter gene expression 
## retained the v genes (default: 1000) with the highest variance 
data<-as.matrix(
  data[rowSums(data)>0,]
)
ncol<-dim(data)[2]
sd<-matrixStats::rowVars(data)
data<-cbind(data,sd)
newdata<-data[order(-sd),]
newdata<-newdata[1:1000,1:ncol] 

## calculate the gene kernel density
den=512
name<-rownames(newdata)
nrow<-dim(newdata)[1]
ncol<-dim(newdata)[2]
d<-list()
for(i in 1:nrow){
  a<-density(as.numeric(newdata[i,]),n=den)
  d[[i]]=a$y
}
d<-as.data.frame(d)
rm(a)
names(d)<-name
d<-as.data.frame(t(d),)
write.table(d,file="Biase_densty.csv",sep = ',',col.names = FALSE)
write.table(newdata,file = "Biase.csv",sep = ',',col.names = TRUE)
