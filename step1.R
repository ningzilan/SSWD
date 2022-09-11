
rm(list = ls())
load("Biase.RData")
Mdata<-max(data)
if(Mdata>10000){
  data<-log(data+1)/log(10)
}else {
  data<-log(data+1)/log(2)
}

library(matrixStats)
data<-as.matrix(
  data[rowSums(data)>0,]
)
sd<-matrixStats::rowVars(data)
data<-cbind(data,sd)
newdata<-data[order(-sd),]
newdata<-newdata[1:1000,1:ncol]


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
