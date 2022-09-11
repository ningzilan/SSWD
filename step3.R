
file1<-read.table(file="Biase_subcluster.txt",sep=",",header=T)
gc<-file1$Y_optimal
rm(file1)
newdata=as.data.frame(newdata)
gc=as.data.frame(gc)
names(gc)="y"
gcluster=cbind(newdata,gc)
gcluster<-gcluster[order(gcluster$y),]
gcluster<-t(gcluster)
write.table(gcluster,file="Biase_subcluster.csv",sep = ',')
