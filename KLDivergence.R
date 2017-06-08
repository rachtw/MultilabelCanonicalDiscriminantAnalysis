args<-commandArgs()
print(args)
#args <- c()
#args[8] <- "yeast.fasta.sgml.freq"
#args[9] <- "yeast.dist.1"
freq<-read.table(args[8],sep=",")
nproteins<-dim(freq)[1]
ncolumns<-dim(freq)[2]
start<-args[10]
end<-args[11]
if (end > nproteins) { end<-nproteins }

library("entropy")
#distance<-array(0,dim=c(nproteins,nproteins))
if (start==1) start=2
for (i in start:end) {
  print(paste("row",i))
  for (j in 1:(i-1)) {
    if (j%%100==0) { print(paste("column",j)) }
    index<-intersect(which(freq[i,]!=0),which(freq[j,]!=0))
    if (length(index)==0) {
    #  distance[i,j] <- ncolumns
      print(paste("d",i,j,ncolumns))
    } else {
    #  distance[i,j] <- (KL.plugin(freq[i,index], freq[j,index])+KL.plugin(freq[j,index], freq[i,index]))/2
    #  distance[j,i] <- distance[i,j]
      print(paste("d",i,j,(KL.plugin(freq[i,index], freq[j,index])+KL.plugin(freq[j,index], freq[i,index]))/2))
    }
  }
}
#write.table(distance,args[9],row.names=FALSE, col.names=FALSE, sep=",")
q()
