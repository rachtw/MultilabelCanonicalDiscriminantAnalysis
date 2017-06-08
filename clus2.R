args<-commandArgs()
print(args)
dist_no<-args[8]
print(length(args))
if (length(args)==9) {
  suffix<<-args[9]
  suffix2<<-""
} else if (length(args)==10) {
  suffix<<-args[9]
  suffix2<<-args[10]
} else {
  suffix<<-""
  suffix2<<-""
}
size<-4292
species_name<-"yeast"
#dist_no<-35
#suffix<-""
#suffix2<-""
#dist_no<<-0
#suffix<<-".0.6"
#suffix2<<-"t"
if (dist_no==35) {
	size2<-size+1
	dist.data<-read.table(paste(species_name,".dist.",dist_no,suffix,sep=""), sep=",", header=T)[2:size2]
} else if (dist_no==40) {
   dist.data<-read.table(paste(species_name,".dist.",dist_no,suffix,sep=""), header=T)
} else {
	dist.data<-read.table(paste(species_name,".dist.",dist_no,suffix,sep=""), sep=",", header=F)
}
distance=as.dist(dist.data)

# library(cluster)

# ## Hierarchical classification
# clusflex <- agnes(distance,method='ward')
# #clus <- as.hclust(clusflex)
# saveRDS(clusflex, file = paste(species_name,".dist.",dist_no,suffix,".Rda",sep=""))

# nums <- c(114, 400, 715, 4644, 2121)
# g <- list()
# for (i in 1:5) {
	# g[[i]]=cutree(clusflex, k=nums[i])
	# write.csv(g[[i]],paste("clus2",suffix2,".",nums[i],".csv",sep=""))
# }

# postscript(paste("clus2",suffix2,"-2.eps",sep=""))
# plot(clusflex, which.plots = 2, cex= 0.6)
# dev.off()
# postscript(paste("clus2",suffix2,"-1.eps",sep=""))
# plot(clusflex, which.plots = 1)
# dev.off()


## Principal Coordinate Analysis
pco <- cmdscale(distance, k=100, eig = TRUE, x.ret = TRUE)
write.csv(pco$points,paste(species_name,".dist.",dist_no,suffix,".pco.points",sep=""));
write.csv(pco$eig,paste(species_name,".dist.",dist_no,suffix,".pco.eig",sep=""));
write.csv(pco$GOF,paste(species_name,".dist.",dist_no,suffix,".pco.gof",sep=""));

# points<-read.csv(paste(species_name,".dist.",dist_no,suffix,".pco.points",sep=""))

# for (i in 1:5) {
	# postscript(paste(species_name,".dist.",dist_no,suffix,".pco.",nums[i],".eps",sep=""))
	# plot(points[, 2:3], asp = 1, col = g[[i]])
	# dev.off()
# }

# library("stats")
# for (i in 1:5) {
  # km<-kmeans(points[,2:101],nums[i])
  # write.csv(km$cluster,paste("clusk",suffix2,".",nums[i],".csv",sep=""))
# }
