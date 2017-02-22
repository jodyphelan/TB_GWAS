library(scales)


args = commandArgs(trailingOnly=TRUE)
pca<-read.table(args[1])
vars <- read.table(args[2])[,1]
meta<-read.table(args[3])

meta<-meta[match(rownames(pca),meta$V1),]

meta.uniq<-unique(meta$V2)
cols.uniq<-alpha(rainbow(length(meta.uniq)),0.5)
cols<-cols.uniq[match(meta$V2,meta.uniq)]

leg_text<-paste(meta.uniq, " (n=", table(meta$V2)[match(meta.uniq,names(table(meta$V2)))], ")", sep="")

xlabs<-paste("PC1 (",vars[1],"%)",sep="")
ylabs<-paste("PC1 (",vars[2],"%)",sep="")
png(args[4],800,600)
plot(pca[,1],pca[,2],col=cols,pch=20,cex=2,xlab=xlabs,ylab=ylabs)
temp<-numeric()
for (x in c("topleft","topright","bottomright","bottomleft")){
	leg<-legend(x,legend=leg_text,fill=cols.uniq,plot=F)
	temp<-c(temp,length(which(pca$points[,1]>leg$rect$left & pca$points[,1]<(leg$rect$left+leg$rect$w) & pca$points[,2]<leg$rect$top & pca$points[,2]>(leg$rect$top-leg$rect$h))))
}
bestplace<-c("topleft","topright","bottomright","bottomleft")[order(temp)[1]]
legend(bestplace,legend=leg_text,fill=cols.uniq)

dev.off()

