
# data processing
data = rnaseq_processing(data,eratio = 0.06)

rnaseq_processing<-function(data,eratio=0.06,log2=T,median=T,expressionThreshold=2){
  #data: genes in rows and samples in column
  
  # rm cnst genes
  kvar = which(apply(data,1,sd)>0)
  data=data[kvar,]
  
  # control for expression ratio
  if (eratio>0){
    minsample = floor(ncol(data)*eratio)
    kgood = rowSums(data > expressionThreshold) > minsample
    data = data[kgood,]
  }
  
  # log2-transform
  if (log2){data=log2(data+1)}
  
  # median-center
  if (median){data=sweep(data,1,apply(data,1,median))}
  
  return(data)
}

# Tsne plot
tsneRes = Rtsne(t(data[genes,]))
df = as.data.frame(tsneRes$Y)
df = cbind(df,label)
colnames(df) = c("tsne_1","tsne_2","label")
pdf("*.pdf")
g = ggplot(df,mapping = aes(tsne_1,tsne_2,colour=label))+
  geom_point()
# +theme(legend.position='none')
plot(g)
dev.off()

# Gene expression heatmap plot
palette=c("white","blue")
thisPal = brewer.pal(10,"Set3")
cbar = thisPal[label]

cols<-brewer.pal(3, "PuBu")
pal<-colorRampPalette(cols)
my_palette<-pal(10)

pdf("*.pdf")
heatmap.2(dat[genes,],Colv=NA, Rowv=TRUE, 
          symm=FALSE,col=my_palette,breaks = 11,
          density.info="none",trace ="none",scale="none",
          ColSideColors=cbar,dendrogram="row",key = FALSE,
          labRow = NA,
          labCol = NA)
dev.off()


