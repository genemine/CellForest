
# load R packages
library(SingleCellExperiment)
library(SC3)
library(scater)
library(igraph)

# load functions

# SC3 algorithm
# data A expression data matrix,gene in rows and samples in columns.
# label The cell labels provided by authors of the original publication.
# genes List of genes selected by different gene selection strategies.
# return clustering label
SC3Cluster<-function(data,label,genes){
  if(length(label)==1)
    label=label[,1]
  
  x <- data[genes,]
  lab <- label
  lab = as.matrix.Vector(lab)
  rownames(lab) = colnames(x)
  cty = "cell_type1"
  cty = as.matrix.Vector(cty)
  colnames(cty)=cty[,1]
  colnames(lab) <- colnames(cty)
  
  # create a SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(x),
      # logcounts = log2(as.matrix(x) + 1)
      logcounts = as.matrix(x)
    ), 
    colData = lab
  )
  
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce <- sc3(sce, ks = kprior,gene_filter = FALSE,biology = FALSE)
  return(eval(parse(text = paste0("colData(sce)$sc3_",kprior,"_clusters"))))
}

# Louvain algorithm
# data A expression data matrix,gene in rows and samples in columns.
# genes List of genes selected by different gene selection strategies.
# return clustering label
SNNLouvain<-function(data,genes){
  cat("SNN graphs construction...\n")
  k = 30
  nsamples = ncol(data)
  dist = as.matrix(dist((t(data[genes,]))))
  IDX = t(apply(dist,1,order)[1:k,])
  rownames(IDX) = c(1:nsamples)
  EdgeMat = matrix(0,ncol = 3,nrow = 0)
  for (x in 1:(nsamples-1)){
    for (y in (x+1):nsamples){
      shared = intersect(IDX[x,], IDX[y,])
      if(length(shared) > 0){
        s = k-0.5*(match(shared, IDX[x,])+match(shared, IDX[y,]))
        strength = max(s)
        if (strength > 0)
          EdgeMat = rbind(EdgeMat, c(rownames(IDX)[x],rownames(IDX)[y],strength))
      }
    }
  }
  EdgeMat = as.data.frame(EdgeMat)
  colnames(EdgeMat) = c("N1","N2","w")
  
  cat("SNN graphs combining...\n")
  comgraph = data.frame(n1 = EdgeMat[,1],
                        n2 = EdgeMat[,2],
                        weight = as.numeric(EdgeMat[,3]))
  comgraph = graph.data.frame(comgraph, directed = F)
  cluster = cluster_louvain(comgraph)$membership
  return(cluster)
}

# load single cell RNA seq data
load("*.Rdata")

# k-means clustering algorithm
cluster = kmeans(t(data[genes,]),centers = kprior,iter.max = 1e+09,nstart = 1000)$cluster
NMI = evalcluster(label,cluster)[[1]]
ARI = evalcluster(label,cluster)[[3]]

# Clustering LARge Applications
cluster = clara(t(data[genes,]),k = kprior)$clustering
NMI = evalcluster(label,cluster)[[1]]
ARI = evalcluster(label,cluster)[[3]]

# Hierarchical Clustering
cluster = cutree(hclust(as.dist(1-cor(data[genes,]))),kprior)
NMI = evalcluster(label,cluster)[[1]]
ARI = evalcluster(label,cluster)[[3]]

# SNN louvain
cluster = SNNLouvain(data,genes)
NMI = evalcluster(label,cluster)[[1]]
ARI = evalcluster(label,cluster)[[3]]

# SC3(Average of ten runs)
NMI = c()
ARI = c()
for (iter in 1:10) {
  cluster = SC3Cluster(data,genes)
  NMI = c(NMI,evalcluster(label,cluster)[[1]])
  ARI = c(ARI,evalcluster(label,cluster)[[3]])
}
NMI = mean(NMI)
ARI = mean(ARI)
