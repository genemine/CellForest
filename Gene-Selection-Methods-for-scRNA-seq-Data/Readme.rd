 
# Expr
## Description
Expr is a gene selection method based on single-cell RNA sequencing data which only retains the genes with the highest average
expression (log-normalized count) value across all cells. Usually 10% of the original gene number is retained (counts in at least one cell are non-zero).

## Usage
Note: The input data is a gene expression matrix, gene in rows and cells in columns.

Step 1: open your R or Rstudio 
Step 2: in the R command window, run the following command to calculate the average expression value of genes across all cells
```
> meanExpr <- apply(data,1,mean)
```
Step 3: in R command window, run the following command to sort the genes according to the average expression value
```
> meanExpr <- meanExpr[order(meanExpr,decreasing = TRUE)]
```
Step 4: At the end of the help page, run the following command to retain 10% of the original number of genes
```
> Exprgenes <- names(meanExpr)[1:floor(0.1*nrow(data))]
```


# HVG
## Description
The gene selection method in R package Seurat which estimate the variability of the features and retain only the most highly variable ones. Usually 10% of the original gene number is retained (counts in at least one cell are non-zero).
Satija R, Farrell JA, Gennert D, et al.: Spatial reconstruction of single-cell gene expression data. Nat Biotechnol. 2015; 33(5): 495–502.

## Usage
Note: The input data is a gene expression matrix, gene in rows and cells in columns. This method has more parameters and can be adjusted according to the data.

Step 1: open your R or Rstudio 
Step 2: in the R command window, run the following command to install and load the R package
```
> install.packages('Seurat')
> library(Seurat)
```
Step 3: in R command window, run the following command to create a Seurat object from a feature (e.g. gene) expression matrix
```
> Sobj <- CreateSeuratObject(data)
```
Step 4: in R command window, run the following command to identify features that are outliers on a 'mean variability plot'
```
> Sobj <- FindVariableFeatures(object = Sobj,nfeatures = floor(0.1*nrow(data)))
```
Step 5:  At the end of the help page, run the following command to obtain the selected genes
```
> HVGgenes <- Sobj[["RNA"]]@var.features
```


# scMarker
## Description
A new marker selection strategy (SCMarker) to accurately delineate cell types in single cell RNA-sequencing data by identifying genes
that have bi/multi-modally distributed expression levels and are co- or mutually-exclusively expressed with some other genes.
Wang, F. et al. (2019). SCMarker: ab initio marker selection for single cell transcriptome profiling PLoS Comput Biol, 15(10), e1007445

## Usage
Note: The input data is a gene expression matrix, gene in rows and cells in columns. 

Step 1: open your R or Rstudio 
Step 2: in the R command window, run the following command to install and load the R package
```
> library(devtools)
> install_github("KChen-lab/SCMarker")
> library(SCMarker)
```
or install through GitHub
```
> install.packages("SCMarker_2.0.tar.gz",repos=NULL,type="source")
> library(SCMarker)
```
Step 3: in R command window, run the following command to filter genes(cells) that expressed (non zero) distribution is similar with normal distribution
```
> Res <- ModalFilter(data = data, geneK = 10, cellK = 10)
```
Step 4: in R command window, run the following command to filter genes which are widely expressed in most of the cell population
```
> Res <- GeneFilter(obj = Res)
```
Step 5: in R command window, run the following command to select markers if gene pairs are mutual nearest coexpressed or exclusivity in a number of cells
```
> Res <- getMarker(obj = Res, k = 300, n = 30)
```
Step 6: At the end of the help page, run the following command to obtain the selected genes
```
> scMarker_genes <- Res$marker
```


# M3Drop
## Description
M3Drop is used to model the dropout rate of the genes as a function of the mean expression level using the Michaelis-Menten equation. The gene-wise Michaelis-Menten constants are computed and log-transformed, and the genes are then ranked by their p-value from a Z-test
comparing the gene-wise constants to a global constant obtained by combining all the genes.
Andrews, T. "M3Drop: Michaelis-Menten Modelling of Dropouts in single-cell RNASeq." R package version 1.0 (2019).

## Usage
Note: The input data is a gene expression matrix, gene in rows and cells in columns. 

Step 1: open your R or Rstudio
Step 2: in the R command window, run the following command to install and load the R package
```
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("M3Drop")
> library(M3Drop)
```
Step 3: in R command window, run the following command to filter and normalize the given expression matrix. Removes low quality cells and undetected genes, and normalizes counts to counts per million
```
> Normalized_data <- M3DropCleanData(data, labels = label, is.counts = FALSE)
```
Step 4: in R command window, run the following command to find differentially expressed (DE) genes by using Michaelis-Menten curve
```
> DE_genes <- M3DropFeatureSelection(Normalized_data$data, mt_method = "fdr", mt_threshold = 0.05)
```
Step 5: in R command window, run the following command to sort the genes according to the p-value
```
> DE_genes = DE_genes[order(DE_genes$p.value, decreasing = TRUE),]
```
