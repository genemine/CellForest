
# 1. Cell Forest
## Description
Cell Forest is a gene selection method based on single-cell RNA sequencing data that is aimed at selecting biologically meaningful subsets of genes through predicted label and random forests. Based on a large-scale analysis of 20 real world scRNA-seq datasets that cover a wide spectrum of biological scenarios, it was found that Cell Forest achieved, on average, superior performances to four gene selection strategies inlcuding Expr, HVG, SCMarker, and M3Drop. 

## Download
Cell Forest is implemented as an R package, which is freely available for non-commercial use. 

Version Changes 
[CellForest_2.0.0.tar.gz](https://github.com/genemine/CellForest/blob/master/CellForest_2.0.0.tar.gz)

# 2. Install

- Step 1: Firstly, install the dependent packages: **ranger**.
```
> install.packages("ranger")
```

- Step 2: Download the above CellForest package and install it in R (tested on version 3.6.2)

```
> install.packages("D:/CSU/cell_forest/20200706/CellForest_*.0.0.tar.gz", repos = NULL, type = "source")
```

# 3. Usage
Notes: CellForest was tested on linux, Mac and Windows; and it runs smoothly on these different systems.

Using CellForest is very simple. Just follow the steps below: 

Step 1: open your R or Rstudio

Step 2: in the R command window, run the following command to load the R package
```
> library(CellForest)
```

Step 3: in R command window, run the following command to see the help document for running Cell Forest. Then, you should be able to see a help page.
```
> ?CellForest
```

At the end of the help page, there is an example code. Copy these codes to command to run as follows:
Step 4: load demo data containing 191 cells and 6278 genes. The cells in the demo data come from 3 different types. The file 'cfdemo' includes a gene expression matrix (gene in rows and samples in columns) and the corresponding cell label.
```
> data(cfdemo)
> # data (6278*191): a expression data matrix,gene in rows and samples in columns
> # label (1*191): corresponding cell label
```

Step 5: Running Cell Forest function  
Parameters:
```
> # data: A expression data matrix,gene in rows and samples in columns.
> # k: The number of clusters to output.
> # ncores (default: -1): The number of cores to be used when the program running in parallel.
```
Return:
```
> # samplegeneweight: The weights of genes calculated by function randomForest based on random permuted predict label.
> # geneImportance: The weights of genes calculated by function randomForest.
> # selgenes: The genes Cell Forest final select.
> # selgeneratio: The proportion of the number of selected genes in the total genes.
> # runningtime: The running time of the program.
> # cutoff: Threshold for selecting genes.
> # plabel: Predicated label.
```
Run code:
```
> result = CellForest(data = data, k = 3)
```

Step 6: Clustering test (Hierarchical Clustering)
```
> # clustering data with all genes
> clusters1 = cutree(hclust(dist(t(data))),k=3)
> 
> # clustering data with selected genes
> clusters2 = cutree(hclust(dist(t(data[result$selgenes,]))),k=3)
```

Step 7: Evaluate clustering performance ([Evalcluster.R](https://github.com/genemine/CellForest/blob/master/code/Evalcluster.R))  
Parameters:
```
> # truelabel: A numeric vector of true labels of each sample.
> # predlabel: A numeric vector of predicted labels of each sample.
```
Return:
```
> # NMI: Value of normalized mutual information.
> # RI: Value of rand index.
> # ARI: Value of adjusted rand index.
```

Run code:
```
> evalcluster(truelabel = label, predlabel = clusters1)
> NMI        RI       ARI 
> 0.6202449 0.7133095 0.4315412 
>
> evalcluster(truelabel = label, predlabel = clusters2)
> NMI  RI ARI 
>  1   1   1 
```

# 4. Contact
If any questions, please do not hesitate to contact us at: 

Hongdong Li, hongdong@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn


# 5. How to cite?
If you use this tool, please cite the following work.

Hongdong Li, Yunpei Xu, Cuixiang Lin, Fangxiang Wu and Jianxin Wang, Cell Forest: gene selection strategy based on
single-cell RNA sequencing analysis, 2020, submitted  
