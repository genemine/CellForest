
# 1. Cell Forest
## Description
Cell Forest is a gene selection method based on single-cell RNA sequencing data that is aimed at selecting biologically meaningful subsets of genes through predicted label and random forests. Based on a large-scale analysis of 20 real world scRNA-seq datasets that cover a wide spectrum of biological scenarios, it was found that Cell Forest achieved, on average, superior performances to four gene selection strategies inlcuding Expr, HVG, SCMarker, and M3Drop. 

## Download
Cell Forest is implemented as an R package, which is freely available for non-commercial use. 

Version Changes 
[CellForest_1.0.0.tar.gz](https://github.com/genemine/CellForest/blob/master/CellForest_1.0.0.tar.gz)

# 2. Install

- Step 1: Firstly, install the dependent packages: **randomForest**.
```
> install.packages("randomForest")
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
Step 4: At the end of the help page, there is an example code. Copy these codes to command to run as follows:
```
data(CFDemo)
result = CellForest(data,kcluster = kprior)
```

# 4. Contact
If any questions, please do not hesitate to contact us at: 

Hongdong Li, hongdong@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn


# 5. How to cite?
If you use this tool, please cite the following work.

Hongdong Li, Yunpei Xu, Cuixiang Lin, Fangxiang Wu and Jianxin Wang, Cell Forest: gene selection strategy based on
single-cell RNA sequencing analysis, 2020, submitted  
