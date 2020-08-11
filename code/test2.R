

# load R packages
library(splatter)
library(Rtsne)
library(ggplot2)
library(gplots)
library(RColorBrewer)

# simulated data
nGroup = 10
set.seed(12345)
sp = sample(30,10)
gprob = sp/sum(sp)

params = newSplatParams(nGenes = 20000,
                        batchCells = 1000,
                        group.prob = gprob,
                        de.prob = 0.3)

sim.dat <- splatSimulate(params,method = "groups", verbose = FALSE)
sim.dat <- normalize(sim.dat)
plotPCA(sim.dat, colour_by = "Group")

data = counts(sim.dat)
label = as.matrix(sim.dat$Group)
rownames(label) = colnames(data)

