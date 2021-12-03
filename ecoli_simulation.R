library(dbscan)
library(speccalt)
library(clusterR)
library(dirichletprocess)
library(mvnfast)
library(mclust)

data("iris")
x <- as.matrix(iris[, 1:4])
db <- dbscan(x, eps = .4, minPts = 4)

kNNdistplot(x, k = 5 - 1)
abline(h=.7, col = "red", lty=2)

ecoli <- read.table("./data/ecoli.data")

# e coli dataset without labels
data <- ecoli[, 2:8]
# labels
labels <- ecoli[,9]

# dbscan
kNNdistplot(data, k = 8)
db <- dbscan(data, eps = 0.3, minPts = 3)
db_assignments <- db$cluster
db_rand <- adjustedRandIndex(labels, db_assignments)

# speccalt
kern <- local.rbfdot(data)
s <- speccalt(kern) # I think these are the labels of the clusters it assigns
s_rand <- adjustedRandIndex(labels, s)

# dpmg
model <- DirichletProcessMvnormal(data)

# fgm
source("./CODE_original.R")
fgm <- FG_mixture(data, M=25, Iter=500, alpha=1, mm=10)
index_fg = rep(1,nrow(data)) # cluster assignments
# assignments are based on what sphere it is assigned 
L = ncol(fgm$inclusion_matrix)
for(j in 1:L) {
  index_fg[which(fgm$inclusion_matrix[,j]==1)]=j
}
fgm_rand <- adjustedRandIndex(labels, index_fg)

