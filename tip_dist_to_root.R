tree <- read.tree("tree.tre")
library(adephylo)
treetipheights <- distRoot(tree, tips = "all")
write.csv(treetipheights, "treetipheights.csv")
