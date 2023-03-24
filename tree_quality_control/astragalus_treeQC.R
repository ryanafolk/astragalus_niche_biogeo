library(phytools)
library(castor)
library(caper)
library(geiger)

tree<-read.tree("astragalus.wBranchLenghts.outliersremoved.renamed.tre")

X<-read.csv("qc_astragalus_rename.csv", header=TRUE, row.names=1)

treedata_object = treedata(tree, X)
tree.reduced <- treedata_object$phy
data.reduced <- treedata_object$data

genes<-as.factor(setNames(data.reduced[,1],rownames(data.reduced)))
paralogs<-as.factor(setNames(data.reduced[,2],rownames(data.reduced)))

pdf(file="genes.at50p.new.pdf", width = 30 , height = 90)
plot((contMap(tree.reduced,genes,plot=FALSE)),legend=0.7*max(nodeHeights(tree.reduced)),sig=2,fsize=c(0.7,0.9),leg.txt="genes at 50%")
dev.off()

pdf(file="paralogs.new.pdf", width = 30 , height = 100)
plot(setMap(contMap(tree.reduced,paralogs,plot=FALSE), invert=TRUE),legend=0.7*max(nodeHeights(tree.reduced)),sig=2,fsize=c(0.7,0.9),leg.txt="paralogs")
dev.off()
