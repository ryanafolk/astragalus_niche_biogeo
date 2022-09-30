

library(ape)
library(phytools)
library(geiger)

##### ploidy from nQuire ####

tree = read.tree(file="astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.tre")
X<-read.csv("ploidy.call.csv", row.names=1)
call<-as.factor(setNames(X[,1],rownames(X)))
co <- c("purple","red","pink")

er.call <- ace(call, tree, type="d", model="ER")
pdf(file="rec.ploidy.ER.pdf", width = 20, height = 20)
plot(tree, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=er.call$lik.anc, piecol=co, cex=0.2)
tiplabels(pie=to.matrix(call[tree$tip.label],levels(call)),piecol=co,cex=0.05)
legend("topleft",levels(call),pch=21,pt.bg=co,pt.cex=2.2)
dev.off()
aic.er.call <- AIC(er.call)  # 1602.433

ard.call <- ace(call, tree, type="d", model="ARD")
pdf(file="rec.ploidy.ARD.pdf", width = 20, height = 20)
plot(tree, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=ard.call$lik.anc, piecol=co, cex=0.2)
tiplabels(pie=to.matrix(call[tree$tip.label],levels(call)),piecol=co,cex=0.05)
legend("topleft",levels(call),pch=21,pt.bg=co,pt.cex=2.2)
dev.off()
aic.ard.call <- AIC(ard.call) # 1521.603

#########

###### FAO soil types ######

tree2 = read.tree(file="tree.trimmed.tre")
Y<-read.csv("soiltype_mostprobable.mode.renamed.csv", row.names=1)
treedata <- treedata(tree2, Y)
tree.reduced <- treedata$phy
data.reduced <- treedata$data
trait.reduced.vector <- data.reduced[,1]

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1","skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6", "#FDBF6F", "gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

write.tree(tree.reduced, file = "tree.trimmed.tre")
write.csv(data.reduced, file = "soil.reduced.csv")
tree3 <- read.tree(file = "tree.trimmed.tre")
Z<-read.csv("soil.reduced.csv", row.names=1)
soil<-as.factor(setNames(Z[,1],rownames(Z)))


### ER model ###
er.soil <- ace(trait.reduced.vector, tree.reduced, type="d", model="ER")
pdf(file="rec.soil.ER.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=er.soil$lik.anc, piecol=c25, cex=0.2)
tiplabels(pie=to.matrix(soil[tree.reduced$tip.label],levels(soil)),piecol=c25,cex=0.05)
legend("topleft",levels(soil),pch=21,pt.bg=c25,pt.cex=2.2)
dev.off()
aic.er.soil <- AIC(er.soil) # 4153.723

### ARD model ###
# Can't plot: Error in floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], col = piecol) : 
# floating.pie: x values must be non-negative
ard.soil <- ace(trait.reduced.vector, tree.reduced, type="d", model="ARD")
pdf(file="rec.soil.ARD.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=ard.soil$lik.anc, piecol=c25, cex=0.2)
tiplabels(pie=to.matrix(soil[tree.reduced$tip.label],levels(soil)),piecol=c25,cex=0.05)
legend("topleft",levels(soil),pch=21,pt.bg=c25,pt.cex=2.2)
dev.off()
aic.ard.soil <- AIC(ard.soil)  # 3873.385

