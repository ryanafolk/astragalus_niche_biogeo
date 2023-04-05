library(ape)
library(phytools)
library(geiger)

##### ploidy from nQuire ####

tree = read.tree(file="astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.mar23.tre")
X<-read.csv("ploidy.call.original.csv", row.names=1)
call<-as.factor(setNames(X[,1],rownames(X)))
co <- c("purple","red","pink")

treedata_object = treedata(tree, X)
tree.reduced <- ladderize(treedata_object$phy)
data.reduced <- as.factor(treedata_object$data)

er.call <- ace(data.reduced, tree.reduced, type="d", model="ER")
pdf(file="rec.ploidy.ER.mar23.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3)
nodelabels(pie=er.call$lik.anc, piecol=co, cex=0.2)
tiplabels(pie=to.matrix(call[tree.reduced$tip.label],levels(call)),piecol=co,cex=0.05)
legend("topleft",levels(data.reduced),pch=21,pt.bg=co,pt.cex=2.2)
dev.off()
aic.er.call <- AIC(er.call)  #  1967.5

ard.call <- ace(data.reduced, tree.reduced, type="d", model="ARD")
pdf(file="rec.ploidy.ARD.mar23.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3)
nodelabels(pie=ard.call$lik.anc, piecol=co, cex=0.2)
tiplabels(pie=to.matrix(call[tree.reduced$tip.label],levels(call)),piecol=co,cex=0.05)
legend("topleft",levels(call),pch=21,pt.bg=co,pt.cex=2.2)
dev.off()
aic.ard.call <- AIC(ard.call) # 1874.789

#########

###### FAO soil types ######

tree2 = read.tree(file="astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.newtrim23.tre")
Y<-read.csv("soiltype_mostprobable.mode.renamed.csv", row.names=1)
treedata2 <- treedata(tree2, Y)
tree.reduced2 <- ladderize(treedata2$phy)
data.reduced2 <- treedata2$data


c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1","skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6", "#FDBF6F", "gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

soil<-as.factor(setNames(Y[,1],rownames(Y)))


### ER model ###
er.soil <- ace(data.reduced2, tree.reduced2, type="d", model="ER")
pdf(file="rec.soil.ER.mar23.pdf", width = 20, height = 20)
plot(tree.reduced2, type="fan", cex = 0.3)
nodelabels(pie=er.soil$lik.anc, piecol=c25, cex=0.2)
tiplabels(pie=to.matrix(soil[tree.reduced2$tip.label],levels(soil)),piecol=c25,cex=0.05)
legend("topleft",levels(soil),pch=21,pt.bg=c25,pt.cex=2.2)
dev.off()
aic.er.soil <- AIC(er.soil) # 4489.094

### ARD model ###
# Can't plot: Error in floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], col = piecol) : 
# floating.pie: x values must be non-negative
ard.soil <- ace(data.reduced2, tree.reduced2, type="d", model="ARD")
pdf(file="rec.soil.ARD.pdf", width = 20, height = 20)
plot(tree.reduced2, type="fan", cex = 0.3)
nodelabels(pie=ard.soil$lik.anc, piecol=c25, cex=0.2)
tiplabels(pie=to.matrix(soil[tree.reduced2$tip.label],levels(soil)),piecol=c25,cex=0.05)
legend("topleft",levels(soil),pch=21,pt.bg=c25,pt.cex=2.2)
dev.off()
aic.ard.soil <- AIC(ard.soil)   # 4186.876


#AIC table
trait <- c("Ploidy" ,"Soil type")
ER.model <- c(aic.er.call, aic.er.soil)
ADR.model <- c(aic.ard.call, aic.ard.soil)
aic.table <- data.frame(trait, ER.model, ADR.model)
write.csv(aic.table, file = "astragalus.aic.table.csv")

#         trait ER.model ADR.model
#	1    Ploidy 1967.500  1874.789
#	2 Soil type 4529.349  4186.876
