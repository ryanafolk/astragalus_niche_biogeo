setwd("~/Dropbox/postdoc_MSU/Astragalus/redoing_march2023/")
getwd()

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
#trait.reduced.vector <- data.reduced[,1]

c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1","skyblue2", "#FB9A99", "palegreen2",
  "#CAB2D6", "#FDBF6F", "gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

#write.tree(tree.reduced, file = "tree.trimmed.tre")
#write.csv(data.reduced, file = "soil.reduced.csv")
#tree3 <- read.tree(file = "tree.trimmed.tre")
#Z<-read.csv("soil.reduced.csv", row.names=1)
soil<-as.factor(setNames(Y[,1],rownames(Y)))
#name.check(tree.reduced, trait.reduced.vector, data.names=NULL)

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
aic.ard.soil <- AIC(ard.soil)  # 3873.385

class(ard.soil$lik.anc[1,1])

##### chromosome number #####
# full tree #
tree.c = read.tree(file="chrom.call.corrected.tre")
C<-read.csv("chrom.states.corrected.csv", row.names=1, header = TRUE)
chrom<-as.factor(setNames(C[,1],rownames(C)))
#name.check(tree.c, chrom, data.names=NULL)

# the tiplabel function doesn't allow NAs, so this file has NA coded as 99, for aesthetic effects only, 
# not used in the reconstruction, only plotting current states
D<-read.csv("chrom.states.corrected.tips.csv", row.names=1, header = TRUE)
chrom.level<-as.factor(setNames(D[,1],rownames(D)))


c17<-c('magenta','forestgreen','cornflowerblue', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
        'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3', "tan1",'darkgray', 
       'wheat4', '#DDAD4B', 'chartreuse', 'black','seagreen1', 'moccasin', 'mediumvioletred', 
       'seagreen','cadetblue1',"darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

### full tree ###
### ER model ####
er.chrom <- ace(chrom, tree.c, type="d", model="ER")
pdf(file="rec.chrom.ER.pdf", width = 20, height = 20)
plot(tree.c, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=er.chrom$lik.anc, piecol=c17, cex=0.2)
tiplabels(pie=to.matrix(chrom.level[tree.c$tip.label],levels(chrom.level)),piecol=c17,cex=0.07)
legend("topleft",levels(chrom),pch=21,pt.bg=c17,pt.cex=2.2)
dev.off()
aic.er.chrom <- AIC(er.chrom) # 1384.383

### ARD model ####
## can't plot pies, complex values in the matrix
# Can't plot: Error in floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], col = piecol) : 
# floating.pie: x values must be non-negative
ard.chrom <- ace(chrom, tree.c, type="d", model="ARD")
pdf(file="rec.chrom.ARD.pdf", width = 20, height = 20)
plot(tree.c, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=ard.chrom.red$lik.anc, piecol=c17, cex=0.2)
tiplabels(pie=to.matrix(chrom.level[tree.c$tip.label],levels(chrom.level)),piecol=c17,cex=0.07)
legend("topleft",levels(chrom),pch=21,pt.bg=c17,pt.cex=2.2)
dev.off()
aic.ard.chrom <- AIC(ard.chrom) #1586.576


# clean NAs #
# run the same analyses but with a clean file without NAs

E <-read.csv("chrom.states.corrected.clean.csv", row.names=1, header = TRUE)
E<-as.factor(setNames(E[,1],rownames(E)))
chrom.match <- treedata(tree.c, E)
tree.reduced <- chrom.match$phy
data.reduced <- chrom.match$data

write.tree(tree.reduced, file = "tree.for.chromevol.tre")
write.csv(data.reduced, file = "data.chromevol.csv")

### ER model ###
er.chrom.red <- ace(data.reduced, tree.reduced, type="d", model="ER")
pdf(file="rec.chrom.clean.ER.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=er.chrom.red$lik.anc, piecol=c17, cex=0.3)
tiplabels(pie=to.matrix(E[tree.reduced$tip.label],levels(E)),piecol=c17,cex=0.1)
legend("topleft",levels(E),pch=21,pt.bg=c17,pt.cex=2.2)
dev.off()
aic.er.chrom.clean <- AIC(er.chrom.red) # 392.3815

### ARD model ###
# Can't plot: Error in floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], col = piecol) : 
# floating.pie: x values must be non-negative

ard.chrom.red <- ace(data.reduced, tree.reduced, type="d", model="ARD")
pdf(file="rec.chrom.clean.ARD.pdf", width = 20, height = 20)
plot(tree.reduced, type="fan", cex = 0.3, legend=TRUE)
nodelabels(pie=ard.chrom.red$lik.anc, piecol=c17, cex=0.3)
tiplabels(pie=to.matrix(E[tree.reduced$tip.label],levels(E)),piecol=c17,cex=0.1)
legend("topleft",levels(E),pch=21,pt.bg=c17,pt.cex=2.2)
dev.off()
aic.ard.chrom.clean <- AIC(ard.chrom.red) # 954.5484


#AIC table
trait <- c("Ploidy" ,"Soil type" ,"Chromosome count" ,"Chromosome count clean")
ER.model <- c(aic.er.call, aic.er.soil, aic.er.chrom, aic.er.chrom.clean)
ADR.model <- c(aic.ard.call, aic.ard.soil, aic.ard.chrom, aic.ard.chrom.clean)
aic.table <- data.frame(trait, ER.model, ADR.model)
write.csv(aic.table, file = "astragalus.aic.table.csv")

# trait  ER.model ADR.model
# 1                 Ploidy 1602.4330 1521.6033
# 2              Soil type 4153.7228 3873.3846
# 3       Chromosome count 1384.3831 1586.5759
# 4 Chromosome count clean  392.3815  954.5484
