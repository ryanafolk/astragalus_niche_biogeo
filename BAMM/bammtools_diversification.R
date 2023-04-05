library("phytools")
library("ape")
library("BAMMtools")

# Both trait and diversification
# check convergence
mcmcout <- read.csv("astrag_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
# burnin 20% of samples -- should be fine based on plots
burnstart <- floor(0.2 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation)

# check the effective sample sizes of the log-likelihood 
# and the number of shift events present in each sample, should be larger than 200
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# the number of macroevolutionary rate regimes on our phylogenetic tree
post_probs <- table(postburn$N_shifts) / nrow(postburn)
#names(post_probs)
post_probs

# to compute the posterior odds ratio for (say) two models 
post_probs["8"] / post_probs["9"] # How much more posterior probability is in 8 shifts than 9

##Alternatively, we can summarize the posterior distribution of the number of shifts 
## using summary methods:
tree.diversification.ladderized = read.tree("/mnt/Heuheu/nitfix/astragalus/BAMM/astragalus.nodups.forcedultra.astragalusonly.tre")
# BAMM IGNORES LADDERIZATION WITHIN R, LADDERIZE FILE INSTEAD


####################################
# Diversification
####################################

edata_diversification_ladderized <- getEventData(tree.diversification.ladderized, eventdata = "astrag_event_data.txt", burnin=0.2)
shift_probs <- summary(edata_diversification_ladderized)
shift_probs
# plot.phylo(tree)

# marginal shift probablities
marg_probs <- marginalShiftProbsTree(edata_diversification_ladderized)
plot.phylo(marg_probs, lwd=0.3, cex=0.02, show.tip.label = TRUE)
branch_priors <- getBranchShiftPriors(tree.diversification.ladderized, expectedNumberOfShifts = 1)
plot(branch_priors, edge.width = 0.3, cex=0.02)
mo <- marginalOddsRatioBranches(edata_diversification_ladderized, branch_priors)


############
# Diversification plots 

##plot phylorate
# For some datasets with large numbers of taxa and rate shifts (e.g., trees with thousands of taxa), all shift configurations may have low probability. There are simply too many parameters in the model to allow a single shift configuration to dominate the credible set. An alternative approach is to extract the shift configuration that maximizes the marginal probability of rate shifts along individual branches. This is very similar to the idea of a maximum clade credibility tree in phylogenetic analysis. BAMM has a function maximumShiftCredibility for extracting this shift configuration:
#for trees with thousands of taxa
## Speciation
#msc.set <- maximumShiftCredibility(edata_diversification, maximize='product')
#msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
#plot.bammdata(best, lwd = 0.5, method ='polar', labels=TRUE, spex = "s", cex=0.01, logcolor = TRUE, breaksmethod = "jenks", legend=TRUE)
#addBAMMshifts(best, method ='polar', cex=0.4)
## add vtheta = 0.5 to plot.bammdata for no gap in circle tree, rbf = 0 for branches drawn to center

# Net diversification
best_diversification_ladderized <- getBestShiftConfiguration(edata_diversification_ladderized, expectedNumberOfShifts=1)
msc.set <- maximumShiftCredibility(edata_diversification_ladderized, maximize='product')
msc.config <- subsetEventData(edata_diversification_ladderized, index = msc.set$sampleindex)
pdf("net_diversification.pdf", width=5, height=5)
plot.bammdata(best_diversification_ladderized, lwd = 0.5, method ='polar', labels=TRUE, spex = "netdiv", cex=0.01, logcolor = TRUE, breaksmethod = "jenks", legend=TRUE)
addBAMMshifts(best_diversification_ladderized, method ='polar', cex=0.4)
dev.off()
# add vtheta = 0.5 to plot.bammdata for no gap in circle tree, rbf = 0 for branches drawn to center

## Extinction
#msc.set <- maximumShiftCredibility(edata_diversification, maximize='product')
#msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
#plot.bammdata(best, lwd = 0.5, method ='polar', labels=TRUE, spex = "e", cex=0.01, logcolor = TRUE, breaksmethod = "jenks", legend=TRUE)
#addBAMMshifts(best, method ='polar', cex=0.4)
## add vtheta = 0.5 to plot.bammdata for no gap in circle tree, rbf = 0 for branches drawn to center

##plot rate through time
ratematrix <- getRateThroughTimeMatrix(edata_diversification_ladderized) # Calculating ahead of time avoids repeating calculations to adjust figure; still need to recalculate for different nodes
#plotRateThroughTime(ratematrix,intervalCol="red", avgCol="red")
pdf("rate_through_time.pdf", width=6, height=5)
plotRateThroughTime(ratematrix,intervalCol="chartreuse4", avgCol="chartreuse4",ratetype="netdiv", ylim = c(0, 1.5))
dev.off()
#plotRateThroughTime(ratematrix,intervalCol="blue", avgCol="blue",ratetype="extinction")

#caculating Bayes Factor
#to return a pairwise matrix of Bayes Factors
bfmat <- computeBayesFactors(postburn, expectedNumberOfShifts=1, burnin=0.1)
bfmat
#

# Evolutionary rates:
allrates <- getCladeRates(edata_diversification_ladderized)
allrates
#compute the mean speciation rate for tree and estimate the 90% highest posterior density (HPD):
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))



#
#BAMM VS SOIL CAT BAMM VS BIO3 BAMM VS BIO4

#########
# Tests of trait-associated diversification

library(geiger)
# BIO2
tempfile = read.csv("/mnt/Heuheu/nitfix/astragalus/BAMM/BIOCLIM_2.average.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
lapply(trait.vector, as.numeric)
names(trait.vector) <- row.names(tempfile)

# If you get a warning about negative log values, you can put logrates = FALSE
library(parallel)
edata.subset <- subtreeBAMM(edata_diversification_ladderized, tips = names(trait.vector))
BIO2_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 
#environment_vs_speciation = traitDependentBAMM(edata, trait.vector, 1000, rate = "speciation", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 
#environment_vs_extinction = traitDependentBAMM(edata, trait.vector, 1000, rate = "extinction", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

# BIO3
tempfile = read.csv("/mnt/Heuheu/nitfix/astragalus/BAMM/BIOCLIM_3.average.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
lapply(trait.vector, as.numeric)
names(trait.vector) <- row.names(tempfile)
BIO3_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "spearman", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

# Soil
tempfile = read.csv("soiltype_mostprobable.mode.renamed.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
names(trait.vector) <- row.names(tempfile)
soil_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

# Ploidy
tempfile = read.csv("ploidy.call.renamed.txt", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
names(trait.vector) <- row.names(tempfile)
edata.subset <- subtreeBAMM(edata_diversification_ladderized, tips = names(trait.vector))
ploidy_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

# Ploidy (based on chromosomes)
tempfile = read.csv("chromosomes_cleaned.txt", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
names(trait.vector) <- row.names(tempfile)
trait.vector[trait.vector < 14] = 0
trait.vector[trait.vector >= 14] = 1
edata.subset <- subtreeBAMM(edata_diversification_ladderized, tips = names(trait.vector))
ploidy_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 

# Chromosomes
tempfile = read.csv("chromosomes_cleaned.txt", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
names(trait.vector) <- row.names(tempfile)
edata.subset <- subtreeBAMM(edata_diversification_ladderized, tips = names(trait.vector))
chromosomes_vs_diversification = traitDependentBAMM(edata.subset, trait.vector, 1000, rate = "net diversification", return.full = FALSE, method = "kruskal", logrates = TRUE, two.tailed = TRUE, traitorder = NA, nthreads = 4) 


####################
# ESSim statistic
library(mvtnorm)
source("essim.R")

tree.reduced <- treedata(tree.diversification.ladderized, tempfile)$phy

# BIO2
tempfile = read.csv("/mnt/Heuheu/nitfix/astragalus/BAMM/BIOCLIM_2.average.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
lapply(trait.vector, as.numeric)
names(trait.vector) <- row.names(tempfile)
essim(tree.reduced, trait.vector, nsim = 1000)

# BIO3
tempfile = read.csv("/mnt/Heuheu/nitfix/astragalus/BAMM/BIOCLIM_3.average.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
trait.vector = tempfile[,1]
lapply(trait.vector, as.numeric)
names(trait.vector) <- row.names(tempfile)
essim(tree.reduced, trait.vector, nsim = 1000)

####################
# FiSSE statistic

#Soil
# requires binary so we focus on calcisol vs. non-calcisol contrast
source("traitDependent_functions.R")
tempfile = read.csv("soiltype.calcisolbinary.csv", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.diversification.ladderized, tempfile)$data
tempfile <- na.omit(tempfile)
trait.vector = as.numeric(tempfile[,1])
names(trait.vector) <- as.character(row.names(tempfile))
#tree.soil <- treedata(tree.reduced, trait.vector)$phy

library(phangorn)
library(diversitree)
res <- FISSE.binary(tree.reduced, trait.vector)

# Ploidy
# requires binary so we focus on diploid vs. tetraploids
tempfile = read.csv("ploidy.call.renamed.binary.txt", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.reduced, tempfile)$data
tempfile <- na.omit(tempfile)
trait.vector = as.numeric(tempfile[,1])
names(trait.vector) <- as.character(row.names(tempfile))
tree.ploidy <- treedata(tree.reduced, trait.vector)$phy
res <- FISSE.binary(tree.reduced, trait.vector)

# Ploidy (based on chromosomes)
# requires binary so we focus on diploid vs. tetraploids
tempfile = read.csv("chromosomes_cleaned.txt", sep = "\t", row.names=1, header = FALSE)
tempfile <- treedata(tree.reduced, tempfile)$data
tempfile <- na.omit(tempfile)
trait.vector = as.numeric(tempfile[,1])
names(trait.vector) <- as.character(row.names(tempfile))
trait.vector[trait.vector < 14] = 0
trait.vector[trait.vector >= 14] = 1
tree.ploidy <- treedata(tree.reduced, trait.vector)$phy
res <- FISSE.binary(tree.reduced, trait.vector, fail_tol = 50000)



####################
## How to get marginal relative density of rates
## This follows Fig. 6 of Rabosky et al. 2014 in Sys. Bio. 
## ALWAYS check node numbering like so:
#plot(tree, show.tip.label = TRUE, cex =0.1)
#nodelabels(cex = 0.05)
## Should be identical to phytools numbering but make sure -- dropping tips would affect this
#
## Example below is 
#rate_item1 = getCladeRates(edata, node = 2200, nodetype = "include")
#rate_item2 = getCladeRates(edata, node = 2200, nodetype = "exclude")
#hist(rate_item1$beta/rate_item2$beta, xlim = c(-2, 10), breaks = 20)
#quantile(rate_item1$beta/rate_item2$beta, 0.025, na.rm = TRUE) # Lower bound of 95% HPD
#quantile(rate_item1$beta/rate_item2$beta, 0.975, na.rm = TRUE) # Upper bound of 95% HPD
#
