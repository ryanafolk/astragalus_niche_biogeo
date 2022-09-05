library(ape)
library("phytools")
tree <- read.tree("/mnt/Heuheu/nitfix/astragalus/BAMM/astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.tre")
is.ultrametric(tree) 
tree2 <- force.ultrametric(tree, method="extend")
is.ultrametric(tree2) 
is.binary.tree(tree)
# Now to check min branch length; must not be zero
min(tree$edge.length)

write.tree(tree2, file= "astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.forcedultra.tre", append = FALSE, digits = 10, tree.names = FALSE)

tree2 <- read.tree("/mnt/Heuheu/nitfix/astragalus/BAMM/astragalus.nodups.forcedultra.astragalusonly.tre")

#set priors
install.packages("BAMMtools")
library(BAMMtools) 
setBAMMpriors(tree2)

## Use this version for a trait analysis -- don't attempt to load trait file as an R object
#setBAMMpriors(phy = tree2, traits = "niche_phylopca_PC1.txt")
#
#tree2 <- read.tree("ultrametric_occur_trait_matched_forcedultra.tre")
#
#setBAMMpriors(phy = tree2, traits = "trait_phylomds_PC1.txt")