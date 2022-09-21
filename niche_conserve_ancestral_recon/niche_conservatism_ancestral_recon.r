library(ape)
library(phytools)
library(geiger)
library(dplyr)

###########
# Combined data file

bio1 <- read.csv("./../environmental_data/BIOCLIM_1.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio1) <- c("species", "bio1")
bio1 <- distinct(bio1, species, .keep_all= TRUE)

bio2 <- read.table("./../environmental_data/BIOCLIM_3.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio2) <- c("species", "bio2")
bio2 <- distinct(bio2, species, .keep_all= TRUE)

bio3 <- read.table("./../environmental_data/BIOCLIM_3.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio3) <- c("species", "bio3")
bio3 <- distinct(bio3, species, .keep_all= TRUE)

bio4 <- read.table("./../environmental_data/BIOCLIM_4.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio4) <- c("species", "bio4")
bio4 <- distinct(bio4, species, .keep_all= TRUE)

bio7 <- read.table("./../environmental_data/BIOCLIM_7.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio7) <- c("species", "bio7")
bio7 <- distinct(bio7, species, .keep_all= TRUE)

bio12 <- read.table("./../environmental_data/BIOCLIM_12.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio12) <- c("species", "bio12")
bio12 <- distinct(bio12, species, .keep_all= TRUE)

bio15 <- read.table("./../environmental_data/BIOCLIM_15.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio15) <- c("species", "bio15")
bio15 <- distinct(bio15, species, .keep_all= TRUE)

bio17 <- read.table("./../environmental_data/BIOCLIM_17.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(bio17) <- c("species", "bio17")
bio17 <- distinct(bio17, species, .keep_all= TRUE)

elevation <- read.table("./../environmental_data/GTOPO30_ELEVATION.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(elevation) <- c("species", "elevation")
elevation <- distinct(elevation, species, .keep_all= TRUE)

nitrogen <- read.table("./../environmental_data/ISRICSOILGRIDS_new_average_nitrogen_reduced.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(nitrogen) <- c("species", "nitrogen")
nitrogen <- distinct(nitrogen, species, .keep_all= TRUE)

carbon <- read.table("./../environmental_data/ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(carbon) <- c("species", "carbon")
carbon <- distinct(carbon, species, .keep_all= TRUE)

ph <- read.table("./../environmental_data/ISRICSOILGRIDS_new_average_phx10percent_reduced.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(ph) <- c("species", "ph")
ph <- distinct(ph, species, .keep_all= TRUE)

coarsefragment <- read.table("./../environmental_data/ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(coarsefragment) <- c("species", "coarsefragment")
coarsefragment <- distinct(coarsefragment, species, .keep_all= TRUE)

sand <- read.table("./../environmental_data/ISRICSOILGRIDS_new_average_sandpercent_reduced.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(sand) <- c("species", "sand")
sand <- distinct(sand, species, .keep_all= TRUE)

needleleaf <- read.table("./../environmental_data/LandCover_1_Needleleaf.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(needleleaf) <- c("species", "needleleaf")
needleleaf <- distinct(needleleaf, species, .keep_all= TRUE)

deciduousbroadleaf <- read.table("./../environmental_data/LandCover_3_Deciduousbroadleaf.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(deciduousbroadleaf) <- c("species", "deciduousbroadleaf")
deciduousbroadleaf <- distinct(deciduousbroadleaf, species, .keep_all= TRUE)

herbaceous <- read.table("./../environmental_data/LandCover_6_Herbaceous.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(herbaceous) <- c("species", "herbaceous")
herbaceous <- distinct(herbaceous, species, .keep_all= TRUE)

aridity <- read.table("./../environmental_data/aridity_index_UNEP.average.csv", header = FALSE, stringsAsFactors = FALSE, row.names = NULL, sep = "\t")
colnames(aridity) <- c("species", "aridity")
aridity <- distinct(aridity, species, .keep_all= TRUE)

biogeo <- read.table("./../biogeography_manuallycorrected/BIOGEOCODING.formatted.biogeobears.treematched.curated.tsv", header = FALSE)
colnames(biogeo) <- c("species", "biogeo")
biogeo <- distinct(biogeo, species, .keep_all= TRUE)

soil <- read.table("./../soil_types/soiltype_mostprobable.mode.fixed.csv", header = FALSE)
colnames(soil) <- c("species", "soil")
soil <- distinct(soil, species, .keep_all= TRUE)

combined = merge(bio1, bio2, by = "species")
combined = merge(combined, bio3, by = "species")
combined = merge(combined, bio4, by = "species")
combined = merge(combined, bio7, by = "species")
combined = merge(combined, bio12, by = "species")
combined = merge(combined, bio15, by = "species")
combined = merge(combined, bio17, by = "species")
combined = merge(combined, elevation, by = "species")
combined = merge(combined, nitrogen, by = "species")
combined = merge(combined, carbon, by = "species")
combined = merge(combined, ph, by = "species")
combined = merge(combined, sand, by = "species")
combined = merge(combined, coarsefragment, by = "species")
combined = merge(combined, needleleaf, by = "species")
combined = merge(combined, deciduousbroadleaf, by = "species")
combined = merge(combined, herbaceous, by = "species")
combined = merge(combined, aridity, by = "species")
combined = merge(combined, biogeo, by = "species")
combined = merge(combined, soil, by = "species")

###
# Add taxonomy

oldworldnames <- read.table("./../astragalus_taxonomy/oldworld_names_final.tsv", header = FALSE, sep = "\t")
colnames(oldworldnames) <- c("group", "section", "species")
oldworldnames <- distinct(oldworldnames, species, .keep_all= TRUE)
oldworldnames$species <- gsub(" ", "_", oldworldnames$species)

newworldnames <- read.table("./../astragalus_taxonomy/newworld_names_final.tsv", header = FALSE, sep = "\t")
colnames(newworldnames) <- c("section", "species", "synonym")
newworldnames$group <- "Neoastragalus"
newworldnames$species <- gsub(" ", "_", newworldnames$species)
newworldnames$section <- gsub("sect. ", "", newworldnames$section)
library(stringr)
newworldnames$species <- str_extract(newworldnames$species, "[^_]+_[^_]+") # Remove authorities
newworldnames <- distinct(newworldnames, species, .keep_all= TRUE)

newworldnames$synonym <- NULL # not needed, need to match column numbers
allnames <- rbind(oldworldnames, newworldnames)
allnames$section <- NULL

combined = merge(combined, allnames, by = "species")


###
# Add diversification

DR <- read.csv("./../DR_rates/tip_DR.csv", header = TRUE)
DR <- distinct(DR, species, .keep_all= TRUE)
combined = merge(combined, DR, by = "species")


###
# Add soil type

soil <- read.csv("./../soil_types/soiltype_mostprobable.mode.renamed.csv", header = FALSE, sep = "\t")
colnames(soil) <- c("species", "soil")
soil <- distinct(soil, species, .keep_all= TRUE)
combined = merge(combined, soil, by = "species")


###
# Add ploidy

ploidy <- read.csv("./../nquire/ploidy.call.csv", header = TRUE)
ploidy <- distinct(ploidy, species, .keep_all= TRUE)
combined = merge(combined, ploidy, by = "species")





#######################
# Ancestral reconstruction

trait.vector = combined$bio1
names(trait.vector) <- combined$species

tree = read.tree("./../dated_trees/astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.tre")
is.ultrametric(tree)
tree = force.ultrametric(tree, method = "extend")
# For phenotypic MDS:
# tree = read.tree("intree.dated.crossvalidated.traitmatched.tre")
tree <- ladderize(tree)

trait.reduced <- treedata(tree, trait.vector)$data
trait.reduced.vector <- trait.reduced[,1]
tree.reduced <- treedata(tree, trait.vector)$phy

# MODEL COMPARISON
fitBM = fitContinuous(tree.reduced, trait.reduced)
fitOU = fitContinuous(tree.reduced, trait.reduced, model = "OU") # Make sure it is truly ultrametric even considering small rounding errors -- even a small discrepancy will cause complaints and longer runtime due to VCF optimization.
fitEB = fitContinuous(tree.reduced, trait.reduced, model = "EB")

fitBM$opt$aicc
fitOU$opt$aicc
fitEB$opt$aicc

# OU favored for annual temp, precipitation, aridity, etc.
 

# Plot mean annual temperature
figure = contMap(tree.reduced, trait.reduced.vector, plot = FALSE, model = "OU")
dev.new(width=6, height=6)
plot(setMap(figure, invert = TRUE), type = "fan", legend=0.7*max(nodeHeights(tree)), fsize=0.07, ftype="i", lwd=c(0.7,1), lwd=c(0.7,1), outline = FALSE)

# Plot annual precipitation, elevation, pH, nitrogen, aridity for color inversion (for bio15 temp seasonality use above)
figure = contMap(tree.reduced, trait.reduced.vector, plot = FALSE, model = "OU")
dev.new(width=6, height=6)
plot(figure, type = "fan", legend=0.7*max(nodeHeights(tree)), fsize=0.07, ftype="i", lwd=c(0.7,1), lwd=c(0.7,1), outline = FALSE)

# Custom colors
# figure = contMap(tree.reduced, trait.reduced.vector, plot = FALSE, model = "OU")
# figure$cols[1:length(figure$cols)] <- colorRampPalette(c("blue", "red"), space="Lab")(length(figure$cols))
# dev.new(width=6, height=6)
# plot(figure, type = "fan", legend=0.7*max(nodeHeights(tree)), fsize=0.07, ftype="i", lwd=c(0.7,1), lwd=c(0.7,1), outline = FALSE, setMap(figure,invert=TRUE))



###########################
## Niche conservatism
###########################

# Bio1
test = phylosig(tree.reduced, trait.reduced.vector, method = "lambda", test = T)
test

# Bio4; change the first line to get the rest of the variables
trait.vector = combined$bio4
names(trait.vector) <- combined$species
trait.reduced <- treedata(tree, trait.vector)$data
trait.reduced.vector <- trait.reduced[,1]
tree.reduced <- treedata(tree, trait.vector)$phy
test = phylosig(tree.reduced, trait.reduced.vector, method = "lambda", test = T)
test



# Biogeographic correlation -- ask whether clade and biogeography are independent in explaining environment

combined.normalized <- rapply(combined, scale, c("numeric","integer"), how="replace")
combined.normalized$biogeo <- as.factor(combined.normalized$biogeo)
combined.normalized$group <- as.factor(combined.normalized$group)
res.man <- manova(cbind(bio1, bio2, bio3, bio4, bio7, bio12, bio15, bio17, elevation, nitrogen, carbon, ph, sand, coarsefragment, needleleaf, deciduousbroadleaf, herbaceous, aridity) ~ group*biogeo, data = combined.normalized)
summary(res.man, tol = 0, test="Pillai")
summary.aov(res.man)

# res.lm <- lm(bio1 + bio2 + bio3 + bio4 + bio7 + bio12 + bio15 + bio17 + elevation + nitrogen + carbon + ph + sand + coarsefragment + needleleaf + deciduousbroadleaf + herbaceous + aridity ~ group*biogeo, data = combined.normalized)
# summary(res.lm)





###########################
## Diversification plots
###########################

# Plot
# Filter data to exclude multi-region species
combined %>% filter(biogeo %in% c("W", "E", "N", "S", "B", "A")) -> combined.filtered

library(ggplot2)
# DR vs biogeo
ggplot(combined.filtered, aes(x = biogeo, y = DR, fill = biogeo)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black")
# DR vs clade
ggplot(combined, aes(x = group, y = DR, fill = group)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black")


###########################
## Diversification plots
###########################

# Plot
# Filter data to exclude multi-region species
combined %>% filter(biogeo %in% c("W", "E", "N", "S", "B", "A")) -> combined.filtered

library(ggplot2)
# DR vs biogeo
ggplot(combined.filtered, aes(x = biogeo, y = DR, fill = biogeo)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black")
# DR vs clade
ggplot(combined, aes(x = group, y = DR, fill = group)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black")

# DR vs soiltype
# Filter data to keep to the most common soil types
combined %>% filter(soil %in% c("calcisols", "cambisols", "chernozems", "kastanozems", "leptosols", "luvisols", "regosols")) -> combined.filtered
ggplot(combined.filtered, aes(x = soil, y = DR, fill = soil)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))

# DR vs ploidy
combined$calls <- as.factor(combined$calls)
ggplot(combined, aes(x = calls, y = DR, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$DR), (max(combined$DR))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))


###########################
## Soil classification
###########################

combined$soil <- as.factor(combined$soil)
combined.reduced <- combined[combined$soil != "",]
# combined.newworld <- combined.reduced[combined.reduced$group == "Neoastragalus",]
# combined.oldworld <- combined.reduced[combined.reduced$group != "Neoastragalus",]
table <- prop.table(table(combined.reduced$soil))
soil_frequencies <- data.frame(state=names(table), proportion=as.numeric(table))
library(ggplot2)
ggplot(soil_frequencies, aes(x = state, y = proportion)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45))


###########################
## Ploidy
###########################

# Bio2 vs. ploidy
combined$calls <- as.factor(combined$calls)
ggplot(combined, aes(x = calls, y = bio2, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$bio2), (max(combined$bio2))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))

# Bio3 vs. ploidy
ggplot(combined, aes(x = calls, y = bio3, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$bio3), (max(combined$bio3))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))

# Other factors
ggplot(combined, aes(x = calls, y = bio1, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$bio1), (max(combined$bio1))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))
ggplot(combined, aes(x = calls, y = bio12, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$bio12), (max(combined$bio12))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))
ggplot(combined, aes(x = calls, y = aridity, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$aridity), (max(combined$aridity))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))
ggplot(combined, aes(x = calls, y = elevation, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$elevation), (max(combined$elevation))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))
ggplot(combined, aes(x = calls, y = nitrogen, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$nitrogen), (max(combined$nitrogen))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))
ggplot(combined, aes(x = calls, y = carbon, fill = calls)) + geom_violin(trim = TRUE) + ylim(min(combined$carbon), (max(combined$carbon))) + scale_fill_brewer(palette="BrBG") + theme(legend.position="none") + geom_boxplot(width=0.1, color="black") + theme(axis.text.x = element_text(angle = 45))

# Ask whether the effect of ploidy and clade are independent
combined.normalized <- rapply(combined, scale, c("numeric","integer"), how="replace")
combined.normalized$group <- as.factor(combined.normalized$group)
combined.normalized$call <- as.factor(combined.normalized$call)

res.man <- manova(cbind(bio1, bio2, bio3, bio4, bio7, bio12, bio15, bio17, elevation, nitrogen, carbon, ph, sand, coarsefragment, needleleaf, deciduousbroadleaf, herbaceous, aridity) ~ group*call*biogeo, data = combined.normalized)
summary(res.man, tol = 0, test="Pillai")
summary.aov(res.man)



###########################
## FiSSE
###########################


library(ape)
library(phangorn)
library(diversitree)
library(geiger)
library(phytools)
library(stringr)


# setwd("/mnt/Botbot/csnsclch/nitfix_pgls_new/fisse")
# getwd()

source("traitDependent_functions.R") # Get from https://github.com/macroevolution/fisse

# Check for ultrametricity and if binary; FiSSE has its own function for ultrametricity but this should be faster and fine with slight precision errors 
if(is.binary(tree.reduced)) {
	print("TRUE")
	} else {
	tree <- multi2di(tree.reduced)
	}

if(is.ultrametric(tree.reduced)) {
	print("TRUE")
	} else {
	tree <- force.ultrametric(tree.reduced, method="extend")
	}

# West Asia vs. all else
combined.filtered$biogeo_binary <- combined.filtered$biogeo
combined.filtered$biogeo_binary <- gsub("E", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("N", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("S", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("B", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("A", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("W", "1", combined.filtered$biogeo_binary)

traits <- as.numeric(combined.filtered$biogeo_binary)
names(traits) <- as.character(combined.filtered$species)

treedata_object <- treedata(tree.reduced, traits)
tree.fisse <- treedata_object$phy

traits <- traits[tree.fisse$tip.label]

res <- FISSE.binary(tree, traits)

pval_1tailed <- min(res$pval, 1-res$pval)
pval_1tailed

# Americas vs. all else
combined.filtered$biogeo_binary <- combined.filtered$biogeo
combined.filtered$biogeo_binary <- gsub("E", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("N", "1", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("S", "1", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("B", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("A", "0", combined.filtered$biogeo_binary)
combined.filtered$biogeo_binary <- gsub("W", "0", combined.filtered$biogeo_binary)

traits <- as.numeric(combined.filtered$biogeo_binary)
names(traits) <- as.character(combined.filtered$species)

treedata_object <- treedata(tree.reduced, traits)
tree.fisse <- treedata_object$phy

traits <- traits[tree.fisse$tip.label]

res <- FISSE.binary(tree, traits)

pval_1tailed <- min(res$pval, 1-res$pval)
pval_1tailed

