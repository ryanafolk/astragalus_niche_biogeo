library(RevGadgets)
library(grDevices)

t <- processAncStates(path = "ChromEvol_simple_final_50K.tree",  labels_as_numbers = T)

pdf("anc_states_comp_numbers_50k.pdf", height = 80, width = 20)
plotAncStatesMAP(t = t, 
                 tree_layout = "rectangular",
                 node_color_as = "state_posterior",
                 node_size_as = "state",
                 tip_labels_states = T,
                 state_transparency = 0.5,
                 node_size = c(1, 10),
                 tip_states = FALSE,
                 tip_labels_size = 3, node_labels_as="state")
dev.off()
