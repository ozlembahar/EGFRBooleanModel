source("http://research.rupertoverall.net/include.R")

require(igraph)
require(BoolNet)
require(readxl)

# Read rules from file
model = loadNetwork("EgfrBoolean.txt", symbolic = F) 

# Double-check gene names
mapping = AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = model$genes, keytype = "SYMBOL", columns = "ENTREZID")
mapping$SYMBOL[is.na(mapping$ENTREZID)] # Only non-gene processes should be here

# Visually check sanity of loaded model (and return iGraph object)
graph = plotNetworkWiring(model, layout_with_fr, grid = "grid", vertex.size = 10)
layout = layout_as_tree(graph, root = "Egf", circular = FALSE)
plot(graph, layout = layout)

# Load discretised expression data from cell culture (DC)
data = read.delim("DC_na_discretised_data.txt", row.names = 1)
# Get the data for the nodes in the model
model.info = data.frame(read_excel('EgfrMap.xlsx', col_types = "text"), row.names = 3, stringsAsFactors = F)
model.geneids = model.info[model$genes, "Gene.ID..Entrez."]
model.data = data[model.geneids, ]
rownames(model.data) = model$genes
####
#
# Set node types for pretty display (should find a more robust way of doing this)
gene.order = setdiff(c(which(model$genes == "Egf"), which(model$genes == "Egfr"), which(model$genes != "Egf" & model$genes != "Egfr")), which(model$genes %in% c("Apoptosis", "CellGrowth", "CellCycleProgression", "diacylglycerol")))
grouping = list(class = c("phen", "chem", "comp", "gene", ""), index = list("phen" = which(model$genes %in% c("Apoptosis", "CellGrowth", "CellCycleProgression")), "chem" = which(model$genes %in% c("diacylglycerol")), "comp" = which(model$genes[gene.order] %in% c("Mtorc1", "Mtorc2")), "gene" = which(model$genes %in% rownames(model.info))))
source("plotSequence.R") # Override function from BoolNet
##

# Proliferating
# model.data$X0 is the data for t = 0 (proliferation conditions)
proliferating.state = as.numeric(model.data$X0)
missing = is.na(proliferating.state)
names(proliferating.state) = model$genes
proliferating.state[missing] = rep(1, length(which(missing))) # Set nodes without expression data to 'on'

###
# Knockouts
# Here is a way to force a knockout regardless of where it is in the network (set its function to be always false)
ko.model = model
proliferating.state[is.na(proliferating.state)] = 0 # NAs switched off

##Set initial conditions: 
ko.model$interactions$Pdpk1$func = rep(1, length(ko.model$interactions$Pdpk1$func))
ko.model$interactions$Egf$func = rep(1, length(ko.model$interactions$Egf$func))
ko.model$interactions$Akt1$func = rep(1, length(ko.model$interactions$Akt1$func))
ko.model$interactions$Nfkb1$func = rep(0, length(ko.model$interactions$Nfkb1$func))

#plot attractor of proliferation conditions
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions" )
proliferating.attractor = path.results[, ncol(path.results)]
#
### Hypothesis 1: c-JUN make up components of the AP-1 (activator protein 1) complex, which functions as a transcription factor that binds to the CYCLIN D1 promoter to activate transcription of the important G1 driver
### In differentiating cells, we expect an increase in the concentration of Jun(needed for induction of proliferation)
## If Jun is knocked-out, no cellcycleprogression
#ko.model$interactions$Jun$func = rep(0, length(ko.model$interactions$Jun$func))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferation Conditions t=0")
###
#Possible attractors showing prolifetation state found by exhaustive search
knockedOut<- fixGenes(model, c("Jun", "Nfkb1", "Akt1"), c(0,0,1))
attractors <- getAttractors(knockedOut)
plotAttractors(attractors)


#Check if Cell cycle arrest=1, what is the phemotype?
knockedOut<- fixGenes(ko.model, c("Egfr", "Egf"), c(-1, -1))
attractors <- getAttractors(knockedOut)
plotAttractors(attractors)

###Check how the attractor changes when CellCycleprogression=0
ko.model = model
proliferating.state[is.na(proliferating.state)] = 0 # NAs switched off

##-Pdpk1 set to 1 since it is always active in physilgical conditions
ko.model$interactions$Pdpk1$func = rep(1, length(ko.model$interactions$Pdpk1$func))
ko.model$interactions$Egf$func = rep(1, length(ko.model$interactions$Egf$func))
ko.model$interactions$CellCycleProgression$func = rep(0, length(ko.model$interactions$CellCycleProgression$func))### Hypothesis 1: c-JUN make up components of the AP-1 (activator protein 1) complex, which functions as a transcription factor that binds to the CYCLIN D1 promoter to activate transcription of the important G1 driver
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions, CCP=0 ")
###
# Find self loops in differentiation states
inputs = sapply(model$interactions, "[[", "expression") == sapply(model, names)$interactions
t(BoolNet::getPathToAttractor(model, as.numeric(inputs)))


#Differentiation state
differentiated.state = as.numeric(model.data$X6)
missing = is.na(differentiated.state)
names(differentiated.state) = model$genes
differentiated.state[missing] = rep(0, length(which(is.na(differentiated.state)))) # Set nodes without expression data to 'off'

# Starting states of differentiation condition
ko.model$interactions$Pdpk1$func = rep(1, length(ko.model$interactions$Pdpk1$func))
ko.model$interactions$Mtorc2$func = rep(1, length(ko.model$interactions$Mtorc2$func))
ko.model$interactions$Nfkb1$func = rep(1, length(ko.model$interactions$Nfkb1$func))
ko.model$interactions$Egf$func = rep(0, length(ko.model$interactions$Egf$func))
ko.model$interactions$Cdk4$func = rep(0, length(ko.model$interactions$Cdk4$func))
path = getPathToAttractor(ko.model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Differentiation conditions, t=12h")
#
##Attractors showing possible differentiation state found by exhaustive search
attractors <- getAttractors(ko.model)
plotAttractors(attractors)

#How the system behaves when Nfkb1 is knocked-out
ko.model$interactions$Nfkb1$func = rep(0, length(ko.model$interactions$Nfkb1$func))
#
##Check the predictions upon knocking-out TFs in proliferation conditions
ko.model = model
proliferating.state[is.na(proliferating.state)] = 0 # NAs switched off
##1-Pdpk1 set to 1 since it is always active in physilgical conditions


ko.model$interactions$Rac1$func = rep(1, length(ko.model$interactions$Rac1$func))### Hypothesis 1: c-JUN make up components of the AP-1 (activator protein 1) complex, which functions as a transcription factor that binds to the CYCLIN D1 promoter to activate transcription of the important G1 driver

path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Knocked-out TFs in proliferating conditions")
###
#
#Check other time points if they mimic differentiation conditions
differentiated.state = as.numeric(model.data$X6)
missing = is.na(differentiated.state)
names(differentiated.state) = model$genes
differentiated.state[missing] = rep(1, length(which(missing))) # Set nodes without expression data to 'on'
#plot attractor of differentiation conditions
ko.model$interactions$Egf$func = rep(0, length(ko.model$interactions$Egf$func))
ko.model$interactions$Mtorc2$func = rep(1, length(ko.model$interactions$Mtorc2$func))
ko.model$interactions$Nfkb1$func = rep(1, length(ko.model$interactions$Nfkb1$func))
ko.model$interactions$Jun$func = rep(1, length(ko.model$interactions$Jun$func))
ko.model$interactions$Stat3$func = rep(1, length(ko.model$interactions$Stat3$func))
path = getPathToAttractor(ko.model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Differentiation conditions t=6h" )



#Differentiation state
differentiated.state = as.numeric(model.data$X96)
missing = is.na(differentiated.state)
names(differentiated.state) = model$genes
differentiated.state[missing] = rep(0, length(which(is.na(differentiated.state)))) # Set nodes without expression data to 'off'
path = getPathToAttractor(model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Differentiation conditions, t=96h")
#

#
###
### 2: Upon activation, mTORC1 signaling drives cell growth and cell proliferation 
### Experiment 1: Fix Mtorc1 and Egf to 1 and plot attractors 
ko.model$interactions$Mtorc1$func = rep(1, length(ko.model$interactions$Mtorc1$func))
ko.model$interactions$Egf$func = rep(1, length(ko.model$interactions$Egf$func))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions, mTORC1=1")

#Hypothesis 1: Upon activation, mTORC1 signaling drives cell growth and cell proliferation
#Hypothesis 2: Mtorc2 signaling regulates apoptosis
#Experiment 2: To test both hypotheses 1 and 2, fix Mtorc2 to 1 and plot attractors
ko.model$interactions$Mtorc2$func = rep(1, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(2, 10, 0.7, 5))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Egf=1, Mtorc2=1")
 

#Hypothesis 3: Both mTOR complexes are crucial for the proper neuronal cell growth
#Experiment 3: Fix Mtorc1 to 1, Mtorc2 to 0 
ko.model$interactions$Mtorc1$func = rep(1, length(ko.model$genes))
ko.model$interactions$Mtorc2$func = rep(0, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions, mTORC1=1")


#Hypothesis 3: Both mTOR complexes are crucial for the proper neuronal cell growth
#Experiment 4: Fix Mtorc1 to 0 , Mtorc2 to 1  
ko.model$interactions$Mtorc1$func = rep(0, length(ko.model$genes))
ko.model$interactions$Mtorc2$func = rep(1, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions, mTORC1=1")

#Hypothesis 3: Both mTOR complexes are crucial for the proper neuronal cell growth
#Experiment 5: Fix Mtorc1 to 0 , Mtorc2 to 0
ko.model$interactions$Mtorc1$func = rep(0, length(ko.model$genes))
ko.model$interactions$Mtorc2$func = rep(0, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions, mTORC1=1")





# 1.Is my model a good model?
# 1.a. Given the binary states for the input variables at t=0, determine the resulting logical steady state. What is the network's response?
# 1.b. How can it be verified that the system behaves consistently at multiple time points?
# Define steady state at t = 0

path = getPathToAttractor(model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions" )
proliferating.attractor = path.results[, ncol(path.results)]
#
#
#
#Perturbed <- perturbNetwork(model, perturb = "transitions", method = "shuffle", excludeFixed = FALSE)
#
#
# 1.c. If you need experimental evidence - what experiments would you need? 
# Check if the model has fragile points where a mutated protein may support uncontrolled growth, proliferation, cell growth
# Set every single gene to 0 at t = 0 to find candidate gene(s) involved in proliferation
# Set every single gene to 0 at t = 96 to determine candidate gene(s) regulating differentiation, apoptosis
# Knock-out genes to perturb the system
# Are there any interesting/unexpected outcomes we can also test experimentally?
#
#Set every single gene to 0 at t = 0
proliferating.state[""] = 1 # Set special condition: (in proliferation conditions, EGF is always present)
path = getPathToAttractor(model, proliferating.state)
par(mar = c(2, 10, 1, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Proliferating conditions" )
proliferating.attractor = path.results[, ncol(path.results)]
#
#
#
## 
#Differentiation state
differentiated.state = as.numeric(model.data$X96)
missing = is.na(differentiated.state)
names(differentiated.state) = model$genes
differentiated.state[missing] = rep(0, length(which(is.na(differentiated.state)))) # Set nodes without expression data to 'off'
#
path = getPathToAttractor(model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Differentiation conditions")
#
#
#
##############Find self loops set their value to 1 others to 0 ###########################
model = loadNetwork("EgfrBoolean.txt")
inputs = sapply(model$interactions, "[[", "expression") == sapply(model, names)$interactions
t(BoolNet::getPathToAttractor(model, as.numeric(inputs)))


##Attractors of Egf withdrawal in differentiation conditions should be the same with the attactor of differentiation condition
##################################################################################################################################
ko.model$interactions$Egf$func = rep(0, length(ko.model$interactions$Egf$func))
path = getPathToAttractor(inputs, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Differentiation Conditions, Egf=0")


#
#
#2) What can your model predict?
#2.a. If you perturb conditions in your model, what information can you get out? 
#2.a.1. Determine which conditions should be perturbed and what you expect to detect 
#
#
#
#
#
#
#
#
#2.b. Can you predict which genes are essential for certain phenotypes? Can you predict the effects of perturbation on a phenotype of interest? Is it likely that you will be able to observe/measure these predicted differences if you perturb the real biological system?
#2.b.1. Set CellCycleProgression (CCP) to 1 and check which genes are active at steady state. Then set CCP to 0. Determine which genes might be needed for cell cycle progression. (The ones expressed when CCP=1 and not expressed when CCP=0 will be determined bearing in mind that genes which are not essential for CCP might be also included.) 
#Given that the activation of Cdk4/6-Ccnd1 complex drive cells though G1 leading cell cycle progression, check whether cyclins will be also activated.
#Find out which time point would be the best to use to reveal candidate genes involved in CCP? 
#CCP=1
ko.model$interactions$CellCycleProgression$func = rep(1, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(0, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Cell cycle progression=1, Egf=0")

#
#CCP=0
ko.model$interactions$CellCycleProgression$func = rep(1, length(ko.model$genes))
ko.model$interactions$Egf$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Cell cycle progression=1, Egf=1")

#2.b.2. Cyclin d1 and Myc are regulators of cell cycle progression.  
ko.model$interactions$Ccnd1$func = rep(0, length(ko.model$genes))
ko.model$interactions$Myc$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, proliferating.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Ccnd1=0, Myc=0")

#2.b.3.The re-expression of cell cycle markers such as Myc has been linked with cell death in certain types of neuronal cells
#Check whether the activation of apoptosis can be detected when Myc is over-expressed at t96

#Differentiation state
differentiated.state = as.numeric(model.data$X96)
missing = is.na(differentiated.state)
names(differentiated.state) = model$genes
differentiated.state[missing] = rep(0, length(which(is.na(differentiated.state)))) # Set nodes without expression data to 'off'

# Define steady state at t = 96
path = getPathToAttractor(model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Differentiation conditions" )


ko.model = model
differentiated.state[is.na(differentiated.state)] = 0 # NAs switched off

#2.b.4. When Myc is over-expressed at differentiation state and cells are forced to die. Check whether apoptosis gets activated.
ko.model$interactions$Ccnd1$func = rep(0, length(ko.model$genes))
ko.model$interactions$Myc$func = rep(1, length(ko.model$genes))
path = getPathToAttractor(ko.model, differentiated.state)
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "Myc=1 at t96")

#2.b.5. Apoptosis (if the aforementioned set-up makes sense to reveal genes that are essential for CCP, the same method can also be used for Apoptosis )

# Check if Apoptosis is regulated by Mtorc2
#model$interactions$Apoptosis$func = rep(0, length(model$genes))
model$interactions$Mtorc2$func = rep(1, length(model$genes))
model$interactions$Egf$func = rep(0, length(model$genes))
path = getPathToAttractor(model, proliferating.state)
par(mar = c(2, 10, 0.7, 5))
plotSequence(sequence = path, offColor = "#F0F0F0", onColor = "#2B2B2B", reverse = F, title = "" )

#
#
#
# 
# 

########################################################################################################
# Is attractor in the absence of EGF the same as the data at 96 h?
attractor.comparison = cbind("attractor" = attractor[rownames(model.data)], model.data * 1)
eqcor(attractor.comparison) # How many of the states are the same?
barplot(eqcor(attractor.comparison)[1,])
plot.heatmap2(attractor.comparison, colour.palette = colour.schemes$warmred.2, cex = 2, y.label.width = 20, x.label.height = 0, x.label.adjust = -2)
#
#
#
#
#
#
#
proliferating.state.matrix[, "Egfr"] = 1

## When Egf = 0, t = 0, 6, 12, 24, 48, 96h
knockedOut <- fixGenes(model, "Egf", 0)
attractors <- getAttractors(knockedOut)


## When Egf = 1, t =0, 6, 12, 24, 48, 96h
overExpressed <- fixGenes(model, "Egf", 1)
attractors <- getAttractors(overExpressed)

# Run simulation using the proliferating culture data as starting state
path = getPathToAttractor(overExpressed, proliferating.state)
par(mar = c(2, 10, 1, 2))
path.results = plotSequence(sequence = path)
proliferating.attractor = path.results[, ncol(path.results)]
plotSequence(sequence = path)


## Single attractor analysis (unknowns are off)
# Run simulation using the proliferating culture data as starting state
proliferating.state["Egf"] = 1 # Set special condition: (in proliferation conditions, EGF is always present)
path = getPathToAttractor(model, proliferating.state)
par(mar = c(2, 10, 1, 2))
path.results = plotSequence(sequence = path)
proliferating.attractor = path.results[, ncol(path.results)]


## Single attractor analysis (unknowns are off) - EGF withdrawal
# Run simulation using the proliferating culture data as starting state, but without EGF
proliferating.state["Egf"] = 0 # Set special condition: (in proliferation conditions, EGF is always present)
path = getPathToAttractor(model, proliferating.state)
par(mar = c(2, 10, 1, 2))
path.results = plotSequence(sequence = path)
attractor = path.results[, ncol(path.results)]


# Is attractor in the absence of EGF the same as 96 h?
attractor.comparison = cbind("attractor" = attractor[rownames(model.data)], model.data * 1)
png(file = "Withdrawal.png", width = 1200, height = 1500)
plot.heatmap2(attractor.comparison, colour.palette = colour.schemes$warmred.2, cex = 2, y.label.width = 20, x.label.height = 0, x.label.adjust = -2)
dev.off()

# Correlate 
eqcor(attractor.comparison) # How many of the states are the same?
barplot(eqcor(attractor.comparison)[1,])
#plot.heatmap2(eqcor(attractor.comparison))





# Run all possible simulations using the proliferating culture matrix as starting states
paths = apply(proliferating.state.matrix, 1, function(state) getPathToAttractor(model, state))
path.sums = Reduce('+', paths)

pdf(file = "Prob_Sim_EGFR.pdf")
plot.heatmap2(t(path.sums), colour.palette = colour.schemes$warmred.2, cex = .7)
dev.off()

##
# Export to different formats
edge.list = as.data.frame(as_edgelist(graph, names = TRUE), stringsAsFactors = F)
colnames(edge.list) = c("Source", "Target")
edge.list$Source = mapping[edge.list$Source, ]$SYMBOL
edge.list$Target = mapping[edge.list$Target, ]$SYMBOL
# Edge sign is hand-coded as the BoolNet graph strips this information
edge.list$Sign = c(-1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
write.tab(edge.list, "EGFRv13_edges.txt")

vertex.list = as.data.frame(cbind(Vertex = union(edge.list$Source, edge.list$Target), Class = "gene"), stringsAsFactors = F)
rownames(vertex.list) = vertex.list$Vertex
vertex.list[c("Mtorc1", "Mtorc2", "Raccdc42"), ]$Class  = "complex"
vertex.list$Identifier = mapping[match(vertex.list$Vertex, mapping$SYMBOL), ]$ENTREZID
write.tab(vertex.list, "EGFRv13_vertices.txt")

