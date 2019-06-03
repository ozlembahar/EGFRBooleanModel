#### Attractor and robustness analysis of Boolean model of the EGFR signaling#####

# Load packages needed for simulations 
require(igraph)
require(BoolNet)
require(readxl)

# Read rules from file
model = loadNetwork("EgfrBoolean.txt", symbolic = F) 

# Double-check gene names
mapping = AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = model$genes, keytype = "SYMBOL", columns = "ENTREZID")
mapping$SYMBOL[is.na(mapping$ENTREZID)] # Only non-gene processes should be here

#Generate a starting state
inputs = sapply(model$interactions, "[[", "expression") == sapply(model, names)$interactions 
path.starting.state = (BoolNet::getPathToAttractor(model, as.numeric(inputs)))

#Get gene names for the nodes in the model
model.info = data.frame(read_excel('EgfrMap.xlsx', col_types = "text"), row.names = 3, stringsAsFactors = F)
names(inputs) = model$genes
model.geneids = model.info[model$genes, "Gene.ID..Entrez."]


# Reorder the nodes for tidy plotting
phen = c("Apoptosis", "Proliferation" ,"Neurogenesis")
comp = c("Mtorc1", "Mtorc2")
chem = "diacylglycerol"
root = "Egf"
gene = setdiff(model$genes, c(phen, comp, chem, root))
node.order = c(root, gene, comp, chem, phen)




################################### Starting state = 0 ########################################

#Generate a starting state in which the input node also is 0

start.state <- getPathToAttractor(model, rep(0, 53))
par(mar = c(1, 10, 2.5, 2))
plotSequence(sequence= start.state, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "An attractor of an un-stimulated cell")


#Simulate when starting conditions set to 0


par(mar = c(1, 10, 2.5, 2))
path.results = plotSequence(sequence = start.state,  offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "An attractor of an un-stimulated cell" )


############################## Proliferation ##################################################
# Simulate proliferation conditions by setting Egf to 1
#assumption: Prcka which is inhibitor of Egfr is initially inactive

ko.state = model
ko.state$interactions$Prkca$func = rep(0, length(ko.state$interactions$Prkca$func))
ko.state$interactions$Egf$func = rep(1, length(ko.state$interactions$Egf$func))
attractorpath = getPathToAttractor(ko.state, as.numeric(inputs))
par(mar = c(1, 10, 2.5, 2))
attractor.results = plotSequence(sequence = attractorpath,  offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Egf=1" )


##################### Robustness analysis ###############################
# Transition robustness is measured.
r  <- perturbTrajectories(model, measure="sensitivity",  numSamples=1000, flipBits=1) # measure=sensitivity,
r$value

#perturbedNet <- perturbNetwork(ko.state,  perturb="functions", method="bitflip")

#Attractor robustness can be measured by using testAttractorRobustness function due to the large scale of the map
#res <- testNetworkProperties(model, numRandomNets=100, testFunction="testAttractorRobustness", testFunctionParams=list(perturb="functions", numSamples=10), accumulation = "characteristic", alternative= "less")

#Calculate normalized hamming distance:the number of genes which differ between the original successor state and the successor state of the perturbed.
# If the distance is zero, the mutation has no effect on evaluated network behavior
par(mar = c(6, 4.1, 4.1, 2.1))
res <- testNetworkProperties(ko.state, numRandomNets=1000, testFunction="testTransitionRobustness", testFunctionParams = list(numSamples=100), alternative= "less")

res$pval
res$significant



##### Knock-outs induced to proliferating cells

ko.state$interactions$Jun$func = rep(0, length(ko.state$interactions$Jun$func))



############################# Differentiation #############################################
#Note: Skip this step and go to next one for attractor analysis where an initial condition is the stable state of proliferating cell and simulation was done upon Egf withdrawal
#Simulate differentiation conditions by setting Egf to 0
#assumption: Prcka which is the inhibitor of Egfr gets activated for inhibition of the excess amount of functional Egfr

ko.state = model

#Knock-out Egf for activation of neurogenesis
ko.state$interactions$Egf$func = rep(0, length(ko.state$interactions$Egf$func))

#Components need to be activated for activation of neurogenesis
ko.state$interactions$Prkca$func = rep(1, length(ko.state$interactions$Prkca$func))

#Attractor analysis
attractorpath = getPathToAttractor(ko.state, as.numeric(inputs))

par(mar = c(1, 10, 2.5, 2))
attractor.results = plotSequence(sequence = attractorpath, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Egf=0" )

########### Use the attractor of proliferation as a starting state and then induce differentiation ################################

# Attractor of a proliferating cell at steady state was used as a starting state
# To induce differentiation, Egf was set to 0 at time step 1

ko.state$interactions$Egf$func = rep(0, length(ko.state$interactions$Egf$func))
ko.state$interactions$Prkca$func = rep(1, length(ko.state$interactions$Prkca$func))
diffpath = getPathToAttractor(ko.state, attractorpath[9, 1:53]) # attractorpath of proliferation conditions
par(mar = c(1, 10, 2.5, 2))
diff.attractor = plotSequence(sequence = diffpath, offColor = "#CCCCCC", onColor = "#006633", reverse = F, title = "Egf=0" )




