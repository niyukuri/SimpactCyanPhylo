

# Before running the code below, make sure you have followed the instructions for installing all the required software,
# as explained in the README file.


# setwd("/path/to/your/working_directory/") # Change this to your own working directory

# Make sure compiled tools (Seq-Gen and FastTree) are in same working directory

work.dir <- "/home/david/SimpactCyanPhylo"

setwd(work.dir)

library(RSimpactCyan)
library(RSimpactHelper)
library(Rcpp)
library(ape)
library(expoTree)
library(data.table)
library(readr)
library(phangorn)
library(lme4)
library(nlme)
library(dplyr)
library(adephylo)
library(treedater)
library(geiger)
library(picante)
library(igraph)
library(phyloTop)
library(phytools)
library(Rsamtools)
library(robustbase)
library(intergraph)
library(lubridate)
library(tidyr)
library(data.table)
library(ggplot2)
library(geomnet)
library(ggtree)


###########################################
# Step 1: Setup and running simpact      #
###########################################



inputvector <- c(777, -0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -3.0) # first element is the SEED (777s)



## Run Simpact for specific parameter combination

age.distr <- agedistr.creator(shape = 5, scale = 65)
#
cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 40, 
                                 population.nummen = 5000, 
                                 population.numwomen = 5000,
                                 hivseed.time = 10, 
                                 hivseed.type = "amount",
                                 hivseed.amount = 10, 
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 debut.debutage = 15
)



# # Assumption of nature of sexual network
# #########################################
#
cfg.list["population.msm"] = "no"


# # Sexual behaviour
# ###################
#
seedid <- inputvector[1]

cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)


# # HIV transmission
# ###################
#

cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)

# Disease progression > may be remove in parameter to estimates

cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)


#
# # Demographic
# ##############
#

cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)


# # Assumptions to avoid negative branch lengths
# ###############################################
# # + sampling == start ART
# # when someone start ART, he/she is sampled and becomes non-infectious

cfg.list["monitoring.fraction.log_viralload"] <- 0 # very important for transmission tree and sequence simulation

# Note: If treatment is started, the person’s set-point viral load value will be lowered according 
# to the setting in monitoring.fraction.log_viralload, if we set it to 0, it means the person is no more infectious

#
# ## Add-ons
#
### BEGIN Add-on
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["dropout.interval.dist.type"] <- "uniform"
cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6


cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45



cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
cfg.list["diagnosis.baseline"] <- -99999 # -2


#### END Add-ons


# # ART intervention
# ###################


# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly

# Introducing ART
art.intro <- list()
art.intro["time"] <- 23    # ~2000
art.intro["diagnosis.baseline"] <- -2
art.intro["monitoring.cd4.threshold"] <- 100

art.intro1 <- list()
art.intro1["time"] <- 25     # ~2002
art.intro1["diagnosis.baseline"] <- -1.8
art.intro1["monitoring.cd4.threshold"] <- 150

art.intro2 <- list()
art.intro2["time"] <- 28     # ~2005
art.intro2["diagnosis.baseline"] <- -1.5
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 33     # ~2010
art.intro3["diagnosis.baseline"] <- -1
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 36     # ~2013
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 39     # ~2016
art.intro5["monitoring.cd4.threshold"] <- 700 # 


interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist

# Events
cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3

# Avoid overlaping in same directory

#creating subfolder with unique name for each simulation
generate.filename <- function(how.long){
  
  rn <- sample(1:100,1)
  t <- as.numeric(Sys.time())
  set.seed((t - floor(t)) * 1e8)
  chars <- c(letters, LETTERS)
  sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
  
  noise.sample1 <- sample(8:15,1, replace = TRUE)
  sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
  noise.sample <- sample(1:1000,1)
  noise.sample2 <- sample(8:17,1, replace = TRUE)
  sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                           paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
  
  return(sub.dir.sim.id)
}



sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))




# Running Simpact 
#################

results <- tryCatch(simpact.run(configParams = cfg.list,
                                destDir = sub.dir.rename,
                                agedist = age.distr,
                                seed = seedid,
                                intervention = intervention),
                    error = simpact.errFunction)



datalist <- readthedata(results)



# Epidemic characteristics: prevalence -----------------


# Women

hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist,
                                             agegroup = c(15, 25),
                                             timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.25.30.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(25, 30),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.30.35.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(30, 35),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.35.40.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(35, 40),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.40.45.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(40, 45),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.45.50.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(45, 50),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()



# Men

hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist,
                                           agegroup = c(15, 25),
                                           timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.25.30.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(25, 30),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.30.35.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(30, 35),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.35.40.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(35, 40),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.40.45.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(40, 45),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.45.50.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(45, 50),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()





agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")

men.prev <- c(hiv.prev.lt25.men, hiv.prev.25.30.men, hiv.prev.30.35.men, hiv.prev.35.40.men, hiv.prev.40.45.men, hiv.prev.45.50.men)
women.prev <- c(hiv.prev.lt25.women, hiv.prev.25.30.women, hiv.prev.30.35.women, hiv.prev.35.40.women, hiv.prev.40.45.women, hiv.prev.45.50.women)


men.prev.dat <- data.frame(agegroup, men.prev)
men.prev.dat$Gender <- "Male"
names(men.prev.dat) <- c("age", "prev", "Gender")

women.prev.dat <- data.frame(agegroup, women.prev)
women.prev.dat$Gender <- "Female"
names(women.prev.dat) <- c("age", "prev", "Gender")

prev_data <- rbind(men.prev.dat, women.prev.dat)

plot.prev.men.women <- ggplot(prev_data, aes(x=age, y=prev, colour=Gender, group = Gender)) + 
  # geom_errorbar(aes(ymin=L, ymax=U), width=.1) +
  geom_line(size=1) +
  geom_point() + 
  xlab("Age Groups") + ylab("HIV Prevalence")


ggsave(filename = "prevalence.pdf",
       plot = plot.prev.men.women,
       path = "/home/david/SimpactCyanPhylo",
       width = 30, height = 20, units = "cm")




###########################################
# Step 2: Construct transmission networks #
###########################################


# Produce a list of transmission networks in epi object format

simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)




###############################
# Step 3: Sequence simulation #
###############################

dirseqgen <- "/home/david/SimpactCyanPhylo"
dirfasttree <- "/home/david/SimpactCyanPhylo" 
sub.dir.rename <- "/home/david/SimpactCyanPhylo"

sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                               sub.dir.rename = sub.dir.rename,
                               simpact.trans.net = simpact.trans.net, # simpact.trans.net,
                               seq.gen.tool = "seq-gen",
                               seeds.num = 777,
                               endpoint = 40,
                               limitTransmEvents = 7, # no less than 7
                               hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                               clust = FALSE) # hiv.seq.file lodged in work.dir

# Output sequence file: C.Epidemic_seed.seq.bis.sim.nwk.fasta

# Transform the sequence format to be handled by ClusterPicker
sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")




#####################################################
# Step 4: Construct time stamped phylogenetic trees #
#####################################################

tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                             sub.dir.rename = sub.dir.rename,
                                             fasttree.tool = "FastTree",
                                             calendar.dates = "samplingtimes.all.csv",
                                             simseqfile = "C.Epidemic.fas",
                                             count.start = 1977,
                                             endsim = 40,
                                             clust = FALSE)

save(tree.calib, file = "tree.calib.RData")

# write.tree(tree.calib, file = paste0(sub.dir.rename, paste0("/calibrated.tree.",simseqfile,".tree")))
# calibrated.tree.C.Epidemic.fas.tree


# Node age with picante package

N <- node.age(tree.calib)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


# Molecular clock parameters

# Estimate of confidence intervals using parametric bootstrap for molecular clock dating.
# This function simulates phylogenies with branch lengths in units of substitutions per site. 
# Simulations are based on a fitted treedater object which provides parameters of 
# the molecular clock model. The treedater method is applied to each simulated tree 
# providing a Monte Carlo estimate of variance in rates and dates.

pbtd <- treedater::parboot(tree.calib) 

save(pbtd, file = "pbtd.RData")



# pseudo                         ML        2.5 %       97.5 %
# Time of common ancestor 1.987779e+03 1.987288e+03 1.988417e+03
# Mean substitution rate  5.550892e-03 5.459126e-03 5.644201e-03

# Confidence interval for time of common ancestor
# > (pbtd$timeOfMRCA_CI)
#       2.5%      97.5% 
#   1987.288   1988.417 



#####################################################################################
# Step 5: Compute transmission events and internal nodes in one year time intervals #
#####################################################################################


# THis help to assess transmission events and internal nodes distribution in one year time intervals

# We choose  transmission network 11 which is IDs.transm[1]


## Entire transmission network 

trans.net <- simpact.trans.net

trans.net <- data.table::rbindlist(simpact.trans.net)

trans.net <- dplyr::filter(trans.net, trans.net$parent != "-1") # remove universal seed individuals

# We considered the epidemic to start at 10 and individuals infected up to 40 simulation time: 30 years of epidemic
# With a seed sequence sampled in 1989, we assume it existed two years before (1987)
# It means that the simulation started in 1977, and the infection in 1987 for 30 years in 2017

trans.net$dtimes <- abs(trans.net$dtimes-40)+1977
trans.net$itimes <- abs(trans.net$itimes-40)+1977

min.val = 1977
max.val = round(max(trans.net$itimes))



# (i) NUmber of internal nodes & transmission events

step.int=1 # in one year

d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(trans.net)
dt.node.age.dt <- int.node.age

numC.tra <- vector() # initialize transmission events
i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  > inf),]
  numC.i <- nrow(dat.f.trans.i)
  numC.tra <- c(numC.tra, numC.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}


# (ii) Transmission network

graph.build <- as.data.frame(trans.net)

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
gag = gag[-1,] # remove universall seed -1
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transNet.yrs.Old <- ga.graph




#########################################################################
# Step 6: Creating the panels A, B and C of the figure for this example #
#########################################################################


# Objects needed for the figure:
SimpactPaperPhyloExample <- list()
SimpactPaperPhyloExample$transNet.yrs.Ord <- transNet.yrs.Old
SimpactPaperPhyloExample$dater.tree <- tree.calib
SimpactPaperPhyloExample$i.vec <- i.vec
SimpactPaperPhyloExample$int.node.vec <- int.node.vec
SimpactPaperPhyloExample$numC.tra <- numC.tra
SimpactPaperPhyloExample$pbtd <- pbtd

save(SimpactPaperPhyloExample, file = "SimpactPaperPhyloExample.RData")

SimpactPaperPhyloExample <- get(load(file = "SimpactPaperPhyloExample.RData"))


# A. Transmission network

network <- SimpactPaperPhyloExample$transNet.yrs.Ord
edges <- as.data.frame(get.edgelist(network))
names(edges) <- c("infector", "infectee")
vertices <- as.data.frame(vertex_attr(network))
network.df <- fortify(as.edgedf(edges), vertices)


transmissionnetwork.plot <- ggplot(data = network.df,
                                   aes(from_id = from_id,
                                       to_id = to_id)) +
  geom_net(directed = TRUE,
           size = 1.5,
           layout.alg = "kamadakawai",
           layout.par = list(niter = 2000,
                              sigma = 200,
                              kkconst = 10),
           alpha = 1,
           arrowsize = 0.7,
           arrowgap = 0.0015,
           ecolour = "darkgrey",
           colour = "black") +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  ylab("")

print(transmissionnetwork.plot)

ggsave(filename = "network.pdf",
       plot = transmissionnetwork.plot,
       path = "/home/david/SimpactCyanPhylo",
       width = 30, height = 20, units = "cm")



# C. Phylogenetic tree

tree <- SimpactPaperPhyloExample$dater.tree
class(tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
sim.start.year <- 1987
first.transmission <- min(SimpactPaperPhyloExample$dater.tree$sts)
mrsd <- max(SimpactPaperPhyloExample$dater.tree$sts)

dates <- format(date_decimal(c(mrsd, first.transmission)), "%Y-%m-%d")
tree$root.edge <-  - sim.start.year
phylotree.plot <- ggtree(tree,
                         mrsd = dates[1],
                         size = 0.1) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time") +
  ylab("")
print(phylotree.plot)

ggsave(filename = "phylotree.pdf",
       plot = phylotree.plot,
       path = "/home/david/SimpactCyanPhylo",
       width = 30, height = 20, units = "cm")




# C. Transmission event versus internal nodes

calendaryear.isrelevant <- SimpactPaperPhyloExample$i.vec >= 1987
calendaryear <- SimpactPaperPhyloExample$i.vec[calendaryear.isrelevant]
intern.nodes <- SimpactPaperPhyloExample$int.node.vec[calendaryear.isrelevant]
trans.events <- SimpactPaperPhyloExample$numC.tra[calendaryear.isrelevant]

trans.and.nodes.df <- data.frame(calendaryear = calendaryear,
                                 intern.nodes = intern.nodes,
                                 trans.events = trans.events)
trans.and.nodes.long.df <- gather(trans.and.nodes.df,
                                  key = "Events",
                                  value = "Count",
                                  intern.nodes:trans.events,
                                  factor_key = TRUE)


transandnodes.plot <- ggplot(data = trans.and.nodes.long.df,
                             aes(x = calendaryear,
                                 y = Count,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Internal nodes",
                                "Transmission events")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time")

print(transandnodes.plot)

ggsave(filename = "events_nodes.pdf",
       plot = transandnodes.plot,
       path = "/home/david/SimpactCyanPhylo",
       width = 30, height = 20, units = "cm")





# Visualize LTT (molecular clock)

pbtd <- SimpactPaperPhyloExample$pbtd




plot(pbtd, show.tree=TRUE, log.lineages=FALSE, main="Temporal trend of viral lineages through time")


ggsave(filename = "LTT.pdf",
       plot = plot(pbtd, show.tree=TRUE, log.lineages=FALSE, main="Temporal trend of viral lineages through time"),
       path = "/home/david/SimpactCyanPhylo",
       width = 30, height = 20, units = "cm")



