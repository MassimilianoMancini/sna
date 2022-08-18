################################################################################
###
### --- RscriptMultipleGroups.R: a script for the introduction to RSiena --- ###
###
###                         version March 7, 2022
################################################################################
#
# This is a script for demonstrating the multiple groups
# and meta-analysis options in RSiena.
# Written by Tom Snijders.
#
########################### ESTIMATION OF PARAMETERS ###########################

        library(RSiena)

# Download the data from http://www.stats.ox.ac.uk/~snijders/siena/CB_data.zip
# and unzip in your working directory.
# Read data for groups 1, 3, 4, 6:
# friendship waves 1 and 2, monadic attributes, dyadic attribute.
# Monadic: gender, delinquent behaviour, importance school friends.
# Dyadic: same ethnicity.

        CB01.w1 <- as.matrix(read.table("N34_1.DAT"))
        CB01.w2 <- as.matrix(read.table("HN34_1.DAT"))
        CB01.m  <- as.matrix(read.table("CBE1.DAT"))
        CB01.d  <- as.matrix(read.table("cbe1.sim"))

        CB03.w1 <- as.matrix(read.table("N34_3.DAT"))
        CB03.w2 <- as.matrix(read.table("HN34_3.DAT"))
        CB03.m  <- as.matrix(read.table("CBE3.DAT"))
        CB03.d  <- as.matrix(read.table("cbe3.sim"))

        CB04.w1 <- as.matrix(read.table("N34_4.DAT"))
        CB04.w2 <- as.matrix(read.table("HN34_4.DAT"))
        CB04.m  <- as.matrix(read.table("CBE4.DAT"))
        CB04.d  <- as.matrix(read.table("cbe4.sim"))

        CB06.w1 <- as.matrix(read.table("N34_6.DAT"))
        CB06.w2 <- as.matrix(read.table("HN34_6.DAT"))
        CB06.m  <- as.matrix(read.table("CBE6.DAT"))
        CB06.d  <- as.matrix(read.table("cbe6.sim"))

# Define Siena variables:
# Friendship
        F1  <- sienaNet(array( c(CB01.w1, CB01.w2), dim = c(45, 45, 2)))
        F3  <- sienaNet(array( c(CB03.w1, CB03.w2), dim = c(37, 37, 2)))
        F4  <- sienaNet(array( c(CB04.w1, CB04.w2), dim = c(33, 33, 2)))
        F6  <- sienaNet(array( c(CB06.w1, CB06.w2), dim = c(36, 36, 2)))

# sex
        sex1 <- coCovar(CB01.m[ ,1])
        sex3 <- coCovar(CB03.m[ ,1])
        sex4 <- coCovar(CB04.m[ ,1])
        sex6 <- coCovar(CB06.m[ ,1])

#delinquency
        del1 <- coCovar(CB01.m[ ,2])
        del3 <- coCovar(CB03.m[ ,2])
        del4 <- coCovar(CB04.m[ ,2])
        del6 <- coCovar(CB06.m[ ,2])

# importance school friends
        impsf1 <- coCovar(CB01.m[ ,3])
        impsf3 <- coCovar(CB03.m[ ,3])
        impsf4 <- coCovar(CB04.m[ ,3])
        impsf6 <- coCovar(CB06.m[ ,3])

# coethnic
        coethn1 <- coDyadCovar(CB01.d)
        coethn3 <- coDyadCovar(CB03.d)
        coethn4 <- coDyadCovar(CB04.d)
        coethn6 <- coDyadCovar(CB06.d)

# Create the four data sets, giving the variables
# the same names.
# They must have the same names to be combined later on.

        dataset.1  <- sienaDataCreate(F = F1,
                                      sex = sex1, del = del1, impsf = impsf1,
                                      coethn = coethn1)
        dataset.3  <- sienaDataCreate(F = F3,
                                      sex = sex3, del = del3, impsf = impsf3,
                                      coethn = coethn3)
        dataset.4  <- sienaDataCreate(F = F4,
                                      sex = sex4, del = del4, impsf = impsf4,
                                      coethn = coethn4)
        dataset.6  <- sienaDataCreate(F = F6,
                                      sex = sex6, del = del6, impsf = impsf6,
                                      coethn = coethn6)

# Put them together
        FourGroups <- sienaGroupCreate(list(dataset.1, dataset.3, dataset.4, dataset.6))
        GroupEffects <- getEffects(FourGroups)

# Get the initial description
        print01Report(FourGroups, modelname = 'CB1346')

# Extend the model
        GroupEffects <- includeEffects( GroupEffects, transTrip, transRecTrip )
        GroupEffects <- includeEffects( GroupEffects, simX, interaction1 = "sex" )
# Inspect current model specification
        GroupEffects

# Run
        GroupsModel <- sienaAlgorithmCreate(projname = 'CB1346')
        model.1 <- siena07(GroupsModel, data = FourGroups, effects = GroupEffects)
        model.1
		
		
########################### CHECKING TIME HOMOGENEITY  #########################

# Check homogeneity
        timetest.1 <- sienaTimeTest(model.1)
        summary(timetest.1)

# The second group, dataset.3, seems to be the villain.
# Try with only datasets 1, 4, 6.

        ThreeGroups <- sienaGroupCreate(list(dataset.1, dataset.4, dataset.6))
        ThreeGroupEffects <- getEffects(ThreeGroups)

# Get the initial description
        print01Report(ThreeGroups,  modelname = 'CB146')

# Extend the model
        ThreeGroupEffects <- includeEffects( ThreeGroupEffects,
                                             transTrip, transRecTrip )
        ThreeGroupEffects <- includeEffects( ThreeGroupEffects,
                                             simX, interaction1 = "sex" )
# Inspect current model specification
        ThreeGroupEffects

# Run
        ThreeGroupsModel <- sienaAlgorithmCreate(projname = 'CB146')
        (model3.1 <- siena07(ThreeGroupsModel, data = ThreeGroups,
                            effects = ThreeGroupEffects))

# Check homogeneity again
        timetest3.1 <- sienaTimeTest(model3.1)
        summary(timetest3.1)

# This seems somewhat better; mainly differences in density parameter.
        ThreeGroupEffects <- includeTimeDummy(ThreeGroupEffects, density, timeDummy="all")
        model3.2 <- siena07(ThreeGroupsModel, data = ThreeGroups,
                            effects = ThreeGroupEffects)

# Another approach, perhaps giving more insight,
# is to define dummy variables explicitly.
        dum3.1  <- c(rep(0,length(dataset.1$nodeSets$Actors)))
        dum3.3  <- c(rep(1,length(dataset.3$nodeSets$Actors)))
        dum3.4  <- c(rep(0,length(dataset.4$nodeSets$Actors)))
        dum3.6  <- c(rep(0,length(dataset.6$nodeSets$Actors)))
        dum4.1  <- c(rep(0,length(dataset.1$nodeSets$Actors)))
        dum4.3  <- c(rep(0,length(dataset.3$nodeSets$Actors)))
        dum4.4  <- c(rep(1,length(dataset.4$nodeSets$Actors)))
        dum4.6  <- c(rep(0,length(dataset.6$nodeSets$Actors)))
        dum6.1  <- c(rep(0,length(dataset.1$nodeSets$Actors)))
        dum6.3  <- c(rep(0,length(dataset.3$nodeSets$Actors)))
        dum6.4  <- c(rep(0,length(dataset.4$nodeSets$Actors)))
        dum6.6  <- c(rep(1,length(dataset.6$nodeSets$Actors)))
        dumm4.1 <- coCovar(dum4.1, warn = FALSE)
        dumm4.4 <- coCovar(dum4.4, warn = FALSE)
        dumm4.6 <- coCovar(dum4.6, warn = FALSE)
        dumm6.1 <- coCovar(dum6.1, warn = FALSE)
        dumm6.4 <- coCovar(dum6.4, warn = FALSE)
        dumm6.6 <- coCovar(dum6.6, warn = FALSE)

        datasetd.1  <- sienaDataCreate(F = F1,
                                      sex = sex1, del = del1, impsf = impsf1,
                                      dum4 = dumm4.1, dum6 = dumm6.1,
                                      coethn = coethn1)
        datasetd.4  <- sienaDataCreate(F = F4,
                                      sex = sex4, del = del4, impsf = impsf4,
                                      dum4 = dumm4.4, dum6 = dumm6.4,
                                      coethn = coethn4)
        datasetd.6  <- sienaDataCreate(F = F6,
                                      sex = sex6, del = del6, impsf = impsf6,
                                      dum4 = dumm4.6, dum6 = dumm6.6,
                                      coethn = coethn6)

        Threed.Groups <- sienaGroupCreate(list(datasetd.1, datasetd.4, datasetd.6))
        Threed.GroupEffects <- getEffects(Threed.Groups)

# Get the initial description
        print01Report(Threed.Groups, modelname = 'CBd146')

# Extend the model
        Threed.GroupEffects <- includeEffects( Threed.GroupEffects,
                                             transTrip, recTransTrip )
        Threed.GroupEffects <- includeEffects( Threed.GroupEffects,
                                             simX, interaction1 = "sex" )

# Since there are three groups in total, we need to include two group dummies
# (the number of degree of freedom at the group level
#  is the number of groups minus 1).

        Threed.GroupEffects <- includeEffects( Threed.GroupEffects,
                                             egoX, interaction1 = "dum4" )
        Threed.GroupEffects <- includeEffects( Threed.GroupEffects,
                                             egoX, interaction1 = "dum6" )
# Inspect current model specification
        Threed.GroupEffects

# Run
        Threed.GroupsModel <- sienaAlgorithmCreate(projname = 'CB146')
        model3d.1 <- siena07(Threed.GroupsModel, data = Threed.Groups,
                            effects = Threed.GroupEffects)
        model3d.1

# Check homogeneity once more.
# We need to request homogeneity tests only for effects that are not
# interacted with time dummies; else there will be an error.

        timetest3d.1 <- sienaTimeTest(model3d.1, effects=c(2:5))
        summary(timetest3d.1)

# This is satisfactory from an overall point of view.
# However, the parameters for transitivity effect, when considered by itself,
# are significantly heterogeneous.
# As an illustration, we specify a model where transitivity is heterogeneous.

        Threed.GroupEffects <- includeInteraction( Threed.GroupEffects,
                                  transTrip, egoX, interaction1 = c("","dum4") )
        Threed.GroupEffects <- includeInteraction( Threed.GroupEffects,
                                  transTrip, egoX, interaction1 = c("","dum6") )
        Threed.GroupEffects

# Run again
        (model3d.2 <- siena07(Threed.GroupsModel, data = Threed.Groups,
                            effects = Threed.GroupEffects))

# Check homogeneity once more
        timetest3d.2 <- sienaTimeTest(model3d.2, effects=c(2,4:5))
        summary(timetest3d.2)
# All good.


################# ALTERNATIVE: META-ANALYSIS  ##################################

################################################################################
# Now let us illustrate the meta analysis, i.e.,
# the analysis for each group separately, followed by a combination (siena08).
# This illustration is for four groups only, just to demonstrate the syntax,
# but the random effects procedure of Snijders & Baerveldt (2003)
# is not suitable for a small number of groups
# (minimum number will be somewhere between 10 and 20)
# because of the assumption that the groups are a sample from a population.
# The Fisher combination of tests, and the heterogeneity test,
# and the overall test, however, are suitable for any number of groups (even 2).
################################################################################

# We first need all four projects.

        effects.1 <- getEffects(dataset.1)
        effects.3 <- getEffects(dataset.3)
        effects.4 <- getEffects(dataset.4)
        effects.6 <- getEffects(dataset.6)

        print01Report(dataset.1, modelname = 'CB1')
        print01Report(dataset.3, modelname = 'CB3')
        print01Report(dataset.4, modelname = 'CB4')
        print01Report(dataset.6, modelname = 'CB6')

        effects.1 <- includeEffects( effects.1, transTrip, transRecTrip)
        effects.1 <- includeEffects( effects.1, sameX, interaction1 = "sex" )
        effects.3 <- includeEffects( effects.3, transTrip, transRecTrip )
        effects.3 <- includeEffects( effects.3, sameX, interaction1 = "sex" )
        effects.4 <- includeEffects( effects.4, transTrip, transRecTrip)
        effects.4 <- includeEffects( effects.4, sameX, interaction1 = "sex" )
        effects.6 <- includeEffects( effects.6, transTrip, transRecTrip )
        effects.6 <- includeEffects( effects.6, sameX, interaction1 = "sex" )

# Then we need to estimate all of them.

        algo.1 <- sienaAlgorithmCreate(projname = 'CB1')
        algo.3 <- sienaAlgorithmCreate(projname = 'CB3')
        algo.4 <- sienaAlgorithmCreate(projname = 'CB4')
        algo.6 <- sienaAlgorithmCreate(projname = 'CB6')
        (model1.1 <- siena07(algo.1, data = dataset.1, effects = effects.1))
        (model3.1 <- siena07(algo.3, data = dataset.3, effects = effects.3))
        (model4.1 <- siena07(algo.4, data = dataset.4, effects = effects.4))
        (model6.1 <- siena07(algo.6, data = dataset.6, effects = effects.6))

# Combine the four projects and do the meta-analysis
        meta.1346 <- siena08(model1.1, model3.1, model4.1, model6.1,
                             projname = "meta1346")
# There is so much output that it may be more convenient to write it
# to a separate output file, as follows:
        sink("meta1346.out")
        meta.1346
        sink()

# Look at
?siena08
# to see the elements produced by siena08.
# All these are reported in the output sent above to file meta1346.out.
# The various elements produced can also be made available in the R session;
# the following gives examples.

# If you are a LaTeX user, you can use the output of
meta.table(meta.1346)
# Look at the different options in
?meta.table

# The number of effects (leaving out the rate effect used for conditioning) in
effects.1
# is
(neff <- sum(effects.1$type[effects.1$include] != "rate"))

# Parameter estimates and confidence intervals:
parameters.1346 <- t(sapply(1:neff, function(i){c(meta.1346[[i]]$mu.ml,
                    meta.1346[[i]]$mu.ml.se, meta.1346[[i]]$mu.confint,
                    meta.1346[[i]]$sigma.ml, meta.1346[[i]]$sigma.confint,
                    meta.1346[[i]]$n1)}))

colnames(parameters.1346) <- c('mu-hat', 'mu-se',
                    'mu-min', 'mu-plus', 'alpha_mu',
                    'sigma-hat', 'sigma-min', 'sigma-plus', 'alpha_sigma', 'N')
# These are: parameter estimate for population mean; s.e.;
# confidence interval for population mean (left, right, significance level);
# estimate for population standard deviation;
# confidence interval for standard deviation (left, right, significance level);
# number of groups on which this is based.

# Construct the names for the effects:
efnames <- names(meta.1346)[1:neff]
# the first 5 list elements of meta.1346 correspond to the tested effects
# The first 7 letters are redundant and will be dropped:
efnames <- substring(efnames, 8)
rownames(parameters.1346) <- efnames
round(parameters.1346, 3)

# Extract heterogeneity tests from siena08 output;
# these are called Q in Snijders & Baerveldt (2003):
hetero.1346 <- t(sapply(1:neff, function(i){c(meta.1346[[i]]$Qstat,
                        meta.1346[[i]]$n1-1, meta.1346[[i]]$pttilde)}))
colnames(hetero.1346) <- c('Q', 'df', 'pQ')
# These are: Q statistic; degrees of freedom; p-value
rownames(hetero.1346) <- efnames
round(hetero.1346, 3)

# Extract overall tests;
Overalls.1346 <- t(sapply(1:neff, function(i){c(meta.1346[[i]]$Tsq,
                        meta.1346[[i]]$n1-1, meta.1346[[i]]$pTsq)}))
colnames(Overalls.1346) <- c('T^2', 'df', 'pT^2')
# These are: Q statistic; degrees of freedom; p-value
rownames(Overalls.1346) <- efnames
round(Overalls.1346, 3)

# We also present the Fisher combinations:
Fishers <- t(sapply(1:neff,
        function(i){c(meta.1346[[i]]$cjplus, meta.1346[[i]]$cjminus,
                        meta.1346[[i]]$cjplusp, meta.1346[[i]]$cjminusp,
                        2*meta.1346[[i]]$n1 )}))
rownames(Fishers) <- efnames
colnames(Fishers) <- c('Fplus', 'Fminus', 'pplus', 'pminus', 'df')
round(Fishers,3)

# If you wish to combine only the three projects:
        meta.146 <- siena08(model1.1, model4.1, model6.1, projname = "meta146")
# and proceed as above.

# The multi-group project has much smaller standard errors
# than the meta analysis.
# This is because the assumptions are different;
# the multi-group project assumes that the parameters for which no dummy
# variables are used to differentiate them between projects,
# are the same across projects.
# The standard errors for the meta analysis, on the other hand,
# take into account not only the standard error of each estimate,
# but also the between-groups variation of the estimates.

# Visual exploration is possible with
?funnelPlot
# The funnel plot for same sex:
funnelPlot(list(model1.1, model3.1, model4.1, model6.1), k=5)

################# GOODNESS OF FIT FOR SIENAGROUP RESULTS #######################

# For checking goodness of fit, function sienaGOF can be used for each group,
# and then the results can be combined.
# This is illustrated here for model3d.1.
# First it is estimated again, to obtain high precision.
# See Section 6.4 of the manual.

        Threed.GroupsModel.x <- sienaAlgorithmCreate(projname = 'CB146',
                                    nsub=1, n2start=2000, n3=2000)
        (model3d.3 <- siena07(Threed.GroupsModel.x, data = Threed.Groups,
                            effects = Threed.GroupEffects, prevAns=model3d.1))

# Now these results are run with a shorter phase 3 
# and returnDeps=TRUE for sienaGOF:

        Threed.GroupsModel.d <- sienaAlgorithmCreate(projname = 'CB146',
                                    nsub=0, n3=1000)
        (model3d.4 <- siena07(Threed.GroupsModel.d, data = Threed.Groups,
                            effects = Threed.GroupEffects, prevAns=model3d.3,
                            returnDeps=TRUE))

# Apply sienaGOF for the four main auxiliary functions for networks.
?sienaGOF
# This requires using the 'groupName' parameter. What are the groupNames?
names(Threed.Groups)

# One way to obtain sienaGOF results for the outdegree distribution
# for each group is
gof.o <- list()
for (k in 1:length(Threed.Groups)){
    gof.o[[k]] <- sienaGOF(model3d.4,
                        OutdegreeDistribution, verbose=TRUE, join=TRUE,
                        varName="F", groupName=names(Threed.Groups)[k])
}

# If you know how to work with lists and lapply, a more efficient way is:

gof.o <- lapply(1:length(Threed.Groups), function(k){sienaGOF(model3d.4,
                        OutdegreeDistribution, verbose=TRUE, join=TRUE,
                        varName="F", groupName=names(Threed.Groups)[k])})

# This is done also for the other three auxiliary statistics:

gof.i <- lapply(1:length(Threed.Groups), function(k){sienaGOF(model3d.4,
                        IndegreeDistribution, verbose=TRUE, join=TRUE,
                        varName="F", groupName=names(Threed.Groups)[k])})

gof.tc <- lapply(1:length(Threed.Groups), function(k){sienaGOF(model3d.4,
                        TriadCensus, verbose=TRUE, join=TRUE,
                        varName="F", groupName=names(Threed.Groups)[k])})

# For the geodesic distribution, we need to take the auxiliary function from
?"sienaGOF-auxiliary"
# (note the quotes)

GeodesicDistribution <- function (i, data, sims, period, groupName,
  varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(network)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative)
  {
    gdi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
    gdi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}
# if necessary:
# install.packages(c("network", "sna"))

gof.gd <- lapply(1:length(Threed.Groups), function(k){sienaGOF(model3d.4,
                        GeodesicDistribution, verbose=TRUE, join=TRUE,
                        varName="F", groupName=names(Threed.Groups)[k])})

# Inspect the results:
sapply(gof.o, function(x){x[[1]]$p})
sapply(gof.i, function(x){x[[1]]$p})
sapply(gof.tc, function(x){x[[1]]$p})
sapply(gof.gd, function(x){x[[1]]$p})

# To get a single p=value combining the groups,
# Fisher's combination of p-values may be used
# (see the RSiena manual, section 5.13.2):

(p.o <- sum(sapply(gof.o, function(x){-2*log(x[[1]]$p)})))
pchisq(p.o, df = 6, lower.tail=FALSE)
(p.i <- sum(sapply(gof.i, function(x){-2*log(x[[1]]$p)})))
pchisq(p.i, df = 6, lower.tail=FALSE)
(p.tc <- sum(sapply(gof.tc, function(x){-2*log(x[[1]]$p)})))
pchisq(p.tc, df = 6, lower.tail=FALSE)
(p.gd <- sum(sapply(gof.gd, function(x){-2*log(x[[1]]$p)})))
pchisq(p.gd, df = 6, lower.tail=FALSE)

# For the degree distributions, the fit is good.
# However, for the other two the fit is inadequate.
# Depending on the random number seed, you may get a straight 0
# for the p-value of the triad census fit, which will lead to Inf for the log.

# It is best to start improving the fit for the triad census,
# and hope that this also will improve the fit for the geodesic distribution.
# The start of this process is to make plots

plot(gof.tc[[1]], center=TRUE, scale=TRUE)
plot(gof.tc[[2]], center=TRUE, scale=TRUE)
plot(gof.tc[[3]], center=TRUE, scale=TRUE)

# The continuation is an exercise.



