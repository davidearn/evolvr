#' Evolve each of a sequence of models to fixation
#'
#' @param simParms a \code{\link[=set_simParms]{simParms}} object
#' @inheritParams dbenefit
#' @export
multiple_evolve_to_fixation <- function(model, simParms, ...) {

    stopifnot("snowdrift.game" %in% class(model))
    stopifnot("simParms" %in% class(simParms))
    s <- simParms

    ## initialize time for this sequence of simulations:
    pt0 <- proc.time()

    ## column names of output matrix:
    outNames <- c("iRes", "iMut", "resident", "mutant",
                  "mutantFixationProbability",
                  "mutantFixationTimeMean",
                  "mutantFixationTimeSD",
                  "mutantFixationTimeMin",
                  "mutantFixationTimeMax",
                  "residentFixationTimeMean",
                  "residentFixationTimeSD",
                  "residentFixationTimeMin",
                  "residentFixationTimeMax"
                  )

    ## numbers of resident and mutant strategies:
    nRes <- s$nResidentStrategies
    nMut <- s$nMutantStrategies

    ## output matrix:
    out <- matrix(data = NA, ncol = length(outNames),
                  nrow = nRes * nMut,
                  dimnames = list(NULL,outNames))
    ## add class:
    class(out) <- c("fixationData", class(out))

    ## save information as attributes of output matrix:
    attr(out,"model") <- model
    attr(out,"simParms") <- simParms
    out <- add_system_attributes( out )

    iRow <- 0
    for (iRes in seq_along(s$residentStrategies)) {
        resident <- s$residentStrategies[iRes]
        cat("resident[", iRes, "] = ", resident, "\n", sep="")
        for (iMut in seq_along(s$mutantStrategies)) {
            iRow <- iRow + 1
            cat(iMut, " ", sep="")
            mutant <- s$mutantStrategies[iMut]
            ## reset for this mutant:
            pop0 <- initialize_population(model, resident, mutant, nMut=1)
            mutFixCount <- 0
            didNotFixCount <- 0
            mutFixTime <- rep(NA,s$nReplicates) # mutant fixation times
            resFixTime <- rep(NA,s$nReplicates) # resident fixation times
            ## run all replicates:
            for (iSim in 1:s$nReplicates) {
                fixationRecord <-
                    evolve_to_fixation(model, pop0, s$maxGen, s$selectionProcess)
                if (is.na(fixationRecord$mutantFixed)) {
                    didNotFixCount <- didNotFixCount + 1
                } else if (fixationRecord$mutantFixed) {
                    mutFixCount <- mutFixCount + 1
                    mutFixTime[iSim] <- fixationRecord$fixTime
                } else {
                    resFixTime[iSim] <- fixationRecord$fixTime
                }
                ##resFixCount <- s$nReplicates - mutFixCount - didNotFixCount
            } # end iSim loop
            if (didNotFixCount > 0)
                warning(didNotFixCount, " of ", s$nReplicates,
                        " replicates did not fix")
            out[iRow,] <- compute_fixation_stats(
                out[iRow,], s$nReplicates,
                mutFixCount, didNotFixCount,
                mutFixTime, resFixTime)
            ## FIX: the following should work and should be in compute_fixation_stats
            ##out[iRow,1:4] <- c(iRes, iMut, resident, mutant)
            out[iRow,1] <- iRes
            out[iRow,2] <- iMut
            out[iRow,3] <- resident
            out[iRow,4] <- mutant
            pt <- summary(proc.time() - pt0)
            attr(out,"procTime") <- pt # update cpu time
            save_simulation_progress(out, iRes, iMut)
        } # end iMut loop
        cat("\n")
    } # end iRes loop
    print(pt)
    return(out)
}

#' Save simulation progress
#'
#' This function is called by
#' \code{\link{multiple_evolve_to_fixation}} to save data periodically
#' during a sequence of simulations, and to write a plain text
#' progress file to help monitor a job's status.
#'
#' @param x a \code{fixationData} object
#' @param iRes index of resident strategy
#' @param iMut index of mutant strategy
#' @note FIX: this function should probably return a status flag,
#'     indicating whether saving was successful.
#' @export
save_simulation_progress <- function( x, iRes=NA, iMut=NA ) {
    stopifnot("fixationData" %in% class(x))
    s <- attr(x, "simParms")
    pid <- attr(x, "pid")
    progressFile <- paste0(pid, ".progress")
    pt <- attr(x,"procTime")
    ptCPU <- pt["user"] + pt["system"]
    cpuTime <- time_string(ptCPU)
    ptElapsed <- pt["elapsed"]
    elapsedTime <- time_string(ptElapsed)
    nodeName <- attr(x,"sysInfo")["nodename"]
    out <- data.frame(nodeName=nodeName, pid=pid,
                      nRes = s$nResidentStrategies,
                      iRes = iRes,
                      nMut = s$nMutantStrategies,
                      iMut = iMut,
                      cpuTime=cpuTime, elapsedTime=elapsedTime)
    write.table(out,file=progressFile,quote=FALSE,row.names=FALSE,sep="\t")
    saveFile <- paste0(pid, ".RData")
    fixationData <- x
    save(iRes, iMut, fixationData, file=saveFile)
}


#' Compute fixation statistics
#'
#' @param x a single fixation record
#' @export
compute_fixation_stats <- function(x, nReplicates,
                                   mutFixCount, didNotFixCount,
                                   mutFixTime, resFixTime) {

    ## FIX: also save quantiles of mutFixTime and resFixTime ?
    ## and residentFixationProbability since this is not 1-mutFixProb
    ## if didNotFixCount > 0.

    x["mutantFixationProbability"] <- mutFixCount / nReplicates

    if (mutFixCount > 0) {
        x["mutantFixationTimeMean"] <- mean(mutFixTime, na.rm=TRUE)
        x["mutantFixationTimeSD"] <- sd(mutFixTime, na.rm=TRUE)
        x["mutantFixationTimeMin"] <- min(mutFixTime, na.rm=TRUE)
        x["mutantFixationTimeMax"] <- max(mutFixTime, na.rm=TRUE)
    } else {
        x["mutantFixationTimeMean"] <- NA
        x["mutantFixationTimeSD"] <- NA
        x["mutantFixationTimeMin"] <- NA
        x["mutantFixationTimeMax"] <- NA
    }

    resFixCount <- nReplicates - mutFixCount - didNotFixCount

    if (resFixCount > 0) {
        x["residentFixationTimeMean"] <- mean(resFixTime, na.rm=TRUE)
        x["residentFixationTimeSD"] <- sd(resFixTime, na.rm=TRUE)
        x["residentFixationTimeMin"] <- min(resFixTime, na.rm=TRUE)
        x["residentFixationTimeMax"] <- max(resFixTime, na.rm=TRUE)
    } else {
        x["residentFixationTimeMean"] <- NA
        x["residentFixationTimeSD"] <- NA
        x["residentFixationTimeMin"] <- NA
        x["residentFixationTimeMax"] <- NA
    }

    return(x)
}

#' Create a population matrix object
#'
#' @inheritParams create_snowdrift_game
#' @return \code{\link{matrix}} with \code{Npop} rows and three
#'     columns: \code{"ID"}, \code{"strategy"}, and \code{"fitness"}
#' @seealso \code{\link{initialize_population}},
#'     \code{\link{initialize_population_randomly}}
#' @examples
#' population_matrix( 10 )
#' @export
#' @import sticky
population_matrix <- function( Npop ) {
    pop <- matrix(data = NA, ncol=3, nrow=Npop,
                  dimnames = list(NULL,c("ID", "strategy","fitness")))
    pop[,"ID"] <- 1:Npop
    pop<-sticky(pop) #causes attributes not to be discarded when operations are performed on pop

    return( pop )
}


#' Initialize population with residents and mutants
#'
#' Population is represented by an array of numbers, with rows
#' representing individual agents.  Columns are \code{"strategy"} and
#' \code{"fitness"}.  Mutants are randomly distributed among groups.
#' @note Agent in row \code{i} is considered to be in group
#'     \code{ceiling(i/nPlayers)}, where groups are numbered
#'     \code{1:nGroups} and rows are numbered \code{1:Npop}.
#' @inheritParams dbenefit
#' @param resident resident strategy
#' @param mutant mutant strategy
#' @param nMut number of mutants
#' @return a \code{\link{population_matrix}}
#' @seealso \code{\link{initialize_population_randomly}}
#' @export
initialize_population <- function(model, resident, mutant, nMut) {
    Npop <- model$parms$Npop
    pop <- population_matrix( Npop )

    ## rep.int is faster than rep; see ?rep
    pop[,"strategy"] <- rep.int(x=c(resident,mutant), times=c(Npop-nMut,nMut))

    #set initial fitnesses to 0 (should this be done in population_matrix?
    # this replicates code in initialize_population_randomly)
    pop[,"fitness"] <-rep.int(x=0,times=Npop)

    #commented the next two lines since population is now randomized in compute_fitnesses
    # ## permute population so mutants are randomly distributed among groups:
    # pop[,"strategy"] <- sample(pop[,"strategy"])

    #No longer necessary because fitness is now computed in the beginning of selection step
    # ## compute fitness of each agent:
    # pop <- compute_fitnesses(model, pop)

    ## save parameters that determined population structure as attributes:
    attr(pop,"resident") <- resident
    attr(pop,"mutant") <- mutant
    attr(pop,"nMut") <- nMut

    return(pop)
}

#' Initialize population via a Normal distribution of strategies
#'
#' @inheritParams dbenefit
#' @param xMean mean strategy
#' @param xSD standard deviation of strategy distribution
#' @param xMin,xMax extreme values of playable strategies (boundaries
#'     of strategy space)
#' @return a \code{\link{population_matrix}}
#' @seealso \code{\link{initialize_population}}, \code{\link{rnorm}}
#' @export
initialize_population_randomly <- function( model, xMean, xSD, mutRange=c(-Inf,Inf), xMin, xMax ) {
    Npop <- model$parms$Npop

    # pop[,"strategy"] <- rnorm(n = Npop, mean = xMean, sd = xSD)
    pop <- initialize_population(model, resident=xMean, mutant=xMean, nMut=0)
    # remove unnecessary attributes added in initialize_population:
    attr(pop,"resident") <- NULL
    attr(pop,"mutant") <- NULL
    attr(pop,"nMut") <- NULL

    whichAgentsMutate <- 1:Npop
    pop <- mutate_agent_strategies(model, pop, whichAgentsMutate, mutProb, mutSD, mutRange, xMin, xMax )

    # now taken care of in mutate_agent_strategies
    #pop[,"strategy"] <- ifelse(pop[,"strategy"] < xMin, xMin, pop[,"strategy"])
    #pop[,"strategy"] <- ifelse(pop[,"strategy"] > xMax, xMax, pop[,"strategy"])

    #now done in initialize_population
    #set initial fitnesses to 0 (should this be done in population_matrix?)
    #pop[,"fitness"] <-rep.int(x=0,times=Npop)


    ## save parameters that determined population structure as attributes:
    attr(pop,"xMean") <- xMean
    attr(pop,"xSD") <- xSD
    attr(pop,"xMin") <- xMin
    attr(pop,"xMax") <- xMax

    return(pop)
}

#' Evolve a population until fixation.
#'
#' @inheritParams dbenefit
#' @param pop population matrix (columns: \code{strategy} and
#'     \code{fitness})
#' @param maxGen maximum number of generations to evolve
#' @param selectionProcess selection process (\code{"WF"} = Wright-
#'     Fisher, \code{"M"} = Moran, \code{"BD"} = birth-death)
#' @return a list containing a logical flag (\code{mutantFixed}) and
#'     the number of generations that were computed (\code{fixTime})
#' @export
evolve_to_fixation <- function(model, pop, maxGen, selectionProcess) {

    ## choose next generation function
    selection_step <- switch(selectionProcess,
                       WF = selection_step_WF,
                       M  = selection_step_M,
                       BD = selection_step_BD)

    iGen <- 0 # generation counter
    ## pop0 <- pop # initial population # not used
    hasFixed <- FALSE
    while ( (!hasFixed) && (iGen < maxGen) ) {
        pop <- selection_step(model, pop)

        ## Flag indicating fixation
        hasFixed <- all_the_same(pop[,"strategy"])
        attr(pop,"hasFixed") <- hasFixed #Is this attribute used elsewhere?
        iGen <- iGen + 1
    }
    if (hasFixed) {
        ## all strategies are identical so we can check the first:
        mutantFixed <- (pop[1,"strategy"] == attr(pop,"mutant"))
    } else  {
        remainingStrategies <- unique(pop[,"strategy"])
        message(paste("evolve_to_fixation:\n ",
                      selectionProcess,
                      "selection process did not terminate in fixation after",
                      iGen, "generations.\n  The",
                      length(remainingStrategies),
                      "remaining strategies were:",
                      paste(signif(remainingStrategies,4), collapse=" ")
                      ))
        ## FIX: any need to return the strategies that remain?
        mutantFixed <- NA
    }
    return(list(mutantFixed = mutantFixed,
                fixTime = iGen))
}

#' Make fitness positive
#'
#' Probabilities of generating an offspring are supposed to be
#' proportional to fitness, but some fitnesses are negative. This
#' makes them positive.
#' @details It is important the \code{fitnessMin} be positive.  If it
#'     is zero, then mutants will never fix if they have unadjusted
#'     fitness less than residents.
#' @inheritParams selection_step_WF
#' @seealso \code{\link{selection_step_WF}}
#' @export
make_fitness_positive <- function( pop, fitnessMin ) {
    ##pop <- shift_fit_1(model, pop)
    ## DE: changing this to adjust by the current minimum, which will
    ## work for any fitness function:
    popOrig <- pop
    pop[,"fitness"] <- popOrig[,"fitness"] - min(popOrig[,"fitness"]) + fitnessMin
    return( pop )
}

#' Advance population by one generation using Wright-Fisher selection process
#'
#' @inheritParams dbenefit
#' @param pop population matrix: named columns (\code{fitness} and
#'     \code{strategy}) and one row for each agent
#' @details FIX: this is mostly old comments.
#' Assumes fitnesses have been calculated and are available
#'     as the second column of the \code{pop} matrix. Next generation
#'     is constructed using N independent draws from the pool of
#'     strategies currently in the population. The probability of
#'     drawing strategy x in the pool for any one of the next
#'     generation's agents is proportional to the total fitness from
#'     all the agents playing x in the current population (that is, to
#'     the mean fitness of x and to the number of agents playing
#'     x). Randomizes group allocation in population, and calculates
#'     individual fitnessess based on payoffs from games played in
#'     these groups.
#' @return population matrix with \code{hasFixed} attribute to
#'     indicate whether or not the population has reached fixation.
#' @export
## FIX: fitnessMin occurs only here: does the user need to be able to change it?
selection_step_WF <- function(model, pop, fitnessMin = 1) {

    Npop <- model$parms$Npop

    ## Note: cbind construction here produces an object of class matrix.
    ## 1st column is a list of the unique strategies in the population
    ## 2nd column will contain the mean fitness of each strategy
    uniqueStrategies <- cbind(unique(pop[,"strategy"]), 0)
    colnames(uniqueStrategies) <- c("strategy", "meanFitness")

    #commented the next two lines since population is now randomized in compute_fitnesses
    # ## permute the individual strategies (randomizes group assignment)
    # pop[,"strategy"] <- sample(pop[,"strategy"])
    pop <- compute_fitnesses(model, pop)

    pop <- make_fitness_positive( pop, fitnessMin )
    ## FIX: it would be more logical to compute_fitnesses immediately
    ##      before shifting the fitnesses, and not compute fitnesses
    ##      below, or to do the shifting below.  Why not do the shift
    ##      inside compute_fitnesses?  Is the shift necessary only for
    ##      WF?

    ## calculate total fitness for each strategy
    for (i in 1:nrow(uniqueStrategies)) {
        ## agents (indices) playing strategy s[i] :
        agents <- which(pop[,"strategy"] == uniqueStrategies[i])
        uniqueStrategies[i,"meanFitness"] <- sum(pop[agents,"fitness"])
    }

    ## Using rmultinom instead of rbinom to allow for multi-morphic
    ## (as opposed to di-morphic) populations. For anything beyond
    ## invasion analysis, we'll need the rmultinom, especially if/when
    ## mutations are possible.
    ##
    ## NOTE: probability is internally normalized by rmultinom.
    ##
    ## offspring = number of offspring for each strategy:
    offspring <- rmultinom(n=1, size = Npop, prob = uniqueStrategies[,"meanFitness"])
    pop[,"strategy"] <- rep(uniqueStrategies[,"strategy"], offspring)

    return(pop)
}

###---------------------------------------------------------------------------------------------------------------

#' Advance population by one generation using Moran selection process
#'
#' Constructs new generation based on frequency-dependent Moran
#' process
#' @inheritParams selection_step_WF
#' @details Assumes fitnesses are calculated. Agent is chosen to
#'     reproduce (with a probability proportional to its fitness) and
#'     another agent is chosen from the population (with uniform
#'     probability) to die. The strategy of the dying agent is
#'     replaced by the strategy of the reproducing one. Randomizes
#'     group allocation in the population, and calculates individual
#'     fitnessess based on payoffs from games played in these groups.
#' @return population matrix with \code{hasFixed} attribute to
#'     indicate whether or not the population has reached fixation.
#' @export
selection_step_M <- function(model, pop, fitnessMin = 1) {

    Npop <- model$parms$Npop

    #commented the next two lines since population is now randomized in compute_fitnesses
    # ## permute the individual strategies (randomizes group assignment)
    # pop[,"strategy"] <- sample(pop[,"strategy"])
    pop <- compute_fitnesses(model, pop)

    pop <- make_fitness_positive( pop, fitnessMin )

    #############################################################
    ###### FIX: not ideal that so much code is replicated  ######
    ###### THIS IS THE ONLY PART THAT IS DIFFERENT FROM WF ######
    reproducingAgent <- which(rmultinom(n=1, size=1, prob = pop[,"fitness"])==1)
    dyingAgent <- sample(Npop, 1)
    pop[dyingAgent,"strategy"] <- pop[reproducingAgent,"strategy"]
    #############################################################

    return(pop)
}
###---------------------------------------------------------------------------------------------------------------


#' Are all the components of an object identical?
#'
#' @details This is used in next generation functions to determine if
#'     an evolving population has fixed.
#' @seealso \code{\link{selection_step_WF}}, \code{\link{selection_step_M}},
#'     \code{\link{selection_step_BD}}
#' @param x an R object
#' @return logical
#' @export
all_the_same <- function( x ) {return(length(unique(x)) == 1)}

#' Advance by one generation of Birth-Death process
#'
#' Constructs new generation based on the transition probabilities of
#' a birth-death process, which are in turn calculated from the mean
#' fitnesses of mutants and residents in the current generation.  From
#' state i, the state can either increase to i+1, decrease to i-1, or
#' remain i.  This selection process is not defined for more than two
#' traits in the population.  Note that the fitnesses of the current
#' generation are not calculated here and are assumed to be input in
#' \code{pop[,"fitness"]}, but the fitnesses of the next generation
#' are output in \code{pop[,"fitness"]}.
#'
#' @inheritParams dbenefit
#' @param pop population matrix
#' @export
## FIX: see comments in selection_step_WF and make structure consistent
## 7 July 2019: added fitnessMin as arg -- was that what I meant to fix?
selection_step_BD <- function(model, pop, fitnessMin = 1) {
  nPlayers<- model$parms$nPlayers
  Npop <- model$parms$Npop

  ## if population is monomorphic then return immediately
  hasFixed <- all_the_same(pop[,"strategy"])
  attr(pop,"hasFixed") <- hasFixed
  if (hasFixed) return(pop)

  #commented the next two lines since population is now randomized in compute_fitnesses
  # ## permute the individual strategies (randomizes group assignment):
  # pop[,"strategy"] <- sample(pop[,"strategy"])
  pop <- compute_fitnesses(model, pop)

  # technically, we don't need to shift individual fitnesses,
  # only the mean fitnesses for each strategy
  pop <- make_fitness_positive( pop, fitnessMin )


  ## 1st col is a list of the unique strategies in the population
  ## 2nd col will contain the mean fitness of each strategy
  uniqueStrategies <- cbind(unique(pop[,"strategy"]),0)
  colnames(uniqueStrategies)<-c("strategy","meanFitness")

  ## calculate mean fitness for each strategy
  agents <- list()
  for (i in 1:dim(uniqueStrategies)[1]) {
    ## agents (indices) playing strategy uniqueStrategies[i]
    agents[[i]]<- which(pop[,"strategy"]==uniqueStrategies[i])
    uniqueStrategies[i,"meanFitness"]<-mean(pop[agents[[i]],"fitness"])
  }

  ## because the fitness difference is calculated as
  ## uniqueStrategies[2,"fitness"]-uniqueStrategies[1,"fitness"]
  ## transprob["up"] is the probability that the number of individuals
  ## playing uniqueStrategies[2,"strategy"] increases by one.
  transProbs <- compute_transition_prob_mat(fitDiff=diff(uniqueStrategies[,"meanFitness"]))
  ## codes for what happens to population state defined by number of
  ## mutants:  1="down", 2="stay", 3="up"
  bdEvent <- which(rmultinom(n=1, size=1, prob = transProbs)[,1]==1)
  if (bdEvent==1) {
    ## an individual playing uniqueStrategies[2,"strategy"]
    ## is replaced by an individual playing uniqueStrategies[1,"strategy"]
    pop[agents[[2]][1],"strategy"] <- uniqueStrategies[1,"strategy"]
  } else if (bdEvent==3) {
    ## an individual playing uniqueStrategies[1,"strategy"]
    ## is replaced by an individual playing uniqueStrategies[2,"strategy"]
    pop[agents[[1]][1],"strategy"] <- uniqueStrategies[2,"strategy"]
  } else {
    stopifnot(bdEvent==2)
    ##do nothing
  }
  return(pop) # 7 July 2019: how did it ever work without this???
}

#' Calculate transition probabilities for birth-death process defined
#' in invasion/fixation ms.
#'
#' If the fitness difference is calculated between type A and type B
#' (\code{W_A-W_B}), then the state of the population is the number of
#' individuals of type A in the population.
#' Output is a matrix of 3 columns and \code{(#states-1)} rows, labelled
#' \code{"down"}, \code{"up"} and \code{"stay"}.
#'     The i'th entry of the column labelled "up" is the transition probability from state i to state i+1.
#'     The i'th entry of the column labelled "down" is the transition probability from state i to state i-1.
#'     The i'th entry of the column labelled "stay" is the probability to remain at state i,
#'        which is 0 for this process.
#' The \code{"stay"} column is included for compatibility with the function NewGenerationBD.
#' This definition of the state of the population means that if \code{fitDiff>0}, then
#' \code{output["up"]>output["down"]}.
#'
#' @param fitDiff vector of fitness differences
#' @details This function can process a vector \code{fitDiff}.  For
#'     one population state, \code{fitDiff} is the mean fitness
#'     difference between mutants and residents, and thus a
#'     scalar. However, inputing the vector of fitness differences for
#'     all mixed states is useful (see
#'     \code{\link{MutantFixProbExactBD}}).
#' @note FIX: this function ought to take another argument, the
#'     transition probability ratio parameter \eqn{\phi}, with default
#'     \code{1}.  See Appendix A of invasion-fixation ms.  More
#'     generally, the implicit assumption that \eqn{F(x)=e^{-\phi x}}
#'     is unnecessary.  This function could take a function as
#'     argument, with \eqn{e^{-x}} as the default:
#'     \code{function(fitDiff, tpr_fun=function(x){exp(-x)})}
#'     
#' @export 
## FIX: unchanged from Chai's version, except name was calc_trans_probs_mat
compute_transition_prob_mat <- function(fitDiff) {
  ## Initialize matrix of transition probabilities with all entries 0.
  transProbs <- matrix(data = 0,ncol=3,nrow=length(fitDiff))

  colnames(transProbs) <- c("down", "stay", "up")

  transprobRatio <- exp(-fitDiff)

  transProbs[,"up"] <- 1/(1+transprobRatio)
  transProbs[,"down"] <- 1-transProbs[,"up"]
  return(transProbs)
}

#' Evolve a population under selection and mutation
#'
#' @inheritParams dbenefit
#' @inheritParams initialize_population_randomly
#' @param mutProb probability that a mutation occurs
#' @param mutSD standard deviation for normal distribution with zero
#'     mean, used to shift offspring strategy if a mutation occurs
#' @param mutRange the minimal and maximal magnitudes of mutations
#'     (=c(-Inf,Inf) by default)
#' @param xInitMean mean of a normally distributed set of initial
#'     strategies for the population
#' @param xInitSD standard deviation of normally distributed initial
#'     strategies with mean \code{xInitMean}
#' @param nGen number of generations to simulate
#' @param saveTimes generations at which to save population for later
#'     analysis
#' @param nSaveTimes number of times in \code{saveTimes}
#' @param writeInterval number of \code{saveTimes} between writing to disk
#' @return A \code{\link{matrix}}, with one column for each individual
#'     in the population, and one row for each saved generation.  The
#'     matrix entries are the strategies of the individuals.
#' @details By default, \code{nSaveTimes} generations are saved,
#'     equally spaced between \code{0} and \code{nGen}, unless
#'     \code{nGen < nSaveTimes}, in which case every generation is
#'     saved.  If a vector of \code{saveTimes} is given then
#'     \code{nSaveTimes} is ignored.
#' @seealso \code{\link{selection_step_WF}}, \code{\link{selection_step_M}},
#'     \code{\link{selection_step_BD}}
#' @export
#' @import sticky
evolve_with_mutation <- function(model, mutProb = 0.01, mutSD = 0.005, mutRange = c(-Inf,Inf),
                                 xMin = 0, xMax = 6,
                                 xInitMean = 5, xInitSD = mutSD,
                                 nGen = 1000,
                                 nSaveTimes = nGen,
                                 writeInterval = 10,
                                 saveTimes = NULL) {

    if (is.null(saveTimes))
        saveTimes <- floor(seq(from=0, to=nGen, length.out = min(nGen, nSaveTimes)))
    writeTimes <- saveTimes[(1:length(saveTimes)) %% writeInterval == 0]
    cat("writeTimes: ", writeTimes, "\n")
    saveFileName <- paste0("selmutWF", Npop, ".RData")

    pop0 <- initialize_population_randomly(
        model, xMean=xInitMean, xSD=xInitSD, mutRange=mutRange, xMin=xMin, xMax=xMax )

    ## row iGen will contain the population strategies at generation iGen
    Npop <- model$parms$Npop
    popArchive <- matrix(data = NA,
                         ncol = Npop + 1,
                         nrow = length(saveTimes),
                         dimnames = list(NULL,c("iGen", paste0("a",1:Npop))))

    attr(popArchive, "model") <- model
    attr(popArchive, "mutProb") <- mutProb
    attr(popArchive, "mutSD") <- mutSD
    attr(popArchive, "mutRange") <- mutRange
    attr(popArchive, "xMin") <- xMin
    attr(popArchive, "xMax") <- xMax
    attr(popArchive, "xInitMean") <- xInitMean
    attr(popArchive, "xInitSD") <- xInitSD
    attr(popArchive, "nGen") <- nGen
    attr(popArchive, "nSaveTimes") <- nSaveTimes
    attr(popArchive, "saveTimes") <- saveTimes
    popArchive <- add_system_attributes( popArchive )
    class(popArchive) <- c("evolutionData", "matrix")

    popArchive <- sticky(popArchive)

    ## initialize CPU time
    pt0 <- proc.time()

    if (saveTimes[1] == 0) {
        popArchive[1,] <- c(0, pop0[,"strategy"]) # save initial state
        iSave <- 2 # index of saveTimes
    } else {
        iSave <- 1
    }
    iWrite <- 1 # index of writeTimes

    pop <- pop0

    cat("iSave\tiGen\n")

    for(iGen in 1:nGen) {
        # To Do: It should be possible to modify the the selection process (not just WF).
        pop <- selection_step_WF( model, pop )

        pop <- mutation_step_nonoverlapping( model, pop, mutProb, mutSD, mutRange, xMin, xMax )

        if (iGen == saveTimes[iSave]) {
            ##cat(iSave, "\t", iGen, "\n")
            popArchive[iSave,] <- c(iGen, pop[,"strategy"])

            pt <- summary(proc.time() - pt0)
            attr(popArchive,"procTime") <- pt # update cpu time

            ##if (nGen == writeTimes[iWrite]) {#FIX: why doesn't this work???
            ## ... I may have 'tried it' without recompiling...
            if (iSave %% writeInterval == 0) {
                cat("WRITE ", iWrite, ", iSave", iSave, ", iGen", iGen, "\n")
                ##save_simulation_progress(out, iRes, iMut) # FIX: need new function here
                save(file=saveFileName, popArchive, iGen, iSave, iWrite)
                iWrite <- iWrite + 1
            }

            iSave <- iSave + 1
        }
    }

    return(popArchive)
}

#' Mutation step for Wright-Fisher simulations
#'
#' @inheritParams dbenefit
#' @inheritParams evolve_with_mutation
mutation_step_WF <- function( model, pop, mutProb, mutSD, mutType="norm", xMin, xMax ) {
    # we might want to call this function "mutation_step_nonoverlapping_generations"
    # or something to that effect, as it works for any model with nonoverlapping generations
    Npop <- model$parms$Npop

    whichAgentsMutate <- (runif(n=Npop, min=0, max=1) < mutProb)

    # This unnecessarily generates Npop random numbers;
    # we only need mutations for those individuals who mutate.
    mutations <- rnorm(n = Npop, mean = 0, sd = mutSD)
    mutations <- ifelse(whichAgentsMutate, mutations, 0)


    pop[,"strategy"] <- pop[,"strategy"] + mutations

    ## FIX: same 2 lines of code appear in initialize_population_randomly:
    ##      maybe better to sample from truncnorm via truncnorm package,
    ##      except that it is a different truncation for each strategy...
    pop[,"strategy"] <- ifelse(pop[,"strategy"] < xMin, xMin, pop[,"strategy"])
    pop[,"strategy"] <- ifelse(pop[,"strategy"] > xMax, xMax, pop[,"strategy"])

    return(pop)
}

#' Mutation step for simulations with non-overlapping generations
#' (e.g., Wright-Fisher).
#'
#' @inheritParams dbenefit
#' @inheritParams evolve_with_mutation
mutation_step_nonoverlapping <- function( model, pop, mutProb, mutSD, mutRange, xMin, xMax ) {
  Npop <- model$parms$Npop

  whichAgentsMutate <- which((runif(n=Npop, min=0, max=1) < mutProb)==TRUE)
  pop <- mutate_agent_strategies(model, pop, whichAgentsMutate, mutProb, mutSD, mutRange, xMin, xMax )
  return(pop)
}

#' Mutates strategies of a set of agents
#'
#'
#' @param whichAgentsMutate vector of agents whose strategies mutate
#' @inheritParams dbenefit
#' @inheritParams evolve_with_mutation
#'
#' @import truncnorm
mutate_agent_strategies <- function(model, pop, whichAgentsMutate, mutProb, mutSD, mutRange, xMin, xMax ) {

  Npop <- model$parms$Npop

  nMutants <- length(whichAgentsMutate)
  if (nMutants==0){return(pop)} #break if no mutants

  initStrats <-pop[whichAgentsMutate,"strategy"]

  #mutRanges will contain the possible mutation range for each of the strategies to be mutated
  mutRanges <- matrix(0, nrow=nMutants, ncol=length(mutRange))
  mutRanges[,1] <- pmax(xMin-initStrats,mutRange[1]) #mutated strategies cannot be smaller than xMin
  mutRanges[,2] <- pmin(xMax-initStrats,mutRange[2]) #mutated strategies cannot be larger than xMax

  mutations <- rtruncnorm(nMutants, a=mutRanges[,1], b=mutRanges[,2], mean = 0, sd = mutSD) #uses "truncnorm" package
  #mutations[agentInd] <- trandn(l=mutRanges[,1],u=mutRanges[,2]) #uses "TruncatedNormal" package

  pop[whichAgentsMutate, "strategy"] <- initStrats + mutations

  return(pop)
}



#' Plot method for \code{evolutionData} objects
#'
#' @param x an \code{evolutionData} object
#' @importFrom sfsmisc eaxis
#' @export
plot.evolutionData <- function(x,
                               pch=".",
                               las=1,
                               col="grey",
                               ... ) {
    generation <- x[,1]
    strategies <- x[,-1]
    matplot(generation, strategies, pch=pch, las=las, col=col,
            xaxt="n", # use eaxis instead
            ...)
    sfsmisc::eaxis(1, at.small=FALSE, lab.type=if (dev_is_tikz()) "latex" else "plotmath")
}
