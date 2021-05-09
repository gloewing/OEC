############################
# multi country mortality
############################

multiStudySim_mort <- function(sim.num = 130,
                          iter = c(),
                          sampSize = NA, #if NA then randomly draw sample sizes according to dates given in mort data
                          min.time = 100, # minimum number of weeks for training dataset
                          beta.var = 1, # variance of betas
                          beta.meanLinear = 1, # mean of betas for linear term
                          beta.varLinear = 1, # variance for beta of linear term
                          minDate = "2003-01-01", # minimum date for a training set to start on
                          clust.num.mu = 80, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                          clust.num.beta = NULL, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                          sigmaBeta = diag(num.covariates + 1),
                          fixed.effects = c(), # indices of fixed effects
                          fixed.effectsX = c(), # indices of "fixed covariates" across studies
                          studyNoise = NA, #c(1, 1), # if NA, use the noise from the mort data, otherwise sample uniformly with specified amout
                          num.studies = 18,
                          testStudy = num.studies, # arbitrarily set test study to the last one for sample size purposes
                          testSampSz = 52, # length of training set for specialist country
                          testSet.length = 52, # length of test set for specialist country
                          exchg.testStudy = 18, # study index of exchangable study
                          exchg.testSet.length = 52, # number of observations in test set of exchangable study
                          testDateStart = "2019-03-01", # start of "test date" as if we were actually using mortality data for hold-test-period-out CV
                          beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                          perturb = beta.var * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                          SB = 1,
                          params = TRUE,
                          sigmaDiag = FALSE, # if true then the covariance matrix of covariates (X and Z) is diagonal
                          sigmaIdentity = FALSE, # if true then covariance matrix is the identity matrix
                          Xmeans_constant = FALSE, # if TRUE then the means of all covariates in a study are shifted by the same amount
                          Xmean0 = FALSE # IF TRUE then the marginal distribution of the Xs is mean 0
){
    
    library(MASS) # for multivariate normal
    library(lubridate)
    source("msMort_functions.R") # use for fourier series
    load("mortality_parameters.rda") # has parameters
    
    # minimum training dates -- remove all dates before minDate reducing the maximum training period
    result$dates <- result$dates[result$dates >= minDate]
    
    # divide up test dates and training dates
    testDates <- result$dates[ result$dates >= testDateStart ] # save the ones after the start of the test dates
    result$dates <- result$dates[ result$dates < testDateStart ] # remove the test dates from training dates
    
    # add dates to test dates in case you want to test more than 50 dates int est set
    testDates <- c( testDates, testDates + years(1:2) ) # allows you to test up to two more years-- should be more than enough
    
    # if(num.studies > ncol(result$betas)){
    #     # use at most as many studies there are in the real data
    #     message("Too many studies: using number in mortality data")
    #     num.studies <- ncol(result$betas)
    #     
    # }     
    n <- num.covariates <- nrow(result$betas) # use real data for p
    p <- qr.Q(qr(matrix(rnorm(n^2), n)))
    sig.vec <- abs(rnorm(n)) # diagonal elements so variances of the covariates
    
    #### Varies between Studies
    # if( length(beta.var) != num.covariates + 1){
    #     # if the length of the vector does not match, just use the first element
    #     message("Length of beta.var not equal to number of model coefficients: only the first element used for all coefficients")
    #     beta.var <- rep(beta.var[1], num.covariates + 1)
    # }
    
    ############################################################################
    # Uniformly draw number of time points for each study/country
    ############################################################################
    # uniformly draw starting week 
    maxW <- length(result$dates) - min.time # leave sufficient time which is full time minus the minumum
    studyRanges <- matrix( nrow = num.studies, ncol = 2 )
    studyRanges[,1] <- sample(1:maxW, 
                              size = num.studies, 
                              replace = TRUE)
    colnames(studyRanges) <- c("start", "end")
    
    for(j in 1:num.studies){
        # ensure theres enough training
        # ensure it doesnt go over limit
        studyRanges[j, 2] <- studyRanges[j, 1] + sample(min.time:(length(result$dates) - studyRanges[j, 1] - 1),
                                                       size = 1
                                                       ) 
    }
    
    # specialist test country's training data
    studyRanges[testStudy, 1] <- (length(result$dates) - testSampSz) + 1  # fixed amount of data for specialist test dataset's training data
    studyRanges[testStudy, 2] <- length(result$dates)
    
    # exchangable test country's training data
    studyRanges[exchg.testStudy, 1] <- 1  # fixed amount of data for exchangable test dataset's training data
    studyRanges[exchg.testStudy, 2] <- exchg.testSet.length 
    
    # number of observations per study (all studies)
    obs.vec <- studyRanges[, 2] - studyRanges[, 1] + 1
    
    if(is.numeric(sampSize)){
        # 
        obs.vec <- rep(sampSize[1], num.studies)
        
    }
    
    # otherwise use the vector of sample sizes provided if it is the correct length (i.e., as long as num.studies)
    
    # row corresponds to study and column corresponds to True Beta for corresponding covaraite
    beta.matrix <- matrix(NA, nrow = num.studies, ncol = num.covariates )
    
    # fixed across covariates and used as a dial for between study variability
    
    # fix the covariance matrix and mean of the hyperdistribution for the Betas constant
    beta.Sigma <- result$betaCov
    
    # generate different mean vector to generate betas from for MVN
    beta.mean.vec <- result$betaMeans # use means from mort data
    beta.mean.vec[2] <- beta.mean.vec[2] * beta.meanLinear # scale the variance of the random effect for the linear term of time specifically
    ###
    # generate matrix of true betas
    ###
    #scale just the variances of the random effects covariance matrix -- this is what determines between study heterogeneity
    beta.Sigma <- beta.var[1] * beta.Sigma
    beta.Sigma[2,2] <- beta.varLinear * beta.Sigma[2,2] # scale variance of linear term of time
    beta.matrix <-  mvrnorm(n = num.studies, 
                                     mu = beta.mean.vec, 
                                     Sigma = beta.Sigma)

    # # random effects
    # perturb.mat <- mvrnorm(n = num.covariates,
    #                        mu = rep(0, num.studies),
    #                        Sigma = result$randomEff
    #                        )
    
    beta.matrix <- beta.matrix # + t( perturb.mat )
    
    
    # if(length(fixed.effects) > 0){
    #     # fix betas across studies for fixed effects
    #     # if fixed.effects = 0, that is intercept
    #     fixed.effects <- fixed.effects + 1
    #     
    #     indxS <- sample.int(nrow(beta.matrix), 1)
    #     beta.matrix[, c(fixed.effects)] <- t( replicate(nrow(beta.matrix), beta.matrix[indxS, c(fixed.effects)]) ) # just arbitrarily choose one of the betas
    # }
    
    # different variance levels of \epsilon for different studies
    if(is.numeric(studyNoise)){
        # if given, sample noise uniformly
        noiseVec <- vector(length = num.studies)
        for(j in 1:num.studies){
            
            noiseVec[j] <- runif(1, studyNoise[1], studyNoise[2])
            
        }
        
    }else{
        # if not specified then use the residuals from mortality data (sampled with replacement)
        noiseVec <- sample(result$residuals, 
                           size = num.studies, 
                           replace = TRUE)
        noiseVec <- sqrt( noiseVec ) # make it so it is on standard deviation scale to sample from rnorm() below
        
    }

    ####################
    # training studies
    ####################
    
    # Simulate data with above parameters
    min.timeLinear <- min( as.numeric( result$date  )   ) # smallest time point of all dates for linear term -- use to subtract of all time points of all studies
    
    datList <- vector(length = num.studies, "list")
    
    for(y in 1:num.studies){
        data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
        
        
        ##############################
        # Generate the design matrix
        ##############################
        # make 
        dateIndices <- (studyRanges[y,1]):(studyRanges[y,2]) # indices of dates to use
        
        if(y == exchg.testStudy){
            # use different dates because exhangable test study is after training dates end
            dateRange <- testDates[dateIndices]  # use specific dates from mort data
        }else{
            dateRange <- result$dates[dateIndices]  # use specific dates from mort data
            
        }

        data.sim <- data.frame(outcome = rep(1000, obs.vec[y]), # outcome is arbitrary and only used for XGen code
                               date =  dateRange,
                               population = rep(1000, obs.vec[y]) # population is arbitrary and only used for XGen code
                               )
        
        X <- suppressMessages( Xgen(data.sim)[,-c(1,2)] )# subtract first two columns because they correspond to splines
        timeLinear <- as.numeric(dateRange)
        timeLinear <- timeLinear - min.timeLinear # make first timne point 0 for ease
        X <- cbind(timeLinear, X) # add linear
        rm(timeLinear)
        
        # generate a vector Y and add noise
        Y <- cbind(1, X) %*% beta.matrix[y, ] + rnorm(obs.vec[y], 0, noiseVec[y]) # noise is mean 0 with study specific noise levels
        
        # bind it to data
        data.sim <- cbind(y, Y, X) # first column is study number, then Y outcome, then design matroix
        colnames(data.sim) <- c("Study", "Y", paste0("V_", 1:(num.covariates - 1) ) )
        
        if(y == exchg.testStudy){
            # save as exhangable test set
            exchngTest <- data.sim
        }else{
            # save to list in order to merge studies
            datList[[y]] <- data.sim
            
            }
            
        
    }
    
    
    ##############
    # test set
    ##############
    
    #****************************
    # Generate the design matrix
    #****************************
    data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
    # make test dataset
    y <- testStudy # set index to test study
    dateIndices <- 1:testSet.length # indices of dates to use
    dateRange <- testDates[dateIndices] #result$dates[dateIndices]  # use specific dates from mort data
    
    data.sim <- data.frame(outcome = rep(1000, testSet.length), # outcome is arbitrary and only used for XGen code
                           date =  dateRange,
                           population = rep(1000, testSet.length) # population is arbitrary and only used for XGen code
    )
    
    X <- suppressWarnings( Xgen(data.sim)[,-c(1,2)] ) # subtract first two columns because they correspond to splines
    timeLinear <- as.numeric(dateRange)
    timeLinear <- timeLinear - min.timeLinear # make first timne point 0 for ease
    X <- cbind(timeLinear, X) # add linear
    rm(timeLinear)
    
    # generate a vector Y and add noise
    Y <- cbind(1, X) %*% beta.matrix[y, ] + rnorm(testSet.length, 0, noiseVec[y]) # noise is mean 0 with study specific noise levels
    
    # bind it to data
    testData <- cbind(y, Y, X) # first column is study number, then Y outcome, then design matroix
    colnames(testData) <- c("Study", "Y", paste0("V_", 1:(num.covariates - 1) ) )
    
    # +++++++++++++++++++++++++++++++++++++++++++
    #######################
    # concatenate data
    #######################
    final.results <- do.call(rbind, datList) # merge studies
    colnames(final.results) <- c("Study", "Y", paste0("V_", 1:(num.covariates - 1) ) )
    
    # write.csv(final.results, paste0(filename, "_Combined"))
    message(paste0("Data Simulation Complete"))
    
    if(params){
        return( list(data = final.results, 
                     test = testData, # specialist test set
                     exchngTest = exchngTest, # exchangable test set
                     betaMeans = beta.mean.vec, 
                     betas = beta.matrix, 
                     xMean = NA, # mean.vec,
                     Xs = NA, #mu.matrix,
                     Sigma = NA, #Sigma,
                     d = NA #d
                     )
                )
    }else{
        return( final.results )
    }
    
}


