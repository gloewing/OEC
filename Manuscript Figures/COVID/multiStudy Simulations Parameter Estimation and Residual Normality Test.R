# Gabe Loewinger
# December 30, 2020
# Multi Study Mortality COVID -- Parameter estimation for simulations

setwd("~/Desktop/Research")
source("msMort_functions.R")
source("Study Strap Functions.R")
source("OEC Functions.R")

pM <- expand.grid( c(FALSE), #nH
             c(100), # minTrain
             c(4), # xPro
             c(1)) # sampSzWt

itr <- 1

sH <- c("New Zealand", "Chile", "Australia DCD")

# test set date ranges
month_start <- "-01-01" #
month_end <-  "-12-31" #
addYear <- FALSE #FALSE # addyear means add a year to the test date--if start and end are the same year then set to FALSE
testYears <- 2020 # use all data

minTrain <- pM[itr, 2] # minimum number of training data points to be considered a training country for any iterations
minTest <- 50 # minimum number of testing data points to be considered a test country for any iterations
testDate_start <- "-01-01" #"-03-01" 
testDate_end <- "-12-31" #"-02-28"
testCountryData <- 12 # this is the number of months that we pretend the test country has (different than minTest)
nH <- pM[itr, 1] # northern hemisphere
# to emulate south africa

full <- read.csv("world-mort.csv")
full$date <- as.Date(full$date)

if(nH){
    # northern hemisphere countries
    countries <- unique(full$country)
    countries <- setdiff(countries, sH)
}else{
    countries <- unique(full$country)
}

K <- length(countries)
# remove countries who dont have TRAINING data on time trained
# ****** these could still be included as test countries but here I didnt for simplicity -- REVISIT with Rolando

OECeta <- 0.5 #0.99 # 0.01
sampSzWeight <- pM[itr, 4]
AvgW <- AvgWEtaTn <- FALSE
etaTune <- "cv"
etaTuneInd <- TRUE
# AvgWEtaTn <- TRUE # whether to tune eta using average weights
# lambdaZero <- TRUE # used to use this to determine whether stacking weights, w, were zero-- now use stackCV to determine if tune with "cv", "hoso" or "zero" -- which just sets it to 0
tuneParamZero <- FALSE # set all study specific alphas to 0
etas <- sort( c(0.001, 0.01, seq(0.05, 0.9, 0.05)) ) #, 0.99, 0.999) ) # hyperparameter for stacking and SSL losses
stackTol <- 1e-9 # tolerance for stacking regression for PGD
oecTune <- "zero"#"cvCFTS" #"cvOEC" # "cvCFTS"#"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
# added in
stdTn <-  "zero"#"cvCFTS" #"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
############
stackCV <- "cv" #"cv" # can be "1", "cv", "hoso" or "zero." If 1, then just use the same parameter as in standard vanilla stacking tuned above. If "zero", then set stacking hyperparamter to 0 if using window stacking
stackCVSpec <- FALSE #"cvSpec"# "cvSpecTS" # "cvSpecTSHOO" # FALSE # "cvSpecTS" # 
stackCVSpec0 <- FALSE # "cvSpec0" #"cvSpecTS0" # "cvSpec0TSHOO"
etaTuneIndSpec <- "cvSpec"# "cvSpecTS" # "cvSpecTSHOO" #TRUE # FALSE
etaTuneIndSpec0 <- "cvSpec0"# "cvSpecTS0" # "cvSpec0TSHOO" # FALSE
horizon <- 10
scaleInd <- TRUE #"spec" # scale covaraites
glmnetOpt <- FALSE
itrs <- 20 # number of runs of OEC
lowLimTune <- FALSE
TScv <- FALSE #TRUE # tune stacking with caret TS instead of cv.glmnet
psiVec <- seq(0,1, length = 11)
pcaInd <- FALSE
tnOrder <- FALSE # determines whether tuning eta or mu first

intercept <- TRUE
xProcess <- 1 # each country is processed separately
psiL <- 0
psiH <- Inf
simplexInd <- FALSE

save.folder <- "/n/home12/gloewinger/mort14"
load.folder <- "~/Desktop/Research"

library(tidyverse)
library(lubridate)
library(splines)
library(glmnet)
library(CVXR)
library(foreach)
library(doParallel)

######### CHANGE NUMBER OF THREADS BEFORE RUNNING REAL THING !!!!!!!!!!!!!!!!!!!!!!!
simNum <- 1

timeStart <- Sys.time()


   cnt <- 1
    print(cnt)
    set.seed(cnt)
    
    library(tidyverse)
    library(lubridate)
    library(splines)
    library(glmnet)
    library(CVXR)
    
    resMat <- matrix(nrow = length(testYears), ncol = 44)
    # rownames(resMat) <- countries
    colnames(resMat) <-  c("merge", "avg", "stacking", "stacking_zeroOut", 
                           "stacking1yr_country", "stacking1yr_zeroOut_country", 
                           "stacking_country", "stacking_country_zeroOut",
                           "country", "oec", "oec2", "stacking1yrALL", "stacking1yrALL_0",
                           "avg", "oecNoInt", "oec_country", "oec_country2", "oec_country3",
                           "oec_country4", "oec_country5",
                           "oec_countryAvg",
                           "eta", "country", "testYear", "K", "n_k", "n_test", "N",
                           "oec_countryAvgStack", "oec_AvgStack", "oec_country0", "oec_country6",
                           "eta1", "etaSpec", "etaSpec0", "mu", "muSpec", "muSpec0", 
                           "muStand", "muStand_Spec", 
                           "oec_SpecUp", "oecWindow", "up", "upWindow")
    
    # 1) chop off all data past covid start
    # 2) chop off all countries with less than 1 year of data
    # 3) DO NOT REMOVE test country
    # 4) fit models with linear model (no ridge penalty)
    # 5) save models
    # 6) estimate mean and variance of model (country)-specific error terms
    # 7) estimate covariance matrix of betas
    # 8) estimate variance of random effects
    # 9) To simulate betas--come up with a mean and covariance matrix for betas and 
    #    simualte study wise (i.e., draw from p-dimensional normal)
    # 10) estimate random effects matrix by centering t(B) matrix and find covariance matrix B B^T so its study-wise covariance
    # 11) use that as the "noise" to add on top of K-dimensonal normal (mean 0) and covariance matrix estimated from centered B B^T
    
    # use only the last test year (do not actually iterate through)

        iterNum <- 1
        full <- read.csv("world-mort.csv")
        full$date <- as.Date(full$date)
        
        
        if(nH){
            # northern hemisphere countries
            countries <- unique(full$country)
            countries <- setdiff(countries, sH )
        }else{
            countries <- unique(full$country)
        }
        
        
        K <- length(countries)
        
        taskNm <- paste0("country_", countries[cnt], "_year: ", iterNum)
        print(paste0("begin ", taskNm))
        
        testDate_start <- paste0(testYears[iterNum], month_start)
        
        if(addYear){
            # add a year onto (if start and end go across calender years )
            testDate_end <- paste0(testYears[iterNum] + 1, month_end)
            # resMat[iterNum, 24] <- testYears[iterNum] + 1 # test year
        }else{
            # start and end are the same calender year
            testDate_end <- paste0(testYears[iterNum], month_end)
            # resMat[iterNum, 24] <- testYears[iterNum] # test year
        }
        

        # iterate through test years
        test_country <- countries[cnt]
        
        # read data
        full <- read.csv("world-mort.csv")
        full$date <- as.Date(full$date)
        
        full <- full[is.element(full$country, countries),] # only include the countries selected above
        
        origCountries <- unique(full$country)
        
        # indices for training and test set
        T_indx <- which( full$date <= testDate_end &  full$date >= testDate_start)
        train_indx <- which(full$date < testDate_start) # all times strictly before start times
        
        # remove countries that do not have enough training data points to be considered a training country
        countriesTrain <- names(which(table(full$country[train_indx]) >= minTrain)) 

        countriesInt <- countriesTrain  
        
        if(length(countriesInt) < 2)   countriesInt <- c()
        # ensure the test country has enough training and testing data
   
        countries <- countriesInt # only include countries that have both
        full <- full[is.element(full$country, countries),] # only include countries with both test and training data
        K <- length(countries) # number of training countries
        countries <- unique(full$country) # resave it so it is not alphabetized but comes in order it appears
        
        # rename outcome
        nm <- which(names(full) == "outcome")
        names(full)[nm] <- "Y"
        
        # test set
        # specific country between certain dates
        indx <- which(full$country == test_country) # test country
        indxDates <- which( full$date <= testDate_end & full$date >= testDate_start) # test dates
        indx <- intersect(indx, indxDates) # test country during test dates
        
        # save a vector of dates for use in simulations
        train_indx <- which(full$date <= testDate_end) # all times strictly before start time of test period
        full <- full[train_indx,] # training data is all countries (with appropriate training and test data) *before* start of test period
        dates <- full$date[full$country=="Finland"] # use Finland because goes baxk the furthest
        
        # training set 
        train_indx <- which(full$date < testDate_start) # all times strictly before start time of test period
        full <- full[train_indx,] # training data is all countries (with appropriate training and test data) *before* start of test period
        
        # save training and test info
        resMat[iterNum, 25] <- K # number of training countries including the test country
        resMat[iterNum, 28] <- nrow(full) # size of training set
        
        if(addYear){
            # add a year onto (if start and end go across calender years )
            resMat[iterNum, 24] <- testYears[iterNum] + 1 # test year
        }else{
            # start and end are the same calender year
            resMat[iterNum, 24] <- testYears[iterNum] # test year
        }
        
        #################
        # Study Labels
        #################
        full$country <- as.numeric(as.factor( full$country) ) # turn numeric
        studyVec <- countries <- unique(full$country)
        num.trainStudy <- length(studyVec)
        indxVec <- vector(length = (nrow(full)))
        
        for(j in 1:K){
            
            indx <- which(full$country == studyVec[j])
            indxVec[indx] <- j
            
        }
        
        full$country <- indxVec # replace with new study labels from 1:K
        countries <- unique(full$country) # number corresponding to it
        
        #################
        # Covaraite Processing
        #################
        
        full <- full[order(full$date),] # reorder by date
        colnames(full)[nm] <- "outcome"
        X <- Xgen(full)[,-c(1,2)] # subtract first two columns because they correspond to splines
        colnames(full)[nm] <- "Y"
        timeLinear <- as.numeric(full$date)
        timeLinear <- timeLinear - min(timeLinear) # make first timne point 0 for ease
        X <- cbind(timeLinear, X) # add linear
        rm(timeLinear)
        
        
        colnames(full)[nm] <- "Y"
        
        # convert outcome to scaled percentage
        full$counts <- full$Y # original outcome
        full$Y  <- full$Y / full$population * 52 * 1000
        
        # indices of rows corresponding to test country but that are in the training set
        trainIndx <- which(full$country == test_country)
        
        ####################################
        # scale covariates
        ####################################
        
        ####################################
        # Make Design Matrix with Full
        ####################################
        if(intercept){
            X <- cbind(1, X)
        }
        
        
        fullDates <- full$date # save dates
        full2 <- as.data.frame( cbind(full$country, full$counts, full$population, X) )
        full <- as.data.frame( cbind(full$country, full$Y, X) )
        if(intercept)           full <- full[,-3] # remove column of ones
        colnames(full) <- c("Study", "Y", 1:(ncol(full) - 2))
        colnames(full2) <- c("Study", "Y", "pop", 1:(ncol(full) - 2))
        
        set.seed(cnt) 
        
        ############
        # stacking 
        ############
        
        betas <- matrix(nrow = ncol(X), ncol = K)
        errorVec <- vector(length = K) # kth entry is the variance of the residual (fitted values)
        timeVec <- vector(length = K) # kth entry is the number of weeks of training data for kth country
        rmse_mat <- matrix(NA, nrow = K, ncol = 3)
        aic_mat <- matrix(NA, nrow = K, ncol = 2)
        
        set.seed(cnt) 
        
        for(k in 1:K){
            # country specific rows
            indx <- which(full$Study == k) # changed for indices 9/11/20
            len <- length(indx)
            
            # linear model
            fit <- lm( as.vector( full$Y[indx] ) ~ as.matrix( X[indx,-1] ) )
            fit2 <- lm( as.vector( full2$Y[indx] ) ~ as.matrix( X[indx,-1] ) )
            fit_pois <- glm( as.vector( full2$Y[indx] ) ~ as.matrix( X[indx,-1] ) + 
                                 offset( log( full2$pop[indx] ) ), 
                             #offset = full2$pop[indx],
                             family = poisson(link = "log"))
            
            # save coefficients
            betas[,k] <- as.numeric( coef(fit) ) # betas of country specific model
            errorVec[k] <- var( fit$residuals ) # time 
            timeVec[k] <- length(indx)
            
            aic_mat[k,1] <- AIC(fit2)
            aic_mat[k,2] <- AIC(fit_pois)
            
            ##############################
            # prediction performance
            ##############################
            # test on last year of data
            
            # test and train indices
            indx2 <- indx[ -seq(1, len - 52) ] # removed last year
            test_indx <- indx[seq(len - 52 + 1, len)]
            
            # linear model
            fit <- lm( as.vector( full$Y[indx2] ) ~ as.matrix( X[indx2,-1] ) )
            fit2 <- lm( as.vector( full2$Y[indx2] ) ~ as.matrix( X[indx2,-1] ) )
            fit_pois <- glm( as.vector( full2$Y[indx2] ) ~ as.matrix( X[indx2,-1] ) + 
                                 offset( log( full2$pop[indx2] ) ), 
                             #offset = full2$pop[indx],
                             family = poisson(link = "log"))
            
            # predictions
            p1 <- predict(fit, as.data.frame( X[test_indx,-1] ) )
            p2 <- predict(fit2,  as.data.frame( X[test_indx,-1] ) )
            p3 <- predict(fit_pois,  as.data.frame( X[test_indx,-1] ) )
            
            p1 <- p1 * full2$pop[test_indx] / ( 52 * 1000 ) # rescale
            
            # rmse
            rmse_mat[k,1] <-  sqrt( mean(  ( p1 -  full2$Y[test_indx] )^2  ) )
            rmse_mat[k,2] <-  sqrt( mean(  ( p2 -  full2$Y[test_indx] )^2  ) )
            rmse_mat[k,3] <-  sqrt( mean(  ( exp(p3) -  full2$Y[test_indx] )^2  ) )
            
        }
        
        colMeans(rmse_mat/rmse_mat[,3])
        
        library(latex2exp)
        par(mfrow = c(1,2))
        boxplot(rmse_mat[,1] / rmse_mat[,3], 
                ylab = TeX('$\\mathbf{RMSE_{OLS}/RMSE_{Poisson}}$')
                )
        
        boxplot(aic_mat[,1] / aic_mat[,2], 
                ylab = TeX('$\\mathbf{AIC_{OLS}/AIC_{Poisson}}$')
        )
        
        # saved as
        # pdf 4.00  7.00
        # ../Final Figures/Final Covid Figures
        # Poisson_vs_ols
        rm(fit, fit2, fit_pois)
        
        ################
        # results
        ################
        result <- list()
        result$betaCov <- cov(t(betas))
        result$randomEff <- cov(betas)
        result$betaMeans <- rowMeans(betas)
        result$residuals <- errorVec
        result$time <- timeVec
        result$betas <- betas
        result$dates <- dates # use for simulating data

        save(result, 
             file = "mortality_parameters_resids.rda")
            

        