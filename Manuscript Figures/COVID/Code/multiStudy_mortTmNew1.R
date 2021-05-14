# Gabe Loewinger
# April 27, 2021
# Multi Study Mortality COVID
# New Mortality dataset

source("msMort_functions.R")
source("Study Strap Functions.R")
source("OEC Functions.R")

pM <- expand.grid( c(FALSE), #nH
                 c(100), # minTrain
                 c(4), # xPro 
                 c(1), # sampSzWt
                 c("cvCF", "zero"), # country specific model tuning 
                 c("cv", "zero") # stackCV
                )

itr <- 1 

# test set date ranges
month_start <-  "-01-01" #
month_end <- "-12-31" # 
addYear <- FALSE # addyear means add a year to the test date--if start and end are the same year then set to FALSE
testYears <- 2003:2019 # chosen because only a couple countries with enough data before 2003 (i.e., 2003 is first year where training K is large enough, before 2003 K < 5)

minTrain <- pM[itr, 2] # minimum number of training data points to be considered a training country for any iterations
minTest <- 50 # minimum number of testing data points to be considered a test country for any iterations
testDate_start <- "-01-01" #"
testDate_end <- "-12-31" #
testCountryData <- 12 # this is the number of months that we pretend the test country has (different than minTest)
nH <- pM[itr, 1] # northern hemisphere

full <- read.csv("world-mort.csv")
full$date <- as.Date(full$date)

# southern hemisphere does not match for estimating seasonal effects
sH <- c("New Zealand", "Chile", "Australia DCD")

if(nH){
    # northern hemisphere countries
    countries <- unique(full$country)
    countries <- setdiff(countries, sH)
}else{
    countries <- unique(full$country)
}

K <- length(countries)

OECeta <- 0.5 # starting oeceta - arbitrary
sampSzWeight <- pM[itr, 4]
AvgW <- AvgWEtaTn <- FALSE
etaTune <- "cv"
etaTuneInd <- TRUE
tuneParamZero <- FALSE # set all study specific alphas to 0
etas <- sort( c(0.001, 0.01, seq(0.05, 0.9, 0.05)) ) # hyperparameter for stacking and SSL losses
stackTol <- 1e-9 # tolerance for stacking regression for PGD
oecTune <- pM[itr, 5] 
# added in
stdTn <-  pM[itr, 5] 
############
stackCV <- pM[itr, 6] 
stackCVSpec <- FALSE 
stackCVSpec0 <-  FALSE 
etaTuneIndSpec <- "cvSpec"
etaTuneIndSpec0 <- "cvSpec0"
horizon <- 10
scaleInd <- TRUE # scale covaraites
glmnetOpt <- FALSE
itrs <- 2 
lowLimTune <- FALSE
TScv <- FALSE #
psiVec <- seq(0,1, length = 11)
pcaInd <- FALSE
tnOrder <- FALSE #
nFolds <- 5
orderRandom <- FALSE

intercept <- TRUE
xProcess <- pM[itr, 3]
psiL <- 0
psiH <- Inf
simplexInd <- FALSE

save.folder <- "/n/home12/gloewinger/mort17"
load.folder <- "~/Desktop/Research"

library(tidyverse)
library(lubridate)
library(splines)
library(glmnet)
library(CVXR)
library(foreach)
library(doParallel)


# tune for stacking
lambdaVec <- tune.grid <- sort( unique( c(0.0001, 0.001, 0.01, 
                             exp(-seq(0,5, length = 50)),
                             seq(5,20, by = 5),
                             seq(30,100, by = 10) ) ) ) # 2:100

tune.grid <- as.data.frame(tune.grid) # tuning parameters to consider
tune.grid <- cbind(0, tune.grid) # Ridge
colnames(tune.grid) <- c("alpha", "lambda")
warmStartSolver <- "glmnet" # NULL will set to my own function

# parallelize
fileNm <- paste0("mrtNew_eta", OECeta, 
                  "_xPro_", xProcess, 
                 "_etaTn_",  
                 etaTune,
                 "_etaSpec_", etaTuneIndSpec,
                 "_Wcv_", stackCV, 
                 "_Wspec_", stackCVSpec,
                  "_Wspec0_", stackCVSpec0,
                 "_TScv_", TScv,
                 "_sclX_", scaleInd,
                 "_glm_", glmnetOpt,
                 "_hr_", horizon,
                 "_oecTn_", oecTune, 
                 "_stdTn_", stdTn,
                 "_tnOr_", tnOrder,
                 "_tstCntMn_", testCountryData,
                 ".nH", nH,
                 "_mnSt_", month_start,
                 "_mnTr_", minTrain,
                 "_fld_", nFolds,
                 "_smpSzWt_", sampSzWeight)
print(fileNm)

######### Parallelize
simNum <- 1
logfile <- paste0("outputFile_mort", simNum,".txt")
writeLines(c(""), file(logfile,'w'))

num.threads <- as.integer( K ) 
threads <- makeCluster(num.threads, outfile=logfile)
registerDoParallel(threads)

setwd(save.folder)

getDoParWorkers()
timeStart <- Sys.time()

results <- foreach(cnt = 1:K, .combine = list, .multicombine = TRUE) %dopar%{
   
    print(cnt)
    set.seed(cnt)
    
    library(tidyverse)
    library(lubridate)
    library(splines)
    library(glmnet)
    library(CVXR)
    
    resMat <- matrix(nrow = length(testYears), ncol = 72)
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
                           "oec_SpecUp", "oecWindow", "up", "upWindow", "country_noLinear",
                           "oec_specAnneal_left", "oec_specAnneal_leftObj", "oec_spec_left", "oec_spec_leftObj",
                           "oec_specAnneal_right", "oec_specAnneal_rightObj", "oec_spec_right", "oec_spec_rightObj",
                           "oec_0Anneal_left", "oec_0Anneal_leftObj", "oec_0_left", "oec_0_leftObj",
                           "oec_0Anneal_right", "oec_0Anneal_rightObj", "oec_0_right", "oec_0_rightObj",
                           "oec_genAnneal_left", "oec_genAnneal_leftObj", "oec_gen_left", "oec_gen_leftObj",
                           "oec_genAnneal_right", "oec_genAnneal_rightObj", "oec_gen_right", "oec_gen_rightObj",
                           "etaGen", "etaSpec", "etaSpec0")
    
    for(iterNum in 1:length(testYears) ){

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
        
        taskNm <- paste0("country_", countries[cnt], "_year: ", iterNum)
        print(paste0("begin ", taskNm))
        
        testDate_start <- paste0(testYears[iterNum], month_start)
        
        if(addYear){
            # add a year onto (if start and end go across calender years )
            testDate_end <- paste0(testYears[iterNum] + 1, month_end)

        }else{
            # start and end are the same calender year
            testDate_end <- paste0(testYears[iterNum], month_end)

        }
        

        # iterate through test years
        test_country <- countries[cnt]
        
        # read data
        full <- read.csv("world-mort.csv")
        full$date <- as.Date(full$date)
        
        full <- full[is.element(full$country, countries),] # only include the countries selected above
        
        # make test country have only the specified amount of training data
        dateStart <- ymd(testDate_start) - months(testCountryData) # this is when training data for test country starts
        dateEnd <- testDate_end
        indx <- which(full$country == test_country) # test country
        
        # remove dates before training period for test country
        indxDatesStart <- which( full$date < dateStart )
        indx1 <- intersect(indx, indxDatesStart) # test country during test dates
        if(length(indx1 > 0))      full <- full[-indx1,]
        rm(indx1)
        
        # remove dates after test period for test country
        indxDatesEnd <- which( full$date > dateEnd ) # test dates
        indx1 <- intersect(indx, indxDatesEnd) # test country during test dates
        if(length(indx1 > 0))      full <- full[-indx1,]
        rm(indx1)
        
        origCountries <- unique(full$country)
        
        # indices for training and test set
        T_indx <- which( full$date <= testDate_end &  full$date >= testDate_start)
        train_indx <- which(full$date < testDate_start) # all times strictly before start times
        
        # remove countries that do not have enough training data points to be considered a training country
        countriesTrain <- names(which(table(full$country[train_indx]) >= minTrain)) 
        
        
        # allow for test country to have as many training observations as we test (in other words, it only 
        # has to have 50, not as many as other training countries which may require more than 50)
        # so just set it to be the minimum of minTest instead of minTrain like other training sets
        test_trainingSet <- which(full$country[train_indx] == test_country) # training set for test country
        if( length(test_trainingSet) >= minTest)       countriesTrain <- unique( c(countriesTrain, as.character(test_country)) ) # it only has to have as many as training set to be included
        rm(test_trainingSet)
        
        
        # remove countries that do not have enough Testing data points to be considered a test country
        countriesTest <- names(which(table(full$country[T_indx]) >= minTest))
        
        countriesInt <- intersect(countriesTest, countriesTrain ) # countries that have enough both training and test data
        if(length(countriesInt) < 2)   countriesInt <- c()
        # ensure the test country has enough training and testing data
        
        if(is.element(test_country, countriesInt) ){
            
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
            
            # test set is the test country during testing period
            test <- full[indx,] 
            
            # training set 
            train_indx <- which(full$date < testDate_start) # all times strictly before start time of test period
            full <- full[train_indx,] # training data is all countries (with appropriate training and test data) *before* start of test period
            
            # indices of rows corresponding to test country but that are in the training set
            trainIndx <- which(full$country == test_country)
            
            # save training and test info
            resMat[iterNum, 23] <- as.character( unique( test$country ) ) # test country
            resMat[iterNum, 25] <- K # number of training countries including the test country
            resMat[iterNum, 26] <- length(trainIndx) # number of observations in training set for test country
            resMat[iterNum, 27] <- nrow(test) # number of observations in training set for test country
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
            test_country <- which(countries == test_country) # number corresponding to test country
            full$country <- as.numeric( as.factor( full$country)  )# turn numeric
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
            # Covaraite Processing - XProcess = 4
            #################
            # generate design matrix
            
            # Linear Trend
            test$country <- test_country # make this into numeric
            
            # order full and test by dates
            full <- rbind(full, test) # recombine into one
            full <- full[order(full$date),] # reorder by date
            colnames(full)[nm] <- "outcome"
            X <- Xgen(full)[,-c(1,2)] # subtract first two columns because they correspond to splines
            colnames(full)[nm] <- "Y"
            timeLinear <- as.numeric(full$date)
            timeLinear <- timeLinear - min(timeLinear) # make first timne point 0 for ease
            X <- cbind(timeLinear, X) # add linear
            rm(timeLinear)
            
            # test set
            # specific country between certain dates
            indx <- which(full$country == test_country) # test country
            indxDates <- which( full$date <= testDate_end & full$date >= testDate_start) # test dates
            indx <- intersect(indx, indxDates) # test country during test dates
            
            # test set is the test country during testing period
            test <- full[indx,] 
            X_test <- X[indx,] # design matrix for test set
            
            full <- full[-indx,] # remove test from full
            X <- X[-indx,] # remove test rows from full design matrix
            #################
            
            # convert outcome to scaled percentage
            full$Y  <- full$Y / full$population * 52 * 1000
            test$Y  <- test$Y / test$population * 52 * 1000
            
            # indices of rows corresponding to test country but that are in the training set
            trainIndx <- which(full$country == test_country)
            
            ####################################
            # scale covariates
            ####################################
            
            if(scaleInd == TRUE){
                
                nFull <- nrow(X) # sample size of merged
                
                # scale Covaraites
                means <- colMeans( as.matrix(X) )
                sds <- sqrt( apply( as.matrix(X), 2, var) *  (nFull - 1) / nFull )  # use mle formula to match with GLMNET
                
                #
                for(column in 1:ncol(X) ){
                    # center scale
                    X[, column] <- (X[, column ] - means[column]) / sds[column]
                    X_test[, column] <- (X_test[, column ] - means[column]) / sds[column]
                }
                
            }else if(scaleInd == "spec"){
                # scale according to rows of training set of test country
                
                
                nFull <- nrow(X[trainIndx,]) # sample size of merged
                
                # scale Covaraites
                means <- colMeans( as.matrix(X[trainIndx,]) )
                sds <- sqrt( apply( as.matrix(X[trainIndx,]), 2, var) *  (nFull - 1) / nFull )  # use mle formula to match with GLMNET
                
                #
                for(column in 1:ncol(X) ){
                    # center scale
                    X[, column] <- (X[, column ] - means[column]) / sds[column]
                    X_test[, column] <- (X_test[, column ] - means[column]) / sds[column]
                }
            }
            
            ####################################
            # Make Design Matrix with Full
            ####################################
            if(intercept){
                X <- cbind(1, X)
                X_test <- cbind(1, X_test)
            }
            
            fullDates <- full$date # save dates
            full <- as.data.frame( cbind(full$country, full$Y, X) )
            if(intercept)           full <- full[,-3] # remove column of ones
            colnames(full) <- c("Study", "Y", 1:(ncol(full) - 2))
            
            set.seed(cnt) # added 9 /3/20
            
            #######################
            # Tuning
            #######################
            # Merge Tune
            set.seed(cnt) # added 9 /3/20
            if(stdTn == "zero" & oecTune == "zero"){
                mergedLambda <- data.frame(alpha = 0, lambda = 0)
            }else{
                mergedLambda <- hosoCV(data = full,
                                       tune.grid,
                                       hoso = "merged",
                                       method = "glmnet",
                                       metric = "RMSE",
                                       nfolds = nFolds, #, "K",
                                       nnlsInd = TRUE,
                                       OECeta = 0.5,
                                       sampSzWeight = sampSzWeight,
                                       weights = NULL)
            }
            
            # OEC Tune
            set.seed(cnt) # added 9 /3/20
            tuneParam <- hosoCV(data = full,
                                tune.grid,
                                hoso = stdTn,
                                method = "glmnet",
                                metric = "RMSE",
                                nfolds = nFolds, #, "K",
                                nnlsInd = TRUE,
                                OECeta = 0.5,
                                sampSzWeight = sampSzWeight,
                                weights = NULL)
            
            if(oecTune != stdTn){
                # OEC Tune
                set.seed(cnt) # added 9 /3/20
                tuneParamOEC <- hosoCV(data = full,
                                       tune.grid,
                                       hoso = oecTune,
                                       method = "glmnet",
                                       lambdaVec = tuneParam$lambda,
                                       metric = "RMSE",
                                       nfolds = nFolds, #, "K",
                                       nnlsInd = TRUE,
                                       OECeta = 0.5,
                                       sampSzWeight = sampSzWeight,
                                       weights = NULL)
                
            }else{
                tuneParamOEC <- tuneParam
            }
            
            
            
            if(oecTune == "oec"){
                e <- tuneParam
                tuneParam <- data.frame(matrix(ncol = 2, nrow = K))
                tuneParam[1:K,] <- e
                colnames(tuneParam) <- colnames(e)
                rm(e)
            }   
            
            #################
            # Merging
            #################

            fit <- glmnet(y = as.vector(full$Y),
                          x = as.matrix(X[,-1]),
                          alpha = 0,
                          lambda = 0,
                          standardize = FALSE,
                          intercept = TRUE, 
                          thresh = 1e-10 )
            
            preds <- X_test %*% coef(fit) # predict(fit, X_test)
            mrgBeta <- coef(fit)
            resMat[iterNum, 1] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            
            ############
            # stacking 
            ############
            
            predsMat <- matrix(nrow = nrow(full), ncol = K )
            predsMat_testCountry <- matrix(nrow = length(trainIndx), ncol = K )
            betas <- matrix(nrow = ncol(X), ncol = K)
            
            set.seed(cnt) # added 9 /3/20
            
            for(k in 1:K){
                # country specific rows

                indx <- which(full$Study == k) # changed for indices 9/11/20

                if(glmnetOpt){
                    fit <- glmnet(y = as.vector(full$Y[indx]),
                                     x = as.matrix(X[indx,-1]),
                                     alpha = 0,
                                     lambda = tuneParam$lambda[k],
                                     standardize = FALSE,
                                     intercept = TRUE, 
                                  thresh = 1e-10)#,
                    b <- as.numeric( coef(fit) )
                    
                }else if(glmnetOpt == 1){
                    # same as above but standardize
                    fit <- glmnet(y = as.vector(full$Y[indx]),
                                  x = as.matrix(X[indx,-1]),
                                  alpha = 0,
                                  lambda = tuneParam$lambda[k],
                                  standardize = TRUE,
                                  intercept = TRUE, 
                                  thresh = 1e-10)#,
                    b <- as.numeric( coef(fit) )
                }else{
                    b <- ridgeEst(y = as.vector(full$Y[indx]),
                                  x = as.matrix(X[indx, -1]), # already includes a column of ones for intercept
                                  intercept = TRUE, 
                                  lambda = tuneParam$lambda[k] ) # standardized like glmnet in cvCF function tuning
                    
                }
                
                # predict on entire training set for standard stacking
                predsMat[,k] <- X %*% b # predict(fit, X)
                
                # predict on entire training set for test country
                predsMat_testCountry[,k] <- X[trainIndx,] %*% b # predict(fit, full[trainIndx,])
                betas[,k] <- b # coefs
                
                # save the model corresponding to test country
                if(k == test_country){
                    
                    ###########################
                    # country specific model
                    ###########################
                    bCountry <- b
                    preds <- X_test %*% b
                    resMat[iterNum, 9] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                    
                    ############################################
                    # country specific model no linear term
                    ############################################
                    if(glmnetOpt){
                        fit <- glmnet(y = as.vector(full$Y[indx]),
                                      x = as.matrix(X[indx, -c(1,2)]),
                                      alpha = 0,
                                      lambda = tuneParam$lambda[k],
                                      standardize = FALSE,
                                      intercept = TRUE, 
                                      thresh = 1e-10)#,
                        b <- as.numeric( coef(fit) )
                        
                    }else if(glmnetOpt == 1){
                        # same as above but standardize
                        fit <- glmnet(y = as.vector(full$Y[indx]),
                                      x = as.matrix(X[indx, -c(1,2)]),
                                      alpha = 0,
                                      lambda = tuneParam$lambda[k],
                                      standardize = TRUE,
                                      intercept = TRUE, 
                                      thresh = 1e-10)#,
                        b <- as.numeric( coef(fit) )
                    }else{
                        b <- ridgeEst(y = as.vector(full$Y[indx]),
                                      x = as.matrix(X[indx, -c(1,2)]), # already includes a column of ones for intercept
                                      intercept = TRUE, 
                                      lambda = tuneParam$lambda[k] ) # standardized like glmnet in cvCF function tuning
                        
                    }
                    
                    # test without the linear term
                    if(xProcess != 5){
                        # do not test removal of linear term for xPro = 5 because that does not have linear term
                        bCountry_noLinear <- b
                        preds <- X_test[,-2] %*% b
                        resMat[iterNum, 45] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                        
                    }
                    
                }   
        
            }
            
            predsMat0 <- predsMat # save copy of predsMat
            rm(fit)
            
            ###################
            # average weights
            ###################
        
            preds <- rowMeans( X_test %*% betas )
            resMat[iterNum, 2] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            
            #######################
            # standard stacking
            ######################
            # tune stacking L2 penalty
            modTune <- cv.glmnet(y = as.vector(full$Y), 
                                 x = as.matrix(predsMat),
                                 alpha = 0, 
                                 lambda = lambdaVec, 
                                 standardize = TRUE,
                                 intercept = TRUE,
                                 lower.limits = 0, 
                                 thresh = 1e-10
            ) 
            
            stackParam <- modTune$lambda.min # replace old value with tuned one
            stackParamOEC <- stackParam # use same paramter for OEC incase OEC is not retuned below
            rm(modTune)

            if(stackCV == "zero")    stackParamOEC <- stackParam <- 0
            
            mod <- glmnet(y = as.vector(full$Y), 
                          x = as.matrix(predsMat), 
                          alpha = 0,         
                          lambda = stackParam, 
                          standardize = TRUE,
                          intercept = TRUE,
                          lower.limits = 0, 
                          thresh = 1e-10 ) 
            
            w <- as.vector( coef( mod , 
                                   exact = TRUE, 
                                   y = as.vector(full$Y), 
                                   x = as.matrix(predsMat)
                                    )
                            )
            
            rm(mod)
            
            # standard stacking predictions
            preds <- w[1] + X_test %*% betas %*% w[-1]
            resMat[iterNum, 3] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            
            stackParamOEC1 <- stackParam # use same paramter for OEC incase OEC is not retuned below
            
            #######################
            # zero out stacking
            ######################
            # zero out
            indxMat <- cbind(trainIndx, test_country)
            predsMat[indxMat] <- 0 # zero out observations corresponding to this study
            
            # tune stacking L2 penalty
            modTune <- cv.glmnet(y = as.vector(full$Y), 
                                 x = as.matrix(predsMat),
                                 alpha = 0, 
                                 lambda = lambdaVec, 
                                 standardize = TRUE,
                                 intercept = TRUE,
                                 lower.limits = 0, 
                                 thresh = 1e-10
            ) 
            
            stackParam <- modTune$lambda.min # replace old value with tuned one
            stackParamOEC <- stackParam # use same paramter for OEC
            rm(modTune)

            if(stackCV == "zero")    stackParamOEC <- stackParam <- 0
            
            mod <- glmnet(y = as.vector(full$Y), 
                          x = as.matrix(predsMat), 
                          alpha = 0,         
                          lambda = stackParam, 
                          standardize = TRUE,
                          intercept = TRUE,
                          lower.limits = 0, 
                          thresh = 1e-10 ) 
            
            w <- as.vector( coef( mod , 
                                  exact = TRUE, 
                                  y = as.vector(full$Y), 
                                  x = as.matrix(predsMat)
            )
            )
            
            rm(mod)
            
            # standard stacking predictions
            preds <- w[1] + X_test %*% betas %*% w[-1]
            resMat[iterNum, 4] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            
            
            #######################
            # standard stacking on year prior
            ######################
            
            # year prior -- dates
            end <- ymd(testDate_end) - years(1)
            start <- ymd(testDate_start) - years(1)
            
            # indices correpsonding to test country in training set for year prior
            indx1 <- which(full$Study == test_country)
            indxDates <- which( fullDates <= end &  fullDates >= start)
            indx1 <- intersect(indx1, indxDates)
            
            # make sure there are observations to stack on 
            if(length(indx1) > 0 ){
                predsMat <- predsMat0[indx1,] # alter stacking matrix for just the same period as test period but 1 year prior
                
                # tune stacking L2 penalty
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)

                if(stackCV == "zero")    stackParam <- 0
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y[indx1,]), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% betas %*% w[-1]
                resMat[iterNum, 5] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                
            }

            #######################
            # zero out stacking on year prior
            ######################
            
            # year prior -- dates
            end <- ymd(testDate_end) - years(1)
            start <- ymd(testDate_start) - years(1)
            
            # indices correpsonding to test country in training set for year prior
            indx1 <- which(full$Study == test_country)
            indxDates <- which( fullDates <= end &  fullDates >= start)
            indx1 <- intersect(indx1, indxDates)
            
            if(length(indx1) > 0 ){
                
                predsMat <- predsMat0[indx1,] # alter stacking matrix for just the same period as test period but 1 year prior
                predsMat[, test_country] <- 0 # zero out entire column because all rows are for test country
                
                # tune stacking L2 penalty
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)
                    

                if(stackCV == "zero")    stackParam <- 0
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y[indx1,]), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% betas %*% w[-1]
                resMat[iterNum, 6] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
        
                
            }
            
            #######################
            # standard stacking on country
            ######################
            # indices correpsonding to test country in training set for year prior
            indx1 <- which(full$Study == test_country)
            
            # tune l2 parameter
            modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                 x = as.matrix(predsMat_testCountry),
                                 alpha = 0, 
                                 lambda = lambdaVec, 
                                 standardize = TRUE,
                                 intercept = TRUE,
                                 lower.limits = 0, 
                                 thresh = 1e-10
            ) 
            
            stackParam <- modTune$lambda.min # replace old value with tuned one
            rm(modTune)
                

            if(stackCV == "zero")    stackParam <- 0
            
            mod <- glmnet(y = as.vector(full$Y[indx1]), 
                          x = as.matrix(predsMat_testCountry), 
                          alpha = 0,         
                          lambda = stackParam, 
                          standardize = TRUE,
                          intercept = TRUE,
                          lower.limits = 0, 
                          thresh = 1e-10 ) 
            
            w <- as.vector( coef( mod , 
                                  exact = TRUE, 
                                  y = as.vector(full$Y[indx1]), 
                                  x = as.matrix(predsMat_testCountry)
            )
            )
            
            wCountry <- w
            
            rm(mod)
            
            # standard stacking predictions
            preds <- w[1] + X_test %*% betas %*% w[-1]
            resMat[iterNum, 7] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            resMat[iterNum, 39] <- stackParam
            
            #######################
            # standard stacking on country zero out
            ######################
            # indices correpsonding to test country in training set for year prior
            indx1 <- which(full$Study == test_country)
            predsMat_testCountry2 <- predsMat_testCountry
            predsMat_testCountry2[,test_country] <- 0 # zero out country thats being tested
            
           # tune l2 parameter
            modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                 x = as.matrix(predsMat_testCountry2),
                                 alpha = 0, 
                                 lambda = lambdaVec, 
                                 standardize = TRUE,
                                 intercept = TRUE,
                                 lower.limits = 0, 
                                 thresh = 1e-10
            ) 
            
            stackParam <- modTune$lambda.min # replace old value with tuned one
            rm(modTune)

            
            if(stackCV == "zero")    stackParam <- 0
            
            mod <- glmnet(y = as.vector(full$Y[indx1]), 
                          x = as.matrix(predsMat_testCountry2), 
                          alpha = 0,         
                          lambda = stackParam, 
                          standardize = TRUE,
                          intercept = TRUE,
                          lower.limits = 0, 
                          thresh = 1e-10 ) 
            
            w <- as.vector( coef( mod , 
                                  exact = TRUE, 
                                  y = as.vector(full$Y[indx1]), 
                                  x = as.matrix(predsMat_testCountry2)
            )
            )
            
            wZeroOut <- w
            
            rm(mod)
            
            # standard stacking predictions
            preds <- w[1] + X_test %*% betas %*% w[-1]
            resMat[iterNum, 8] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
            resMat[iterNum, 40] <- stackParam
           ############################################################################

            #######################
            # standard stacking on year prior ALL COUNTRIES
            ######################
            
            # year prior -- dates
            end <- ymd(testDate_end) - years(1)
            start <- ymd(testDate_start) - years(1)
            
            # indices correpsonding to test country in training set for year prior
            indx1 <-  which( fullDates <= end &  fullDates >= start)
            
            # make sure there are observations to stack on 
            if(length(indx1) > 0 ){
                predsMat <- predsMat0[indx1,] # alter stacking matrix for just the same period as test period but 1 year prior
                
                # tune l2 parameter
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)   
                    

                if(stackCV == "zero")    stackParam <- 0
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y[indx1,]), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% betas %*% w[-1]
                resMat[iterNum, 12] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                
            }
            
            
            #######################
            # zero out stacking on year prior
            ######################
            # year prior -- dates
            end <- ymd(testDate_end) - years(1)
            start <- ymd(testDate_start) - years(1)
        
            # indices correpsonding to test country in training set for year prior
            indx1 <-  which( fullDates <= end &  fullDates >= start)
            
            
            if(length(indx1) > 0 ){
                
                predsMat <- predsMat0[indx1,] # alter stacking matrix for just the same period as test period but 1 year prior
                predsMat[, test_country] <- 0 # zero out entire column because all rows are for test country
                
                # tune stacking L2 penalty
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)

                if(stackCV == "zero")    stackParam <- 0
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y[indx1,]), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% betas %*% w[-1]
                resMat[iterNum, 13] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
        
            }
            rm(predsMat, predsMat0)
            
            #########################
            # OEC
            ########################
            
            ##########################################
            # Tune OEC stacking hyperparameter (for w)
            ##########################################
            # if stackCV == "1", then just use the value above
            # if "zero" then just no penalization for OEC stacking regression
            # if "cv" or "hoso" then tune it in the OEC framework

                if(stackCV == "zero"){
                    
                    stackParamOEC <- 0 # zero out

                }else if(stackCV == "cv" | stackCV == "hoso"){
                    ############################################
                    # Tune Stacking Hyperparameter again with tuned eta
                    ############################################
                    stackParamOEC <- oecW_CV(data = full,
                                             tune.grid = tune.grid,
                                             sslLambdas = tuneParamOEC$lambda,
                                             method = "glmnet",
                                             nfolds = nFolds, # K,
                                             nnlsInd = TRUE,
                                             OECeta = OECeta,
                                             sampSzWeight = sampSzWeight,
                                             glmnetOpt = glmnetOpt,
                                             xStandardize = FALSE, # I just chose this as FALSE 9/2/20
                                             cv = stackCV,
                                             standardize = TRUE,
                                             weights = Diagonal( nrow(full) ),
                                             horizon = horizon )
                    
                    stackParamOEC <- stackParamOEC$lambda
                }
                
                if(stackCVSpec == "cvSpecTS"| stackCVSpec == "hosoSpecTS" | stackCVSpec == "cvSpec" |
                   stackCVSpec == "cvSpecTSHOO"){
                    ############################################
                    # Tune Stacking Hyperparameter again with tuned eta
                    ############################################
                    stackParamOECSpec <- oecW_CV(data = full,
                                                 tune.grid = tune.grid,
                                                 sslLambdas = tuneParamOEC$lambda,
                                                 method = "glmnet",
                                                 nfolds = nFolds, 
                                                 nnlsInd = TRUE,
                                                 OECeta = OECeta,
                                                 sampSzWeight = sampSzWeight,
                                                 glmnetOpt = glmnetOpt,
                                                 xStandardize = FALSE, 
                                                 cv = stackCVSpec,
                                                 standardize = TRUE,
                                                 weights = Diagonal( nrow(full) ),
                                                 SpecID = test_country,
                                                 horizon = horizon )
                    
                    stackParamOECSpec <- stackParamOECSpec$lambda
                }else if(stackCVSpec == "standard"){
                    
                    # use from above (standard stacking) -- tuned above in standard stacking section
                    stackParamOECSpec <- stackParamSpec
                    
                }else{
                    
                    #otherwise use from above generalist OEC 
                    stackParamOECSpec <- stackParamOEC
                }
                
                # zero out specialist
                if(stackCVSpec0 == "cvSpecTS0"| stackCVSpec0 == "hosoSpecTS0" | stackCVSpec0 == "cvSpec0" |
                   stackCVSpec == "cvSpec0TSHOO"){
                    
                    ############################################
                    # Tune Stacking Hyperparameter again with tuned eta
                    ############################################
                    stackParamOECSpec0 <- oecW_CV(data = full,
                                                  tune.grid = tune.grid,
                                                  sslLambdas = tuneParamOEC$lambda,
                                                  method = "glmnet",
                                                  nfolds = nFolds,
                                                  nnlsInd = TRUE,
                                                  OECeta = OECeta,
                                                  sampSzWeight = sampSzWeight,
                                                  glmnetOpt = glmnetOpt,
                                                  xStandardize = FALSE,
                                                  cv = stackCVSpec0,
                                                  standardize = TRUE,
                                                  weights = Diagonal( nrow(full) ),
                                                  SpecID = test_country,
                                                  horizon = horizon )
                    
                    stackParamOECSpec0 <- stackParamOECSpec0$lambda
                    
                }else if(stackCVSpec == "standard"){
                    
                    # use from above (standard stacking) -- tuned above in standard stacking section
                    stackParamOECSpec0 <- stackParamSpec0 # zero out
                    
                }else{
                    
                    #otherwise use from above generalist OEC 
                    stackParamOECSpec0 <- stackParamOEC
                }
                
                
                if(etaTuneInd){
                    
                    OECeta <-   oecEta_CV(data = full,
                                          tune.grid = etas,
                                          sslLambdas = tuneParam$lambda,
                                          stackParam = stackParamOEC,
                                          method = "glmnet",
                                          nfolds = nFolds, 
                                          nnlsInd = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          cv = etaTune,
                                          standardize = TRUE,
                                          weights = NULL,
                                          glmnetOpt = glmnetOpt,
                                          xStandardize = FALSE,
                                          AvgW = AvgWEtaTn,
                                          horizon = horizon)
                    
                }
                
                if(etaTuneIndSpec != FALSE){
                    
                    OECetaSpec <-   oecEta_CV(data = full,
                                              tune.grid = etas,
                                              sslLambdas = tuneParam$lambda,
                                              stackParam = stackParamOECSpec,
                                              method = "glmnet",
                                              nfolds = nFolds,
                                              nnlsInd = TRUE,
                                              sampSzWeight = sampSzWeight,
                                              cv = etaTuneIndSpec, 
                                              standardize = TRUE,
                                              weights = NULL,
                                              glmnetOpt = glmnetOpt, 
                                              xStandardize = FALSE,
                                              AvgW = AvgWEtaTn,
                                              SpecID = test_country,
                                              horizon = horizon)
                    
                }else{
                    OECetaSpec <- OECeta
                }
                
                
                if(etaTuneIndSpec0 != FALSE){
                    
                    OECetaSpec0 <-   oecEta_CV(data = full,
                                               tune.grid = etas,
                                               sslLambdas = tuneParamOEC$lambda,
                                               stackParam = stackParamOECSpec0,
                                               method = "glmnet",
                                               nfolds = nFolds,
                                               nnlsInd = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               cv = etaTuneIndSpec0, 
                                               standardize = TRUE,
                                               weights = NULL,
                                               glmnetOpt = glmnetOpt, 
                                               xStandardize = FALSE,
                                               AvgW = AvgW,
                                               SpecID = test_country,
                                               horizon = horizon
                    )
                    
                }else{
                    OECetaSpec0 <- OECeta
                }

            if(tuneParamZero)   tuneParam$lambda <- rep(0, K)
            
            
            warmRes <- ridgeWS(data = full, 
                                tuneParam = tuneParam$lambda, 
                                stackParam = stackParamOEC, 
                                nnlsInd = TRUE,
                                stackTune = TRUE, 
                                modelStandardize = TRUE,
                                stackStandardize = TRUE,
                                glmnetOpt = glmnetOpt, # used closed form
                                xStandardize = FALSE,
                                sampSzWeight = 1, 
                                weights = Diagonal( nrow(full) ),
                                lambdaGrid = lambdaVec,
                                AvgW = FALSE) 
            
            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack
           
            
                resMat[iterNum, 22] <- OECeta
                
                if(lowLimTune){
                    psiL <- oecLow_CV(data = full,
                                      tune.grid = psiVec,
                                      sslLambdas = tuneParam$lambda,
                                      stackParam = stackParamOEC,
                                      oecEta = OECeta,
                                      method = "glmnet",
                                      nfolds = nFolds, 
                                      nnlsInd = TRUE,
                                      sampSzWeight = sampSzWeight,
                                      cv = etaTune,
                                      standardize = TRUE,
                                      weights = NULL,
                                      glmnetOpt = glmnetOpt, 
                                      xStandardize = FALSE,
                                      AvgW = AvgWEtaTn)
                    
                }
                
                ### OEC
                if(sampSzWeight < 6){
                    pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                    for(it in 1:itrs){
                        set.seed(it)
                        
                        oecRidgeWS <- ridgeAltFix(data = full, 
                                                  betaStart = warmRes$beta, 
                                                  wStart = warmRes$w, 
                                                  lambdaVec = tuneParam$lambda, 
                                                  mu = stackParamOEC, 
                                                  nnlsInd = psiL / K,
                                                  low = psiL / K,
                                                  up = Inf,
                                                  tol = 0.001,
                                                  objCriteria = TRUE,
                                                  eta = OECeta,
                                                  dataSplit = 1,
                                                  projs = 0,
                                                  Avg = AvgW,
                                                  wUpdate = "glmnet",
                                                  sigK = sigK,
                                                  sigStack = sigStack,
                                                  standardize = TRUE,
                                                  sampSzWeight = sampSzWeight,
                                                  weights = Diagonal( nrow(full)),
                                                  simplex = simplexInd,
                                                  stackTol = stackTol)
                        
                        pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                                  mod = oecRidgeWS)
                        
                    }
                    
                    predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                    
                    resMat[iterNum, 10] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                }
                    
                    ############################ # update w 's only after all K beta_k have been updated
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    
                    oecRidgeWS <- ridgeAlt(data = full, 
                                              betaStart = warmRes$beta, 
                                              wStart = warmRes$w, 
                                              lambdaVec = tuneParam$lambda, 
                                              mu = stackParamOEC, 
                                              nnlsInd = TRUE,
                                              tol = 0.001,
                                              objCriteria = TRUE,
                                              eta = OECeta,
                                              dataSplit = 1,
                                              projs = 0,
                                              Avg = AvgW,
                                              wUpdate = "glmnet",
                                              sigK = sigK,
                                              sigStack = sigStack,
                                              standardize = TRUE,
                                              sampSzWeight = sampSzWeight,
                                              weights = Diagonal( nrow(full))
                                              )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                
                resMat[iterNum, 11] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        
                #################################################################################################### 
                ##########
                # left generalist
                ##########
                etaVec <- c( etas[etas <= OECeta] ) # last eta is the OECetaSpec
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on generalist stacking
                    if(it == 1){
                        B_itr <- betas
                        w_itr <- warmRes$w # generalist stacking weights
                        
                        set.seed(1)
                        # run initialization from left using OECeta as eta
                        oecRidgeWS <- ridgeAlt(data = full, 
                                               betaStart = B_itr, 
                                               wStart = w_itr, 
                                               lambdaVec = tuneParamOEC$lambda, 
                                               mu = stackParamOEC, 
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               orderRandom = orderRandom)
                        
                        # use this as just initializing on stacking
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 64] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 65] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- ridgeAlt(data = full, 
                                           betaStart = B_itr, 
                                           wStart = w_itr, 
                                           lambdaVec = tuneParamOEC$lambda, 
                                           mu = stackParamOEC, 
                                           nnlsInd = TRUE,
                                           tol = 0.001,
                                           objCriteria = TRUE,
                                           eta = eta_itr,
                                           dataSplit = 1,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glmnet",
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           weights = Diagonal( nrow(full)),
                                           orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 62] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 63] <- min(oecRidgeWS$objImp) # best objective achieved
                resMat[iterNum, 70] <- OECeta
                
                ##########
                # right
                ##########
                etaVec <- c( etas[etas >= OECeta] ) # last eta is the OECetaSpec
                etaVec <- sort(etaVec, decreasing = TRUE) # start big
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on specialist stacking
                    if(it == 1){
                        # initalize on all being equal to merged model (i.e., all betas are same with average weights)
                        B_itr <- replicate( ncol(betas), as.vector(mrgBeta) )
                        w_itr <- c(0, rep(1/ K, K) )
                        
                        # run as just initializing on country-specific weights using OECetaSpec as eta
                        set.seed(1)
                        oecRidgeWS <- ridgeAlt(data = full, 
                                               betaStart = B_itr, 
                                               wStart = w_itr, 
                                               lambdaVec = tuneParamOEC$lambda, 
                                               mu = stackParamOEC, 
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               orderRandom = orderRandom)
                        
                        # last run is the tuned eta
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 68] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 69] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- ridgeAlt(data = full, 
                                           betaStart = B_itr, 
                                           wStart = w_itr, 
                                           lambdaVec = tuneParamOEC$lambda, 
                                           mu = stackParamOEC, 
                                           nnlsInd = TRUE,
                                           tol = 0.001,
                                           objCriteria = TRUE,
                                           eta = eta_itr,
                                           dataSplit = 1,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glmnet",
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           weights = Diagonal( nrow(full)),
                                           orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 66] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 67] <- min(oecRidgeWS$objImp) # best objective achieved
                
                #################################################################################################### 
                    # ----------------------------
                    # ridgeAltFixNoIntercept
                oecRidgeWS <- ridgeAlt(data = full, 
                                       betaStart = warmRes$beta, 
                                       wStart = warmRes$w, 
                                       lambdaVec = tuneParam$lambda, 
                                       mu = stackParamOEC, 
                                       nnlsInd = TRUE,
                                       tol = 0.001,
                                       objCriteria = TRUE,
                                       eta = OECeta,
                                       dataSplit = 1,
                                       projs = 0,
                                       Avg = TRUE,
                                       wUpdate = "glmnet",
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       weights = Diagonal( nrow(full))
                )
                
                predsMat <- as.vector(oecRidgeWS$w[1]) + X %*% as.matrix( oecRidgeWS$beta )
                
                modTune <- cv.glmnet(y = as.vector(full$Y), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)
                
                if(stackCV == "zero")    stackParam <- 0
                
                mod <- glmnet(y = as.vector(full$Y), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% as.matrix( oecRidgeWS$beta ) %*% w[-1]
                resMat[iterNum, 15] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                
                
                ######################
                # Specialist OEC -- Country specific OEC -- initialize on country-specific stacking
                ######################
                # rerun Warm Start with different \mu (stackParam) for specialist
                warmRes <- ridgeWS(data = full, 
                                   tuneParam = tuneParam$lambda, 
                                   stackParam = stackParamOECSpec, 
                                   nnlsInd = TRUE,
                                   stackTune = TRUE, 
                                   modelStandardize = TRUE,
                                   stackStandardize = TRUE,
                                   glmnetOpt = glmnetOpt, 
                                   xStandardize = FALSE,
                                   sampSzWeight = 1, 
                                   weights = Diagonal( nrow(full) ),
                                   lambdaGrid = lambdaVec,
                                   AvgW = FALSE) 
                
                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack
                
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                                     betaStart = warmRes$beta, 
                                                     wStart = wCountry,
                                                     lambdaVec = tuneParam$lambda, 
                                                     mu = stackParamOECSpec, 
                                                     Stackindx = indx1, # rows corresponding to test country
                                                     nnlsInd = TRUE,
                                                     tol = 0.001,
                                                     objCriteria = TRUE,
                                                     eta = OECetaSpec,
                                                     dataSplit = 1,
                                                     projs = 0,
                                                     Avg = FALSE,
                                                     wUpdate = "glmnet",
                                                     sigK = sigK,
                                                     sigStack = sigStack,
                                                     standardize = TRUE,
                                                     sampSzWeight = sampSzWeight,
                                                     weights = Diagonal( nrow(full)),
                                                    low = psiL / K,
                                                    up = Inf,
                                                    stackTol = 1e-9
                                                    )
                
                        pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                                  mod = oecRidgeWS)
                        
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 16] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                
                ######################
                # Specialist OEC -- Country specific OEC-- initialize on country-specific stacking
                ######################
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wZeroOut,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 17] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )

####################################################################################################                
                ##########
                # left
                ##########
                etaVec <- c( etas[etas <= OECetaSpec] ) # last eta is the OECetaSpec
                
                indx1 <- which(full$Study == test_country)
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on specialist stacking
                    if(it == 1){
                        B_itr <- betas
                        w_itr <- wCountry
                        
                        set.seed(1)
                        
                        # run initialization from left using OECetaSpec as eta
                        oecRidgeWS <- ridgeAltSpec(data = full, 
                                                   betaStart = B_itr, 
                                                   wStart = w_itr,
                                                   lambdaVec = tuneParamOEC$lambda, 
                                                   mu = stackParamOECSpec, 
                                                   Stackindx = indx1, # rows corresponding to test country
                                                   nnlsInd = TRUE,
                                                   tol = 0.001,
                                                   objCriteria = TRUE,
                                                   eta = OECetaSpec,
                                                   dataSplit = 1,
                                                   projs = 0,
                                                   Avg = FALSE,
                                                   wUpdate = "glmnet",
                                                   sigK = sigK,
                                                   sigStack = sigStack,
                                                   standardize = TRUE,
                                                   sampSzWeight = sampSzWeight,
                                                   weights = Diagonal( nrow(full)),
                                                   low = psiL / K,
                                                   up = Inf,
                                                   stackTol = 1e-9,
                                                   orderRandom = orderRandom)
                        
                        # use this as just initializing on stacking
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 48] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 49] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = B_itr, 
                                               wStart = w_itr,
                                               lambdaVec = tuneParamOEC$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = eta_itr,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9,
                                               orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 46] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 47] <- min(oecRidgeWS$objImp) # best objective achieved
                
                ##########
                # right
                ##########
                etaVec <- c( etas[etas >= OECetaSpec] ) # last eta is the OECetaSpec
                etaVec <- sort(etaVec, decreasing = TRUE) # start big
                indx1 <- which(full$Study == test_country)
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on specialist stacking
                    if(it == 1){
                        # initalize on all being equal to country-specific model (i.e., all betas are same with average weights)
                        B_itr <- replicate( ncol(betas), as.vector(bCountry) )
                        w_itr <- c(0, rep(1/ K, K) )
                        
                        # run as just initializing on country-specific weights using OECetaSpec as eta
                        set.seed(1)
                        oecRidgeWS <- ridgeAltSpec(data = full, 
                                                   betaStart = B_itr, 
                                                   wStart = w_itr,
                                                   lambdaVec = tuneParamOEC$lambda, 
                                                   mu = stackParamOECSpec, 
                                                   Stackindx = indx1, # rows corresponding to test country
                                                   nnlsInd = TRUE,
                                                   tol = 0.001,
                                                   objCriteria = TRUE,
                                                   eta = OECetaSpec,
                                                   dataSplit = 1,
                                                   projs = 0,
                                                   Avg = FALSE,
                                                   wUpdate = "glmnet",
                                                   sigK = sigK,
                                                   sigStack = sigStack,
                                                   standardize = TRUE,
                                                   sampSzWeight = sampSzWeight,
                                                   weights = Diagonal( nrow(full)),
                                                   low = psiL / K,
                                                   up = Inf,
                                                   stackTol = 1e-9,
                                                   orderRandom = orderRandom)
                        
                        # last run is the tuned eta
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 52] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 53] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = B_itr, 
                                               wStart = w_itr,
                                               lambdaVec = tuneParamOEC$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = eta_itr,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9,
                                               orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 50] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 51] <- min(oecRidgeWS$objImp) # best objective achieved
                resMat[iterNum, 71] <- OECetaSpec
                ####################################################################################################      
                ######################
                # Specialist OEC -- Country specific OEC -- initialize on avg weights
                ######################
                indx1 <- which(full$Study == test_country)
                
                # intercept from country-specific avg weights
                w0 <- ( sum(full$Y[indx1]) - sum( cbind(1, as.matrix(full[indx1,-c(1,2)])) %*% 
                                                      betas %*% wCountry[-1]) ) / length(indx1)
                
                wA <- c(w0, rep(1 / K, K)) # average weights
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wA,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 18] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                ######################
                # Specialist OEC -- Country specific OEC -- initialize on avg weights with different intercept
                ######################
                indx1 <- which(full$Study == test_country)
                
                # intercept from country-specific avg weights
                w0 <- ( sum(full$Y[indx1]) - sum( cbind(1, as.matrix(full[indx1,-c(1,2)])) %*% 
                                                      betas %*% wZeroOut[-1]) ) / length(indx1)
                
                wA <- c(w0, rep(1 / K, K)) # average weights
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wA,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 19] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                ######################
                # Specialist OEC -- Country specific OEC -- initialize on avg weights with different intercept
                ######################
                indx1 <- which(full$Study == test_country)
                
                # intercept from country-specific avg weights
                w0 <- ( sum(full$Y[indx1]) - sum( cbind(1, as.matrix(full[indx1,-c(1,2)])) %*% 
                                                      betas %*% rep(1/K, K)) ) / length(indx1)
                
                wA <- c(w0, rep(1 / K, K)) # average weights
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wA,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 20] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                
                ######################
                # Specialist OEC -- Country specific OEC AVG Weights
                ######################
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = TRUE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 21] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                ###################
                # average weights
                ###################
                w0 <- ( sum(full$Y) - sum( cbind(1, as.matrix(full[,-c(1,2)])) %*% betas %*%  warmRes$w[-1]) ) / nrow(full)
                
                warmRes$w <- c(w0, rep(1 / K, K)) # average weights
                predsVecOEC <- oecRidgePred(data = X_test[,-1],
                                            mod = warmRes)
                #}
                resMat[iterNum, 14] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                rm(oecRidgeWS, predsVecOEC)
                
                #######################################################
                # average specialist and then do a final stacking stage
                #######################################################
                indx1 <- which(full$Study == test_country)
                
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = TRUE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                predsMat <- cbind(1, as.matrix(full[indx1,-c(1,2)] ) ) %*% as.matrix( oecRidgeWS$beta )
                
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam1 <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)
                
                if(stackCV == "zero")    stackParam1 <- 0
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam1, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% as.matrix( oecRidgeWS$beta ) %*% w[-1]
                resMat[iterNum, 29] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                
                ###################
                
                #######################################################
                # average OEC and then do a final stacking stage
                #######################################################
              
                oecRidgeWS <- ridgeAlt(data = full, 
                                       betaStart = warmRes$beta, 
                                       wStart = warmRes$w, 
                                       lambdaVec = tuneParam$lambda, 
                                       mu = stackParamOEC, 
                                       nnlsInd = TRUE,
                                       tol = 0.001,
                                       objCriteria = TRUE,
                                       eta = OECeta,
                                       dataSplit = 1,
                                       projs = 0,
                                       Avg = TRUE,
                                       wUpdate = "glmnet",
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       weights = Diagonal( nrow(full))
                )
                
                predsMat <-  cbind(1, as.matrix(full[,-c(1,2)] ) ) %*% as.matrix( oecRidgeWS$beta )
                
                modTune <- cv.glmnet(y = as.vector(full$Y), 
                                     x = as.matrix(predsMat),
                                     alpha = 0, 
                                     lambda = lambdaVec, 
                                     standardize = TRUE,
                                     intercept = TRUE,
                                     lower.limits = 0, 
                                     thresh = 1e-10
                ) 
                
                stackParam1 <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)
                
                if(stackCV == "zero")    stackParam1 <- 0
                
                mod <- glmnet(y = as.vector(full$Y), 
                              x = as.matrix(predsMat), 
                              alpha = 0,         
                              lambda = stackParam1, 
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = 0, 
                              thresh = 1e-10 ) 
                
                w <- as.vector( coef( mod , 
                                      exact = TRUE, 
                                      y = as.vector(full$Y), 
                                      x = as.matrix(predsMat)
                )
                )
                
                rm(mod)
                
                # standard stacking predictions
                preds <- w[1] + X_test %*% as.matrix( oecRidgeWS$beta ) %*% w[-1]
                resMat[iterNum, 30] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                
                
                ###################
                
                ######################
                # Specialist OEC -- Country specific Zero Out OEC-- initialize on country-specific stacking
                ######################
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec0(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wZeroOut,
                                               lambdaVec = tuneParam$lambda, 
                                               mu = stackParamOECSpec0, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec0,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = psiL / K,
                                               up = Inf,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 31] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )

####################################################################################################                
                ##########
                # left Zero Out
                ##########
                etaVec <- c( etas[etas <= OECetaSpec0] ) # last eta is the OECetaSpec
                
                indx1 <- which(full$Study == test_country)
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on zero out stacking
                    if(it == 1){
                        B_itr <- betas
                        w_itr <- wZeroOut
                        
                        set.seed(1)
                        # run initialization from left using OECetaSpec as eta
                        oecRidgeWS <- ridgeAltSpec0(data = full, 
                                                    betaStart = B_itr, 
                                                    wStart = w_itr,
                                                    lambdaVec = tuneParamOEC$lambda, 
                                                    mu = stackParamOECSpec0, 
                                                    Stackindx = indx1, # rows corresponding to test country
                                                    nnlsInd = TRUE,
                                                    tol = 0.001,
                                                    objCriteria = TRUE,
                                                    eta = OECetaSpec0,
                                                    dataSplit = 1,
                                                    projs = 0,
                                                    Avg = FALSE,
                                                    wUpdate = "glmnet",
                                                    sigK = sigK,
                                                    sigStack = sigStack,
                                                    standardize = TRUE,
                                                    sampSzWeight = sampSzWeight,
                                                    weights = Diagonal( nrow(full)),
                                                    low = psiL / K,
                                                    up = Inf,
                                                    stackTol = 1e-9,
                                                    orderRandom = orderRandom)
                        
                        # use this as just initializing on stacking
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 56] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 57] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec0(data = full, 
                                                betaStart = B_itr, 
                                                wStart = w_itr,
                                                lambdaVec = tuneParamOEC$lambda, 
                                                mu = stackParamOECSpec0, 
                                                Stackindx = indx1, # rows corresponding to test country
                                                nnlsInd = TRUE,
                                                tol = 0.001,
                                                objCriteria = TRUE,
                                                eta = eta_itr,
                                                dataSplit = 1,
                                                projs = 0,
                                                Avg = FALSE,
                                                wUpdate = "glmnet",
                                                sigK = sigK,
                                                sigStack = sigStack,
                                                standardize = TRUE,
                                                sampSzWeight = sampSzWeight,
                                                weights = Diagonal( nrow(full)),
                                                low = psiL / K,
                                                up = Inf,
                                                stackTol = 1e-9,
                                                orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 54] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 55] <- min(oecRidgeWS$objImp) # best objective achieved
                
                ##########
                # right
                ##########
                etaVec <- c( etas[etas >= OECetaSpec0] ) # last eta is the OECetaSpec
                etaVec <- sort(etaVec, decreasing = TRUE) # start big
                indx1 <- which(full$Study == test_country)
                
                for(it in 1:length(etaVec) ){
                    
                    # if first pass the initialize on specialist stacking
                    if(it == 1){
                        # initalize on all being equal to country-specific model (i.e., all betas are same with average weights)
                        B_itr <- replicate( ncol(betas), as.vector(bCountry) )
                        w_itr <- c(0, rep(1/ K, K) )
                        
                        # run as just initializing on country-specific weights using OECetaSpec as eta
                        set.seed(1)
                        oecRidgeWS <- oecRidgeWS <- ridgeAltSpec0(data = full, 
                                                                  betaStart = B_itr, 
                                                                  wStart = w_itr,
                                                                  lambdaVec = tuneParamOEC$lambda, 
                                                                  mu = stackParamOECSpec0, 
                                                                  Stackindx = indx1, # rows corresponding to test country
                                                                  nnlsInd = TRUE,
                                                                  tol = 0.001,
                                                                  objCriteria = TRUE,
                                                                  eta = OECetaSpec0,
                                                                  dataSplit = 1,
                                                                  projs = 0,
                                                                  Avg = FALSE,
                                                                  wUpdate = "glmnet",
                                                                  sigK = sigK,
                                                                  sigStack = sigStack,
                                                                  standardize = TRUE,
                                                                  sampSzWeight = sampSzWeight,
                                                                  weights = Diagonal( nrow(full)),
                                                                  low = psiL / K,
                                                                  up = Inf,
                                                                  stackTol = 1e-9,
                                                                  orderRandom = orderRandom)
                        
                        # last run is the tuned eta
                        predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                                     mod = oecRidgeWS)
                        
                        resMat[iterNum, 60] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                        resMat[iterNum, 61] <- min(oecRidgeWS$objImp) # best objective achieved
                    }
                    
                    eta_itr <- etaVec[it]
                    set.seed(1)
                    oecRidgeWS <- oecRidgeWS <- ridgeAltSpec0(data = full, 
                                                              betaStart = B_itr, 
                                                              wStart = w_itr,
                                                              lambdaVec = tuneParamOEC$lambda, 
                                                              mu = stackParamOECSpec0, 
                                                              Stackindx = indx1, # rows corresponding to test country
                                                              nnlsInd = TRUE,
                                                              tol = 0.001,
                                                              objCriteria = TRUE,
                                                              eta = eta_itr,
                                                              dataSplit = 1,
                                                              projs = 0,
                                                              Avg = FALSE,
                                                              wUpdate = "glmnet",
                                                              sigK = sigK,
                                                              sigStack = sigStack,
                                                              standardize = TRUE,
                                                              sampSzWeight = sampSzWeight,
                                                              weights = Diagonal( nrow(full)),
                                                              low = psiL / K,
                                                              up = Inf,
                                                              stackTol = 1e-9,
                                                              orderRandom = orderRandom)
                    
                    # update betas and w
                    B_itr <- oecRidgeWS$beta
                    w_itr <- oecRidgeWS$w
                    
                    
                }
                
                # last run is the tuned eta
                predsVecOEC1 <- oecRidgePred(data = X_test[,-1],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 58] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 59] <- min(oecRidgeWS$objImp) # best objective achieved
                resMat[iterNum, 72] <- OECetaSpec0
####################################################################################################    
                
                ######################
                # Specialist OEC -- Country specific Zero Out OEC-- initialize on country-specific stacking
                ######################
                # constrain from above
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                                betaStart = warmRes$beta, 
                                                wStart = wZeroOut,
                                                lambdaVec = tuneParam$lambda, 
                                                mu = stackParamOECSpec, 
                                                Stackindx = indx1, # rows corresponding to test country
                                                nnlsInd = TRUE,
                                                tol = 0.001,
                                                objCriteria = TRUE,
                                                eta = OECetaSpec,
                                                dataSplit = 1,
                                                projs = 0,
                                                Avg = FALSE,
                                                wUpdate = "glmnet",
                                                sigK = sigK,
                                                sigStack = sigStack,
                                                standardize = TRUE,
                                                sampSzWeight = sampSzWeight,
                                                weights = Diagonal( nrow(full)),
                                                low = psiL / K,
                                                up = psiH,
                                                stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 32] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                
                
                ####################
                # Upper Limits Specialist 
                ####################
                
                psiH <- oecHigh_CV(data = full,
                                   tune.grid = seq(1 / K ,1, length = 20), 
                                   sslLambdas = tuneParamOEC$lambda,
                                   stackParam = stackParamOEC,
                                   oecEta = OECeta,
                                   method = "glmnet",
                                   nfolds = nFolds, 
                                   nnlsInd = TRUE,
                                   sampSzWeight = sampSzWeight,
                                   cv = "cv",
                                   SpecID = test_country,
                                   standardize = TRUE,
                                   weights = NULL,
                                   glmnetOpt = glmnetOpt, 
                                   xStandardize = FALSE,
                                   AvgW = AvgW)
                
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpec(data = full, 
                                               betaStart = warmRes$beta, 
                                               wStart = wZeroOut,
                                               lambdaVec = tuneParamOEC$lambda, 
                                               mu = stackParamOECSpec0, 
                                               Stackindx = indx1, # rows corresponding to test country
                                               nnlsInd = TRUE,
                                               tol = 0.001,
                                               objCriteria = TRUE,
                                               eta = OECetaSpec0,
                                               dataSplit = 1,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glmnet",
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               standardize = TRUE,
                                               xStandardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               weights = Diagonal( nrow(full)),
                                               low = 0,
                                               up = psiH,
                                               stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                #}
                resMat[iterNum, 41] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                resMat[iterNum, 43] <- psiH
                
                
                ####################
                # Window Specialist 
                ####################
                
                psiH <- oecHigh_CV(data = full,
                                   tune.grid = seq(0, 1, length = 20),
                                   sslLambdas = tuneParamOEC$lambda,
                                   stackParam = stackParamOEC,
                                   oecEta = OECeta,
                                   method = "glmnet",
                                   nfolds = nFolds, 
                                   nnlsInd = TRUE,
                                   sampSzWeight = sampSzWeight,
                                   cv = "cvWindow",
                                   SpecID = test_country,
                                   standardize = TRUE,
                                   weights = NULL,
                                   glmnetOpt = glmnetOpt, 
                                   xStandardize = FALSE,
                                   AvgW = AvgW)
                
                indx1 <- which(full$Study == test_country)
                
                pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
                for(it in 1:itrs){
                    set.seed(it)
                    oecRidgeWS <- ridgeAltSpecWindow(data = full, 
                                                     betaStart = warmRes$beta, 
                                                     wStart = wZeroOut,
                                                     lambdaVec = tuneParamOEC$lambda, 
                                                     mu = 0,
                                                     Stackindx = indx1, # rows corresponding to test country
                                                     nnlsInd = TRUE,
                                                     tol = 0.001,
                                                     objCriteria = TRUE,
                                                     eta = OECetaSpec0,
                                                     dataSplit = 1,
                                                     projs = 0,
                                                     Avg = FALSE,
                                                     wUpdate = "glmnet",
                                                     sigK = sigK,
                                                     sigStack = sigStack,
                                                     standardize = TRUE,
                                                     sampSzWeight = sampSzWeight,
                                                     weights = Diagonal( nrow(full)),
                                                     low = 0,
                                                     up = psiH,
                                                     stackTol = 1e-9
                    )
                    
                    pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                              mod = oecRidgeWS)
                    
                }
                
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
                resMat[iterNum, 42] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                resMat[iterNum, 44] <- psiH
                
                # -----------------------------------
                resMat[iterNum, 33] <- OECeta
                resMat[iterNum, 34] <- OECetaSpec
                resMat[iterNum, 35] <- OECetaSpec0
                resMat[iterNum, 36] <- stackParamOEC
                resMat[iterNum, 37] <- stackParamOECSpec
                resMat[iterNum, 38] <- stackParamOECSpec0
                
                
                rm(warmRes, oecRidgeWS, preds, predsMat)
                
                print(paste0("complete ", taskNm))
                print(resMat[iterNum,])
            }
        }
        
        # remove rows with only NAs (except K the last column which is filled in even when no models fit)
        ind <- apply(resMat[,-ncol(resMat)], 1, function(x) !all(is.na(x)))
        
        resMat <- resMat[ind,]
        
        
    
    return(resMat)
        
}

timeEnd <- Sys.time()
print(difftime(timeEnd, timeStart, units='mins'))

resMat <- do.call(rbind, results)
colnames(resMat) <- c("merge", "avg", "stacking", "stacking_zeroOut", 
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
                      "oec_SpecUp", "oecWindow", "up", "upWindow", "country_noLinear",
                      "oec_specAnneal_left", "oec_specAnneal_leftObj", "oec_spec_left", "oec_spec_leftObj",
                      "oec_specAnneal_right", "oec_specAnneal_rightObj", "oec_spec_right", "oec_spec_rightObj",
                      "oec_0Anneal_left", "oec_0Anneal_leftObj", "oec_0_left", "oec_0_leftObj",
                      "oec_0Anneal_right", "oec_0Anneal_rightObj", "oec_0_right", "oec_0_rightObj",
                      "oec_genAnneal_left", "oec_genAnneal_leftObj", "oec_gen_left", "oec_gen_leftObj",
                      "oec_genAnneal_right", "oec_genAnneal_rightObj", "oec_gen_right", "oec_gen_rightObj",
                      "etaGen", "etaSpec", "etaSpec0")


print("save file")
setwd(save.folder)
write.csv(resMat, fileNm)
