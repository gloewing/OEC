# Gabe Loewinger
# August 10, 2020
# Multi Study Mortality COVID

source("msMort_functions.R")
source("Study Strap Functions.R")
source("OEC Functions.R")

pM <- expand.grid( c(FALSE), #nH
             c(100,150), # minTrain
             c(1,4), # xPro
             c(1)) # sampSzWt

itr <- 2

# test set date ranges
month_start <- "-03-01"  # "-01-01" #
month_end <- "-02-28" # "-12-31" #
addYear <- TRUE #FALSE # addyear means add a year to the test date--if start and end are the same year then set to FALSE
testYears <- 2011:2019 # 2010:2019 #chosen because many countries in March 2012 have just over 2 years worth of data, and
# also allows us to go up until March 2020, right as pandemic was beginning
# 2019 goes up to 2020 because the test year starts in March

minTrain <- pM[itr, 2] # minimum number of training data points to be considered a training country for any iterations
minTest <- 50 # minimum number of testing data points to be considered a test country for any iterations
testDate_start <- "-03-01" #"-01-01" #"-03-01" 
testDate_end <- "-02-28"  #"-12-31" #"-02-28"
testCountryData <- 12 # this is the number of months that we pretend the test country has (different than minTest)
nH <- pM[itr, 1] # northern hemisphere
# to emulate south africa

full <- read.csv("world-mort.csv")
full$date <- as.Date(full$date)

if(nH){
    # northern hemisphere countries
    countries <- unique(full$country)
    countries <- setdiff(countries, c("Ecuador", "Peru", "Chile", "South Africa"))
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
oecTune <- "cvCFTS" #"cvOEC" # "cvCFTS"#"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
# added in
stdTn <- "cvCFTS" #"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
############
stackCV <- "cv" #"cv" # can be "1", "cv", "hoso" or "zero." If 1, then just use the same parameter as in standard vanilla stacking tuned above. If "zero", then set stacking hyperparamter to 0 if using window stacking
stackCVSpec <- FALSE #"cvSpec"# "cvSpecTS" # "cvSpecTSHOO" # FALSE # "cvSpecTS" # 
stackCVSpec0 <-  FALSE #"cvSpec0" #"cvSpecTS0" # "cvSpec0TSHOO"
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
xProcess <- pM[itr, 3]
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

knots <- 15

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
fileNm <- paste0("mrtTmA_eta", OECeta, 
                 #"_spInt_", intercept,
                  "_xPro_", xProcess, 
                 # "_psi_", 
                 # psiL, ".", psiH,
                 #"_Avg_", AvgW, 
                 "_etaTn_",  
                 # etaTuneInd, "_", 
                 etaTune,
                 "_etaSpec_", etaTuneIndSpec,
                 # "_etaSpec0_", etaTuneIndSpec0,
                 # "_eTnAvgW_", AvgWEtaTn, 
                 "_Wcv_", stackCV, 
                 "_Wspec_", stackCVSpec,
                  "_Wspec0_", stackCVSpec0,
                 "_TScv_", TScv,
                 "_sclX_", scaleInd,
                 #"_alph0_", tuneParamZero,
                 "_glm_", glmnetOpt,
                 "_hr_", horizon,
                 # ".pca.",pcaInd,
                 "_oecTn_", oecTune, 
                 "_stdTn_", stdTn,
                 "_tnOr_", tnOrder,
                 # "_itr_", itrs,
                 "_tstCntMn_", testCountryData,
                 ".nH", nH,
                 # "_lwLmTn_", lowLimTune,
                 "_mnSt_", month_start,
                 "_mnTr_", minTrain,
                 "_smpSzWt_", sampSzWeight)
print(fileNm)

######### CHANGE NUMBER OF THREADS BEFORE RUNNING REAL THING !!!!!!!!!!!!!!!!!!!!!!!
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
    
    for(iterNum in 1:length(testYears) ){

        full <- read.csv("world-mort.csv")
        full$date <- as.Date(full$date)
        
        
        if(nH){
            # northern hemisphere countries
            countries <- unique(full$country)
            countries <- setdiff(countries, c("Ecuador", "Peru", "Chile", "South Africa"))
        }else{
            countries <- unique(full$country)
        }
        
        
        K <- length(countries)
        
        taskNm <- paste0("country_", countries[cnt], "_year: ", iterNum)
        print(paste0("begin ", taskNm))
        # month_start <- "-03-01" 
        # month_end <- "-02-29"
        # testYears <- 2012:2020 # chosen because many countries in March 2012 have just over 2 years worth of data, and
        # # also allows us to go up until March 2020, right as pandemic was beginning
        # 
        # testDate_start <- "-03-01" 
        # testDate_end <- "-02-29"
        
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
            
            # countries with training data for this specific test year
            #countriesTrain <- unique( full$country[train_indx]  ) # training data 
            #countriesTest <- unique( full$country[T_indx]  ) # testing data
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
            full$country <- as.numeric(full$country) # turn numeric
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
            # generate spline design matrix
            # X  <- as.matrix( ns(as.numeric(full$date), df = knots, intercept = TRUE) )
            # X_test <- as.matrix( ns(as.numeric(test$date), df = knots, intercept = TRUE) )
            
            if(xProcess == 1){
                
                # process each country and the test set separately
                
                ### xProcess == 1
                cList <- infoList <- vector("list", length = K)
                # Yvec <- c() #vector(length = nrow(full)) # store Y
                
                for(k in 1:K){
                    # process each country separately
                    
                    # country specific rows
                    indx <- which(full$country == countries[k])
                    full2 <- full[indx,]
                    full2 <- full2[order(full2$date),] # reorder by date
                    colnames(full2)[nm] <- "outcome"
                    # Yvec <- c(Yvec, full2$outcome)
                    
                    infoList[[k]] <- full2
                    
                    cList[[k]] <- Xgen(full2)
                    rm(full2)
                    
                    # make design matrices match-- extremely ad hoc!!! Change
                    # if(ncol(cList[[k]]) > 6){
                    #     cList[[k]] <- cList[[k]][,-3]
                    # }
        
                }
                
                # concatenate countries
                X <- do.call(rbind, cList)
                rm(cList)
                
                # update these since they are reordered
                full <- do.call(rbind, infoList)
                colnames(full)[nm] <- "Y"
                rm(infoList)
                
                # test2 <- test
                # colnames(test2)[nm] <- "outcome"
                # test2 <- test2[order(test2$date),]
                # X_test = Xgen(test2)
                # rm(test2)

                colnames(test)[nm] <- "outcome"
                test <- test[order(test$date),]
                X_test <- Xgen(test)
                colnames(test)[nm] <- "Y"
                
                # if(ncol(X_test) > 6){
                #     X_test <- X_test[,-3]
                # }
                
            }else if(xProcess == 2){
                # process test set and all countries together as one universal design matrix
                
                test$country <- test_country # make this into numeric
                
                # order full and test by dates
                full <- rbind(full, test) # recombine into one
                full <- full[order(full$date),] # reorder by date
                colnames(full)[nm] <- "outcome"
                X <- Xgen(full)
                colnames(full)[nm] <- "Y"
                
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
            }else if(xProcess == 3){
                
                # process each country separately but include the test set (of the test country) in each one
                test$country <- test_country # make this into numeric
                
                ### xProcess == 3
                cList <- infoList <- vector("list", length = K)
                # dateVec <- Yvec <- countryVec <- c() #vector(length = nrow(full)) # store Y, country and date
                
                for(k in 1:K){
                    # process each country separately
                    
                    # country specific rows
                    indx <- which(full$country == countries[k])
                    full2 <- rbind(full[indx,], test) # concatenate test set to each country-specific design matrix
                    full2 <- full2[order(full2$date),] # reorder by date
                    
                    # which rows correspond to test set after ordering by date
                    # specific country between certain dates
                    Countindx <- which(full2$country == test_country) # test country
                    CindxDates <- which( full2$date <= testDate_end & full2$date >= testDate_start) # test dates
                    Countindx <- intersect(Countindx, CindxDates) # test country during test dates
                    
                    colnames(full2)[nm] <- "outcome"
                    # Yvec <- c(Yvec, full2$outcome[-Countindx]) # remove rows corresponding to test set
                    # dateVec <- c(dateVec, full2$date[-Countindx]) # remove rows corresponding to test set
                    # countryVec <- c(countryVec, full2$country[-Countindx]) # remove rows corresponding to test set
                    infoList[[k]] <- full2[-Countindx, ] # save matrix of everything now that it is reordered, but remove country-specific stuff
                    
                    if(k == test_country){
                        
                        X_tmp <-  Xgen(full2)
                        cList[[k]] <- X_tmp[-Countindx,] # training set for test country
                        X_test <- X_tmp[Countindx,]  # test set for test country
                        # test$Y <- full2$outcome[Countindx] # outcome for test set
                        # test$date <- full2$date[Countindx] # date for test set
                        # test$population <- full2$population[Countindx] # population
                        test <- full2[Countindx, ] # this is reordered so save it in new order
                        colnames(test)[nm] <- "Y" # needed for for loop
                        
                        rm(X_tmp)
                        
                    }else{
                        
                        cList[[k]] <- Xgen(full2)[-Countindx,] # remove rows corresponding to test set
                        
                    }
                    
                    rm(full2)
                    
                }
                
                # concatenate countries
                X <- do.call(rbind, cList)
                
                # update these since they are reordered
                full <- do.call(rbind, infoList)
                
                colnames(full)[nm] <- "Y" # rename outcome
                # full$Y <- Yvec 
                # full$date <- dateVec
                # full$country <- countryVec
                
                rm(cList, Countindx, infoList)
                #rm(Yvec, dateVec, countryVec, )
            }else if(xProcess == 4){
                # Linear Trend
                # process each country and the test set separately 
                
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
                
        }
            
            # convert outcome to scaled percentage
            full$Y  <- full$Y / full$population * 100000
            test$Y  <- test$Y / test$population * 100000
            
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
            
            if(pcaInd){
                pcs <- prcomp(X)
                # put last PC first so intercept (column of ones) is still first
                pcs$rotation <- cbind(pcs$rotation[, ncol(pcs$rotation)], pcs$rotation[, -ncol(pcs$rotation)])
                X <- cbind(1, pcs$x[,-ncol(pcs$x)] ) # remove last column because just 0s
                X_test <- X_test %*% pcs$rotation
                rm(pcs)
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
            mergedLambda <- hosoCV(data = full,
                                   tune.grid,
                                   hoso = "merged",
                                   method = "glmnet",
                                   metric = "RMSE",
                                   nfolds = "K",
                                   nnlsInd = TRUE,
                                   OECeta = 0.5,
                                   sampSzWeight = sampSzWeight,
                                   weights = NULL)
            # OEC Tune
            set.seed(cnt) # added 9 /3/20
            tuneParam <- hosoCV(data = full,
                                tune.grid,
                                hoso = stdTn,
                                method = "glmnet",
                                metric = "RMSE",
                                nfolds = "K",
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
                                       nfolds = 5,
                                       nnlsInd = TRUE,
                                       OECeta = 0.5,
                                       sampSzWeight = sampSzWeight,
                                       weights = NULL)
                
                # # OEC Tune 0
                # set.seed(cnt) # added 9 /3/20
                # tuneParamOEC0 <- hosoCV(data = full,
                #                        tune.grid,
                #                        hoso = oecTune,
                #                        method = "glmnet",
                #                        lambdaVec = tuneParamOEC$lambda,
                #                        metric = "RMSE",
                #                        nfolds = 5,
                #                        nnlsInd = TRUE,
                #                        OECeta = 0.5,
                #                        sampSzWeight = sampSzWeight,
                #                        weights = NULL)
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
            
            # fit <- lm(full$Y ~ X - 1,
            #           model = FALSE)#)
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
                
                # indx <- which(full$Study == countries[k])
                indx <- which(full$Study == k) # changed for indices 9/11/20
                
                # country specific model
                # fit <- lm(full$Y[indx] ~ X[indx,] - 1,
                #           model = FALSE)
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
        
        
                # remove NAs from coefficient estimates if there are any
                # b <- as.numeric( coef(fit) )
                # b <- ifelse(is.na(b), 0 , b)
                
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
                    #b <- fitTest # coefficients
                    bCountry <- b
                    preds <- X_test %*% b
                    resMat[iterNum, 9] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
                    # fitTest <- b
                    
                    
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
            if(TScv){
                
                
                mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y, predsMat),
                                      tune.grid = tune.grid,
                                      nnls = TRUE,
                                      nfolds = 10)
                    
                stackParamOEC <- stackParam <- mod_stackCVTS$lambda # use same for OEC (incase its not retuned below for OEC)
                
                
            }else{
                
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
                
            }

            
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
            if(TScv){
                
                
                mod_mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y, predsMat),
                                       tune.grid = tune.grid,
                                       nnls = TRUE)
                
                stackParam <- mod_stackCVTS$lambda
                
            }else{
                
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
                
            }
            
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
                if(TScv){
                    
                    
                    mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat),
                                           tune.grid = tune.grid,
                                           nnls = TRUE)
                    
                    stackParam <- mod_stackCVTS$lambda
                    
                }else{
                    
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
                    
                }
                
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
                if(TScv){
                    
                    
                    mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat),
                                           tune.grid = tune.grid,
                                           nnls = TRUE)
                    
                    stackParam <- mod_stackCVTS$lambda
                    
                }else{
                    
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
                    
                }
                
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
            
            if(TScv){
                
                
                mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat_testCountry),
                                       tune.grid = tune.grid,
                                       nnls = TRUE)
                
                stackParam <- mod_stackCVTS$lambda
                
            }else{
                
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
                
            }
            
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
            
            if(TScv){
                
                
                mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat_testCountry2),
                                       tune.grid = tune.grid,
                                       nnls = TRUE)
                
                stackParam <- mod_stackCVTS$lambda
                
            }else{
                
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
                
            }
            
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
                
                if(TScv){
                    
                    
                    mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat),
                                           tune.grid = tune.grid,
                                           nnls = TRUE)
                    
                    stackParam <- mod_stackCVTS$lambda
                    
                }else{
                    
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
                    
                }
                
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
                if(TScv){
                    
                    
                    mod_stackCVTS <- stackCVTS(data = data.frame(Y = full$Y[indx1], predsMat),
                                           tune.grid = tune.grid,
                                           nnls = TRUE)
                    
                    stackParam <- mod_stackCVTS$lambda
                    
                }else{
                    
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
                    
                }
                
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
            # tuneParam <- rep(0, K) # no study-specific regularization
            
            # needed for OEC
            # 
            ##################
            # Order study labels
            ##################
            # rename studies from 1:K 
            # colnames(full)[1] <- "Study"
            # full$Study <- as.numeric(full$Study) # turn numeric
            # studyVec <- studies <- sort( unique(full$Study) )
            # num.trainStudy <- length(studyVec)
            # 
            # for(j in 1:K){
            #     
            #     indx <- which(full$Study == studyVec[j])
            #     full$Study[indx] <- j
            #     
            # }
            # 
            # rm(indx, studyVec)
            # 
            # full <- as.data.frame( cbind(full$Study, full$Y, X) )
            # if(intercept)           full <- full[,-3] # remove column of ones
            # colnames(full) <- c("Study", "Y", 1:(ncol(full) - 2))
            # 
            
            ##########################################
            # Tune OEC stacking hyperparameter (for w)
            ##########################################
            # if stackCV == "1", then just use the value above
            # if "zero" then just no penalization for OEC stacking regression
            # if "cv" or "hoso" then tune it in the OEC framework
            ##########################################
            # Tune OEC stacking hyperparameter (for w)
            ##########################################
            # if stackCV == "1", then just use the value above
            # if "zero" then just no penalization for OEC stacking regression
            # if "cv" or "hoso" then tune it in the OEC framework
            
            if(tnOrder){
                
                
                if(etaTuneInd){
                    
                    OECeta <-   oecEta_CV(data = full,
                                          tune.grid = etas,
                                          sslLambdas = tuneParam$lambda,
                                          stackParam = stackParamOEC1,
                                          method = "glmnet",
                                          nfolds = "K",
                                          nnlsInd = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          cv = etaTune,
                                          standardize = TRUE,
                                          weights = NULL,
                                          glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
                                          xStandardize = FALSE,
                                          AvgW = AvgWEtaTn,
                                          horizon = horizon)
                    
                }
                
                if(etaTuneIndSpec != FALSE){
                    
                    OECetaSpec <-   oecEta_CV(data = full,
                                              tune.grid = etas,
                                              sslLambdas = tuneParam$lambda,
                                              stackParam = as.numeric( resMat[iterNum, 39] ),
                                              method = "glmnet",
                                              nfolds = "K",
                                              nnlsInd = TRUE,
                                              sampSzWeight = sampSzWeight,
                                              cv = etaTuneIndSpec, #"cvSpecTS",
                                              standardize = TRUE,
                                              weights = NULL,
                                              glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
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
                                               stackParam = as.numeric( resMat[iterNum, 40] ),
                                               method = "glmnet",
                                               nfolds = "K",
                                               nnlsInd = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               cv = etaTuneIndSpec0, #"cvSpecTS0",
                                               standardize = TRUE,
                                               weights = NULL,
                                               glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
                                               xStandardize = FALSE,
                                               AvgW = AvgW,
                                               SpecID = test_country,
                                               horizon = horizon
                    )
                    
                }else{
                    OECetaSpec0 <- OECeta
                }
                
                if(stackCV == "zero"){
                    
                    stackParamOEC <- 0 # zero out
                    # if(psiL + psiH > 0){
                    #     # if either are nonzero 
                    # }
                }else if(stackCV == "cv" | stackCV == "hoso"){
                    ############################################
                    # Tune Stacking Hyperparameter again with tuned eta
                    ############################################
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOEC <- oecW_CV(data = full,
                                             tune.grid = tune.grid,
                                             sslLambdas = tuneParamOEC$lambda,
                                             method = "glmnet",
                                             nfolds = K,
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
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOECSpec <- oecW_CV(data = full,
                                                 tune.grid = tune.grid,
                                                 sslLambdas = tuneParamOEC$lambda,
                                                 method = "glmnet",
                                                 nfolds = K,
                                                 nnlsInd = TRUE,
                                                 OECeta = OECetaSpec,
                                                 sampSzWeight = sampSzWeight,
                                                 glmnetOpt = glmnetOpt,
                                                 xStandardize = FALSE, # I just chose this as FALSE 9/2/20
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
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOECSpec0 <- oecW_CV(data = full,
                                                  tune.grid = tune.grid,
                                                  sslLambdas = tuneParamOEC$lambda,
                                                  method = "glmnet",
                                                  nfolds = K,
                                                  nnlsInd = TRUE,
                                                  OECeta = OECetaSpec0,
                                                  sampSzWeight = sampSzWeight,
                                                  glmnetOpt = glmnetOpt,
                                                  xStandardize = FALSE, # I just chose this as FALSE 9/2/20
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
                
                
                
            }else{
                # dont do tune order
                
                
                if(stackCV == "zero"){
                    
                    stackParamOEC <- 0 # zero out
                    # if(psiL + psiH > 0){
                    #     # if either are nonzero 
                    # }
                }else if(stackCV == "cv" | stackCV == "hoso"){
                    ############################################
                    # Tune Stacking Hyperparameter again with tuned eta
                    ############################################
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOEC <- oecW_CV(data = full,
                                             tune.grid = tune.grid,
                                             sslLambdas = tuneParamOEC$lambda,
                                             method = "glmnet",
                                             nfolds = K,
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
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOECSpec <- oecW_CV(data = full,
                                                 tune.grid = tune.grid,
                                                 sslLambdas = tuneParamOEC$lambda,
                                                 method = "glmnet",
                                                 nfolds = K,
                                                 nnlsInd = TRUE,
                                                 OECeta = OECeta,
                                                 sampSzWeight = sampSzWeight,
                                                 glmnetOpt = glmnetOpt,
                                                 xStandardize = FALSE, # I just chose this as FALSE 9/2/20
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
                    # consider switching order to be before tuneParam in reTuning
                    stackParamOECSpec0 <- oecW_CV(data = full,
                                                  tune.grid = tune.grid,
                                                  sslLambdas = tuneParamOEC$lambda,
                                                  method = "glmnet",
                                                  nfolds = K,
                                                  nnlsInd = TRUE,
                                                  OECeta = OECeta,
                                                  sampSzWeight = sampSzWeight,
                                                  glmnetOpt = glmnetOpt,
                                                  xStandardize = FALSE, # I just chose this as FALSE 9/2/20
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
                                          nfolds = "K",
                                          nnlsInd = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          cv = etaTune,
                                          standardize = TRUE,
                                          weights = NULL,
                                          glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
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
                                              nfolds = "K",
                                              nnlsInd = TRUE,
                                              sampSzWeight = sampSzWeight,
                                              cv = etaTuneIndSpec, #"cvSpecTS",
                                              standardize = TRUE,
                                              weights = NULL,
                                              glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
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
                                               nfolds = "K",
                                               nnlsInd = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               cv = etaTuneIndSpec0, #"cvSpecTS0",
                                               standardize = TRUE,
                                               weights = NULL,
                                               glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
                                               xStandardize = FALSE,
                                               AvgW = AvgW,
                                               SpecID = test_country,
                                               horizon = horizon
                    )
                    
                }else{
                    OECetaSpec0 <- OECeta
                }
                
                
            }
            
            # 
            # if(stackCV == "zero"){
            #     
            #     stackParamOEC <- 0 # zero out
            #     # if(psiL + psiH > 0){
            #     #     # if either are nonzero 
            #     # }
            # }else if(stackCV == "cv" | stackCV == "hoso"){
            #     ############################################
            #     # Tune Stacking Hyperparameter again with tuned eta
            #     ############################################
            #     # consider switching order to be before tuneParam in reTuning
            #     stackParamOEC <- oecW_CV(data = full,
            #                              tune.grid = tune.grid,
            #                              sslLambdas = tuneParamOEC$lambda,
            #                              method = "glmnet",
            #                              nfolds = K,
            #                              nnlsInd = TRUE,
            #                              OECeta = OECeta,
            #                              sampSzWeight = sampSzWeight,
            #                              glmnetOpt = glmnetOpt,
            #                              xStandardize = FALSE, # I just chose this as FALSE 9/2/20
            #                              cv = stackCV,
            #                              standardize = TRUE,
            #                              weights = Diagonal( nrow(full) ),
            #                              horizon = horizon )
            #     
            #     stackParamOEC <- stackParamOEC$lambda
            # }
            # 
            # if(stackCVSpec == "cvSpecTS"| stackCVSpec == "hosoSpecTS" | stackCVSpec == "cvSpec" |
            #    stackCVSpec == "cvSpecTSHOO"){
            #     ############################################
            #     # Tune Stacking Hyperparameter again with tuned eta
            #     ############################################
            #     # consider switching order to be before tuneParam in reTuning
            #     stackParamOECSpec <- oecW_CV(data = full,
            #                                  tune.grid = tune.grid,
            #                                  sslLambdas = tuneParamOEC$lambda,
            #                                  method = "glmnet",
            #                                  nfolds = K,
            #                                  nnlsInd = TRUE,
            #                                  OECeta = OECeta,
            #                                  sampSzWeight = sampSzWeight,
            #                                  glmnetOpt = glmnetOpt,
            #                                  xStandardize = FALSE, # I just chose this as FALSE 9/2/20
            #                                  cv = stackCVSpec,
            #                                  standardize = TRUE,
            #                                  weights = Diagonal( nrow(full) ),
            #                                  SpecID = test_country,
            #                                  horizon = horizon )
            #     
            #     stackParamOECSpec <- stackParamOECSpec$lambda
            # }else if(stackCVSpec == "standard"){
            #     
            #     # use from above (standard stacking) -- tuned above in standard stacking section
            #     stackParamOECSpec <- stackParamSpec
            #     
            # }else{
            #     
            #     #otherwise use from above generalist OEC 
            #     stackParamOECSpec <- stackParamOEC
            # }
            # 
            # # zero out specialist
            # if(stackCVSpec0 == "cvSpecTS0"| stackCVSpec0 == "hosoSpecTS0" | stackCVSpec0 == "cvSpec0" |
            #    stackCVSpec == "cvSpec0TSHOO"){
            #     
            #     ############################################
            #     # Tune Stacking Hyperparameter again with tuned eta
            #     ############################################
            #     # consider switching order to be before tuneParam in reTuning
            #     stackParamOECSpec0 <- oecW_CV(data = full,
            #                                   tune.grid = tune.grid,
            #                                   sslLambdas = tuneParamOEC$lambda,
            #                                   method = "glmnet",
            #                                   nfolds = K,
            #                                   nnlsInd = TRUE,
            #                                   OECeta = OECeta,
            #                                   sampSzWeight = sampSzWeight,
            #                                   glmnetOpt = glmnetOpt,
            #                                   xStandardize = FALSE, # I just chose this as FALSE 9/2/20
            #                                   cv = stackCVSpec0,
            #                                   standardize = TRUE,
            #                                   weights = Diagonal( nrow(full) ),
            #                                   SpecID = test_country,
            #                                   horizon = horizon )
            #     
            #     stackParamOECSpec0 <- stackParamOECSpec0$lambda
            # }else if(stackCVSpec == "standard"){
            #     
            #     # use from above (standard stacking) -- tuned above in standard stacking section
            #     stackParamOECSpec0 <- stackParamSpec0 # zero out
            #     
            # }else{
            #     
            #     #otherwise use from above generalist OEC 
            #     stackParamOECSpec0 <- stackParamOEC
            # }
            
            
            if(tuneParamZero)   tuneParam$lambda <- rep(0, K)
            
            
            warmRes <- ridgeWS(data = full, 
                                tuneParam = tuneParam$lambda, 
                                stackParam = stackParamOEC, 
                                nnlsInd = TRUE,
                                stackTune = TRUE, # very important keep at TRUE 8/18/20
                                modelStandardize = TRUE,
                                stackStandardize = TRUE,
                                glmnetOpt = glmnetOpt, # used closed form
                                xStandardize = FALSE,
                                sampSzWeight = 1, 
                                weights = Diagonal( nrow(full) ),
                                lambdaGrid = lambdaVec,
                                AvgW = FALSE) # positive results achieved at FALSE 8/18/20 -- do not change
            
            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack
            
            # if(etaTuneInd){
            #     
            #     OECeta <-   oecEta_CV(data = full,
            #                           tune.grid = etas,
            #                           sslLambdas = tuneParam$lambda,
            #                           stackParam = stackParamOEC,
            #                           method = "glmnet",
            #                           nfolds = "K",
            #                           nnlsInd = TRUE,
            #                           sampSzWeight = sampSzWeight,
            #                           cv = etaTune,
            #                           standardize = TRUE,
            #                           weights = NULL,
            #                           glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
            #                           xStandardize = FALSE,
            #                           AvgW = AvgWEtaTn,
            #                           horizon = horizon)
            #     
            # }
            # 
            # if(etaTuneIndSpec != FALSE){
            #     
            #     OECetaSpec <-   oecEta_CV(data = full,
            #                           tune.grid = etas,
            #                           sslLambdas = tuneParam$lambda,
            #                           stackParam = stackParamOECSpec,
            #                           method = "glmnet",
            #                           nfolds = "K",
            #                           nnlsInd = TRUE,
            #                           sampSzWeight = sampSzWeight,
            #                           cv = etaTuneIndSpec, #"cvSpecTS",
            #                           standardize = TRUE,
            #                           weights = NULL,
            #                           glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
            #                           xStandardize = FALSE,
            #                           AvgW = AvgWEtaTn,
            #                           SpecID = test_country,
            #                           horizon = horizon)
            #     
            # }else{
            #     OECetaSpec <- OECeta
            # }
            # 
            # 
            # if(etaTuneIndSpec0 != FALSE){
            #     
            #     OECetaSpec0 <-   oecEta_CV(data = full,
            #                                tune.grid = etas,
            #                                sslLambdas = tuneParamOEC$lambda,
            #                                stackParam = stackParamOECSpec0,
            #                                method = "glmnet",
            #                                nfolds = "K",
            #                                nnlsInd = TRUE,
            #                                sampSzWeight = sampSzWeight,
            #                                cv = etaTuneIndSpec0, #"cvSpecTS0",
            #                                standardize = TRUE,
            #                                weights = NULL,
            #                                glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
            #                                xStandardize = FALSE,
            #                                AvgW = AvgW,
            #                                SpecID = test_country,
            #                                horizon = horizon
            #                                )
            #     
            # }else{
            #     OECetaSpec0 <- OECeta
            # }
            
            
                resMat[iterNum, 22] <- OECeta
                
                if(lowLimTune){
                    psiL <- oecLow_CV(data = full,
                                      tune.grid = psiVec,
                                      sslLambdas = tuneParam$lambda,
                                      stackParam = stackParamOEC,
                                      oecEta = OECeta,
                                      method = "glmnet",
                                      nfolds = "K",
                                      nnlsInd = TRUE,
                                      sampSzWeight = sampSzWeight,
                                      cv = etaTune,
                                      standardize = TRUE,
                                      weights = NULL,
                                      glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
                                      xStandardize = FALSE,
                                      AvgW = AvgWEtaTn)
                    
                }
                
                ### OEC
                if(sampSzWeight < 6){
                    # giving me problems for 6 and 7
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
                
                    
                    # if(intercept){
                    #     predsVecOEC <- oecRidgePred(data = X_test[,-c(1,2)],
                    #                                 mod = oecRidgeWS)
                    # }else{
                    # predsVecOEC <- oecRidgePred(data = X_test[,-1],
                    #                            mod = oecRidgeWS)
                    #}
                    # resMat[iterNum, 15] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                    
                
                
                ######################
                # Specialist OEC -- Country specific OEC -- initialize on country-specific stacking
                ######################
                # rerun Warm Start with different \mu (stackParam) for specialist
                warmRes <- ridgeWS(data = full, 
                                   tuneParam = tuneParam$lambda, 
                                   stackParam = stackParamOECSpec, 
                                   nnlsInd = TRUE,
                                   stackTune = TRUE, # very important keep at TRUE 8/18/20
                                   modelStandardize = TRUE,
                                   stackStandardize = TRUE,
                                   glmnetOpt = glmnetOpt, # used closed form
                                   xStandardize = FALSE,
                                   sampSzWeight = 1, 
                                   weights = Diagonal( nrow(full) ),
                                   lambdaGrid = lambdaVec,
                                   AvgW = FALSE) # positive results achieved at FALSE 8/18/20 -- do not change
                
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
                                   tune.grid = seq(1 / K ,1, length = 20), # psiVec,
                                   sslLambdas = tuneParamOEC$lambda,
                                   stackParam = stackParamOEC,
                                   oecEta = OECeta,
                                   method = "glmnet",
                                   nfolds = "K",
                                   nnlsInd = TRUE,
                                   sampSzWeight = sampSzWeight,
                                   cv = "cv",
                                   SpecID = test_country,
                                   standardize = TRUE,
                                   weights = NULL,
                                   glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
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
                                               xStandardize = TRUE, #************
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
                                   nfolds = "K",
                                   nnlsInd = TRUE,
                                   sampSzWeight = sampSzWeight,
                                   cv = "cvWindow",
                                   SpecID = test_country,
                                   standardize = TRUE,
                                   weights = NULL,
                                   glmnetOpt = glmnetOpt, # was TRUE until 8/16/20
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
                                                     mu = 0,#stackParamOECSpec0, 
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
                                                     # xStandardize = TRUE, #************
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
                resMat[iterNum, 42] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
                resMat[iterNum, 44] <- psiH
                
                # -----------------------------------
                #####
                # "eta1", "etaSpec", "etaSpec0", "mu", "muSpec", "muSpec0"
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
        
    ## next part, only for test country
    ##### next train up to 2 years prior to test, then stack on 1 year before
    
    
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
                       "oec_SpecUp", "oecWindow", "up", "upWindow")



# rownames(resMat) <- countries

print("save file")
setwd(save.folder)
write.csv(resMat, fileNm) #, row.names = FALSE)

# colMeans(resMat, na.rm = TRUE) # write files

# boxplot(resMat[,3] / resMat[,1], resMat[,10] / resMat[,1],
#         main = paste0("SampSzWt: ",sampSzWeight),
#         names =c("stacking", "oec"),
#         ylab = "RMSE / RMSE_Merging")
# abline(h = 1, col = "blue")
# 
# mean(resMat[,10] / resMat[,3])


# 2.7 RMSE: avg weights, sampSzWt =1 , OECeta = 0.9
# 2.73 but 75% quantile below Merging, SampSzWt = 6