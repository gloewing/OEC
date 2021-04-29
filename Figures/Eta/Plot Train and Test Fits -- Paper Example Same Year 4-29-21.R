cnt = 2 # netherlands
iterNum = 8 # test year 2015

# norway (cnt = 17) 2015 (iterNum = 4), -- not to bad but zeroOut looks better
# iceland (19) 2019
# netherlands (2) 2015
# original was netherlands 2016 zero out
source("msMort_functions.R")
source("Study Strap Functions.R")
source("OEC Functions.R")

pM <- expand.grid( c(FALSE), #nH
                   c(100), # minTrain
                   c(4), # xPro
                   c(1)) # sampSzWt

itr <- 1

# test set date ranges
month_start <- "-01-01" #
month_end <-"-12-31" #
addYear <- FALSE # addyear means add a year to the test date--if start and end are the same year then set to FALSE
testYears <- 2011:2019 # 2010:2019 #chosen because many countries in March 2012 have just over 2 years worth of data, and
# also allows us to go up until March 2020, right as pandemic was beginning
# 2019 goes up to 2020 because the test year starts in March

minTrain <- pM[itr, 2] # minimum number of training data points to be considered a training country for any iterations
minTest <- 50 # minimum number of testing data points to be considered a test country for any iterations
testDate_start <- "-01-01" #"-03-01" 
testDate_end <- "-12-31" #"-02-28"
testCountryData <- 12 # this is the number of months that we pretend the test country has (different than minTest)
nH <- pM[itr, 1] # northern hemisphere
# to emulate south africa

dataNm <- "world-mort_old.csv"
full <- read.csv( dataNm )
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
oecTune <- "zero"#"cvCFTS" #"cvOEC" # "cvCFTS"#"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
# added in
stdTn <-  "zero"#"cvCFTS" #"zero"#"cvCFTS" #"sseCF"# "cvCFTS"# "cvCF"  # "sseCF" # "cvCFTS" # zero" #"
############
stackCV <- "zero" #"cv" # can be "1", "cv", "hoso" or "zero." If 1, then just use the same parameter as in standard vanilla stacking tuned above. If "zero", then set stacking hyperparamter to 0 if using window stacking
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
lambdaVec <- tune.grid <- c(0, 1e-14)#sort( unique( c(0.0001, 0.001, 0.01,
                                          # exp(-seq(0,5, length = 50)),
                                          # seq(5,20, by = 5),
                                          # seq(30,100, by = 10) ) ) ) # 2:100

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


full <- read.csv( dataNm )
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
full <- read.csv( dataNm )
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
        
        
    }
    
    # concatenate countries
    X <- do.call(rbind, cList)
    rm(cList)
    
    # update these since they are reordered
    full <- do.call(rbind, infoList)
    colnames(full)[nm] <- "Y"
    rm(infoList)
    
    colnames(test)[nm] <- "outcome"
    test <- test[order(test$date),]
    X_test <- Xgen(test)
    colnames(test)[nm] <- "Y"

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
preds0 <- w[1] + X_test %*% betas %*% w[-1]
resMat[iterNum, 8] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
resMat[iterNum, 40] <- stackParam
############################################################################

# 
# # country specific model
# plot(c(full$Y[trainIndx], test$Y),
#      ylab = "Population Adjusted Deaths",
#      xlab = "Months",
#      pch = 19,
#      lwd = 0.1)
# 
# abline(v = length(trainIndx), lwd = 2, lty = 2)
# lines(rbind(X[trainIndx,], X_test) %*% bCountry, lwd = 1.5)
# 
# # zero out Country specific stacking
# w <- wZeroOut
# preds <- w[1] + rbind(X[trainIndx,], X_test) %*% betas %*% w[-1]
# abline(v = length(trainIndx))
# lines(preds, col = "darkgreen", lwd = 1.5)

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

resMat[iterNum, 22] <- OECeta

    # giving me problems for 6 and 7
   

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

######################
# Specialist OEC -- Country2 specific OEC-- initialize on country-specific stacking
######################
etasVec <- c(1e-14, 1e-13,  1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 
             0.9, 0.99, 0.999, 0.9999, 0.99999, OECetaSpec)
etaMat <- matrix(nr = length(etasVec), nc = nrow(rbind(X[trainIndx,-1], X_test[,-1])) + 2 )

indx1 <- which(full$Study == test_country)

for(j in 1:length(etasVec)){
    e <- etasVec[j]
    etaMat[j,1] <- e
    
    pMat <- matrix(ncol = itrs, nrow = nrow(rbind(X[trainIndx,-1], X_test[,-1])) ) # store predictions in each column 
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
                                   eta = e,
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
                                   up = Inf,
                                   stackTol = 1e-9
        )
        
        pMat[,it] <- oecRidgePred(data = rbind(X[trainIndx,-1], X_test[,-1]),
                                  mod = oecRidgeWS)
        
    }
    
    # save predictions
    etaMat[j, -c(1,2)] <- predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
    # save RMSE at value of eta
    
    #}
}

################
# zero out
################
# etasVec <- c( 0.1,#1e-14, 1e-6, 
#              0.9, 0.99, 0.995, 0.9975, 0.999, 0.9995, 0.9999, 0.99999, OECetaSpec0)
etasVec <- c(1e-14, 0.1, 0.25, 0.5, 0.75, 0.8, 0.825, 0.85, 0.875,
             0.9, 0.9275, 0.93, 0.9325, 0.933, 0.934, 0.935, 0.936, 0.9375, 0.94, 0.95, 0.99, 0.999, 0.9999, OECetaSpec0)

etaMat0 <- as.data.frame( matrix(nr = length(etasVec) + 3, nc = nrow(rbind(X[trainIndx,-1], X_test[,-1])) + 2 ) )
etaMat0[,1] <- c(etasVec, NA, NA, NA) # first column are eta values
etaMat0[,2] <- "preds" # second column says whether predicted or observed. last row is observed but thats changed below

for(j in 1:length(etasVec)){
    e <- etasVec[j]
    # etaMat0[j,1] <- e
    
    pMat <- matrix(ncol = itrs, nrow = nrow(rbind(X[trainIndx,-1], X_test[,-1])) ) # store predictions in each column 
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
                                    eta = e,
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
        
        pMat[,it] <- oecRidgePred(data = rbind(X[trainIndx,-1], X_test[,-1]),
                                  mod = oecRidgeWS)
        
    }
    
    # save predictions
    etaMat0[j, -c(1,2)] <- predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
    # save RMSE at value of eta
    
    #}
}
# last row are actual observed values
etaMat0[nrow(etaMat0), -c(1,2)] <- c(full$Y[trainIndx], test$Y) # actual observed values
etaMat0[nrow(etaMat0), 2] <- "obs"


# second to last row are the study specific model
etaMat0[nrow(etaMat0) - 1, -c(1,2)] <- rbind(X[trainIndx,], X_test) %*% bCountry
etaMat0[nrow(etaMat0) - 1, c(1,2)] <- c(NA, "country_specific_model")

# 3rd to last row is zero out stacking
w <- wZeroOut # wZeroOut
preds <- w[1] + rbind(X[trainIndx,], X_test) %*% betas %*% w[-1]
etaMat0[nrow(etaMat0) - 2, -c(1,2)] <- preds
etaMat0[nrow(etaMat0) - 2, c(1,2)] <- c(NA, "stacking")

# 3rd

colnames(etaMat0) <- c("eta", 
                       "type", 
                       paste0("week_", 1:(ncol(etaMat0) - 2) ) )
write.csv(etaMat0, "~/Desktop/Research Final/Mortality/Figures/Final Figures/Eta Plots/etaData")
par(mfrow = c(1,1))
# etaMat[j, 2] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
############################################################################

############################################################################

# country specific model
plot(c(full$Y[trainIndx], test$Y),
     x = 1:length( c(full$Y[trainIndx], test$Y) ),
     ylab = "Population Adjusted Deaths",
     xlab = "Months",
     pch = 19,
     lwd = 0.1,
     main = "OEC vs. Stacking (Netherlands 2015)")

abline(v = length(trainIndx), lwd = 2, lty = 2, col = "darkgray")
lines( rbind(X[trainIndx,], X_test) %*% bCountry, lwd = 1.5, col = "green")

# use last years fitted values for predictions as reference
# lines( y = rbind(X[trainIndx,]) %*% bCountry, 
#        x = 51:100,
#        lwd = 1.5, col = "darkgreen")


# zero out Country specific stacking
w <- wZeroOut # wZeroOut
preds <- w[1] + rbind(X[trainIndx,], X_test) %*% betas %*% w[-1]
lines(y = preds, 
      x = 1:length( c(full$Y[trainIndx], test$Y) ),
      col = "red", lwd = 1.5)

# points from last year as a reference plotted as if they happened in test year
# points(y = full$Y[trainIndx],
#        x = 51:100)

# rmses
sqrt(mean( (test$Y - preds)^2 )) # zero out
sqrt(mean( (test$Y - full$Y[trainIndx])^2 )) # last years points as "predictions"
sqrt(mean( (test$Y - X[trainIndx,] %*% bCountry)^2 )) # last years fitted values

# zero out OEC
for(j in 1:(nrow(etaMat0) -2)){
    preds <- etaMat0[j, -c(1,2)] # last one is model with tunes eta
    lines(y = preds, 
          x = 1:length( c(full$Y[trainIndx], test$Y) ),
          col = "blue", lwd = 1.5)
}

# do again so can be seen
abline(v = length(trainIndx), lwd = 2, lty = 2, col = "darkgray")
lines( rbind(X[trainIndx,], X_test) %*% bCountry, lwd = 1.5, col = "green")

# zero out Country specific stacking
w <- wZeroOut # wZeroOut
preds <- w[1] + rbind(X[trainIndx,], X_test) %*% betas %*% w[-1]
lines(y = preds, 
      x = 1:length( c(full$Y[trainIndx], test$Y) ),
      col = "red", lwd = 1.5)
