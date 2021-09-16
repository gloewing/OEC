library(doParallel)
library(caret)
library(glmnet)
library(foreach)

source("Study Strap Functions.R")
source("SimFn_multiCountry.R")
source("OEC Functions.R")

# 
# # 
# sims6 <- cbind(325:378, expand.grid(c(0.25, 1, 4),
#                                     c(3, 6, 39), # 3 ,7, 11, 26
#                                     c("cvCF", "zero"), # "cvCF", "zero", "sseCF"  "cvCF"
#                                     c("100"),
#                                     c("2003", "2010", "2015")
#                                     )
#                )
# colnames(sims6) <- c("simNum", "betaVar", "K", "tune", "minTrain", "yr")
# #
# write.csv(sims6, "paramMat", row.names = FALSE)

sims6 <- read.csv("paramMat")

save.folder <- "/n/home12/gloewinger/simsMort4"
load.folder <- "~/Desktop/Research"

# sim params
simNum <- 511
totalSims <- 100

runNum <- 1 

sampSzWeight <- 1

# sim parameters
totalStudies <- sims6[runNum, 3]  # need number that is divisible by 2, 4, 8 for cluster purposes and structure of 
clusts <- totalStudies #sims6[runNum, 4] # totalStudies
betaVar <- sims6[runNum, 2] 
covVar <- 0 #sqrt( sims6[runNum, 3] )
trainStudies <- 1:(totalStudies - 1)
K <- length(trainStudies)
num.trainStudy <- length(trainStudies)
sample.size <- NA # if NA then randomly determine (uniformly) sample sizes in the function-- keep like this
testSetLength <- 52 # number of actual test set dates
testTrain <- 52 # length of test country's training set
minTime <- sims6[runNum, 5] # minimum amount of training dates for a training country 
nLow <- 50 # samp size lower
nHigh <- 150 # samp size upper
epsLow <- 1 # noise lower
epsHigh <- 2 # noise upper

# OEC Params
OECeta <- 0.5
pcaInd <- FALSE
TScv <- FALSE
etaTune <- "cv"
etaTuneInd <- TRUE
tuneParamZero <- FALSE # set all study specific alphas to 0
stackTol <- 1e-9 # tolerance for stacking regression for PGD
oecTune <- as.character( sims6[runNum,4] ) 
stdTn <- as.character( sims6[runNum,4] ) # this is hyperparameters for standard stacking
stackCV <- "cv" 
stackCVSpec <- FALSE # 
stackCVSpec0 <- FALSE # 
etaTuneIndSpec <- TRUE # tune eta for specialists (and 0 out specialists)
stackTune <- "cv" 
stackTuneInd <- TRUE
scaleInd <- TRUE # scale covaraites
glmnetOpt <- FALSE # 1
itrs <- 20 # number of runs of OEC
lowLimTune <- FALSE
TScv <- FALSE 
psiVec <- seq(0,1, length = 11)
pcaInd <- FALSE
AvgW <- FALSE
psiL <- 0
psiH <- Inf
reTune <- FALSE
etaSet <- TRUE # include extra etas
orderRandom <- FALSE # randomize order of blocks
avgOpt <- FALSE # use average or optimal value of multiple runs of OEC
yr <- 2003 # start of training year

nfolds <- "K"

test_study <- max(trainStudies) + 1 # arbitrarily set to study but any of the non-training studies is fine (11-24 are all random test studies)

lamb <- 0.5
fileNm <- paste0("oecMort_",
                 OECeta, ".etTn_",
                 etaTuneInd, "_", etaTune, 
                "_stTn_", stackTune, 
                "_oecTn_", oecTune,
                 "_stTnIn_", stackTuneInd, 
                "_stCV_", stackCV,
                "_etaSpc_", etaTuneIndSpec,
                "_Wcv_", stackCV, 
                "_Wspc_", stackCVSpec,
                "_sclX_", scaleInd,
                "_glm_", glmnetOpt,
                "_sampSz_", sample.size,
                  "_smpSzW_", sampSzWeight, 
                "_tstTrn_", testTrain, # length of test country's training set
               "_minTrn_", minTime, # min training
                "_bVar_", betaVar,
                 "_reTn_", reTune, 
               "_eta_", etaSet, # set of eta
               "_yr_", yr,
               "_K_", K) 
print(fileNm)

if(!etaSet){
    etas <- sort( c(0.001, 0.01, seq(0.05, 0.95, 0.05), 0.99, 0.999) ) # hyperparameter for stacking and SSL losses
    
}else{
    etas <- sort( c(0.001, 0.01, seq(0.05, 0.95, 0.05)) ) # , 0.99, 0.999 # hyperparameter for stacking and SSL losses
    etas <- sort( c(etas, exp(-seq(5,40, length = 40))
              ) )
}

tune.grid <- sort( unique( c(0.0001, 0.001, 0.01, 
                             exp(-seq(0,5, length = 50)),
                             seq(2,4),
                             seq(5,20, by = 5),
                             seq(30,100, by = 10) ) ) ) # 2:100

tune.grid <- as.data.frame(tune.grid) # tuning parameters to consider
tune.grid <- cbind(0, tune.grid) # Ridge
colnames(tune.grid) <- c("alpha", "lambda")
warmStartSolver <- "glmnet" # NULL will set to my own function
# parallelize

######### Parallelize
logfile <- paste0("outputFile", simNum,".txt")
writeLines(c(""), file(logfile,'w'))

num.threads <- as.integer( ceiling( totalSims / 3 ) )# round up
threads <- makeCluster(num.threads, outfile=logfile)
registerDoParallel(threads)

setwd(save.folder)

getDoParWorkers()
timeStart <- Sys.time()

results <- foreach(iterNum = 1:totalSims, .combine = list, .multicombine = TRUE) %dopar%{
    
    print(paste0("start: ", iterNum, "_lambda: ", lamb))
    
    library(CVXR)
    library(caret)
    library(glmnet)
    library(dplyr)
    resMat <- matrix(nrow = totalSims, ncol = 65)

    colnames(resMat) <- c("merge", "avg", "stacking", "stacking_zeroOut", 
                          "merge2", "stack2", 
                          "stacking_country", "stacking_country_zeroOut",
                          "country", "oec", "oec2", "oec_test2", "oecSpec_test2",
                          "avg", "oecNoInt", "oec_country", "oec_country2", "oec_country3",
                          "oec_country4", "oec_country5",
                          "oec_countryAvg",
                          "eta", "nnet", "nnetAvg", "K", "n_k", "n_test", "N",
                          "oec_countryAvgStack", "oec_AvgStack", "oec_country0", "oec_country6",
                          "oec_countryWindow",
                          "oec_SpecPCA", "objSpec", "objZero", "objGen",
                          "oec_specAnneal_left", "oec_specAnneal_leftObj", "oec_spec_left", "oec_spec_leftObj",
                          "oec_specAnneal_right", "oec_specAnneal_rightObj", "oec_spec_right", "oec_spec_rightObj",
                          "oec_0Anneal_left", "oec_0Anneal_leftObj", "oec_0_left", "oec_0_leftObj",
                          "oec_0Anneal_right", "oec_0Anneal_rightObj", "oec_0_right", "oec_0_rightObj",
                          "oec_genAnneal_left", "oec_genAnneal_leftObj", "oec_gen_left", "oec_gen_leftObj",
                          "oec_genAnneal_right", "oec_genAnneal_rightObj", "oec_gen_right", "oec_gen_rightObj",
                          "etaGen", "etaSpec", "etaSpec0", "linearOnly")
    
        
        # save results
        studyNum <- test_study
        cnt <- seedSet <- iterNum + 400 # ensures no repeats
        set.seed(seedSet)
        
        # simulate data
        full <- multiStudySim_mort(sim.num = simNum,
                                  iter = iterNum,
                                  sampSize = sample.size, #if NA then randomly draw sample sizes according to dates given in mort data
                                  beta.var = betaVar, # variance of betas
                                  min.time = minTime, # minimum number of weeks for training dataset
                                  minDate = paste0(yr, "-01-01"),
                                  num.studies = totalStudies,
                                  testStudy = totalStudies - 1, # arbitrarily set test study to the last one for sample size purposes
                                  testSampSz = testTrain, # number of observations in training set for test country
                                  testSet.length = testSetLength, # length of test set
                                  exchg.testStudy = totalStudies, # study index of exchangable test study
                                  exchg.testSet.length = testSetLength, # number of observations in test set of exchangable study
                                  testDateStart = "2019-01-01",
                                  studyNoise = NA #c(1, 1), # if NA, use the noise from the mort data, otherwise sample uniformly with specified amout c(epsLow, epsHigh), # range of noises for the different studies    
                                  )
        
        # use just the datasets for test set and full merged datasets
        test <- as.data.frame( full$test ) # specialist test set
        test2 <- as.data.frame( full$exchngTest ) # exchangabletest set
        full <- as.data.frame( full$data ) # merged training dataset
        
        ## specialist
        test_country <- totalStudies - 1 # indx of test country is just the last study (arbitrary)
    
        # # rename the specialist rows as the last training study
        trainIndx <- specIndx <- which(full$Study == test_country)
    
        countries <- unique(full$Study) # only include countries that have both
        K <- length(countries) # number of training countries

        resMat[iterNum, 25] <- K # number of training countries including the test country
        resMat[iterNum, 26] <- testTrain # number of observations in training set for test country
        resMat[iterNum, 27] <- nrow(test) # number of observations in training set for test country
        resMat[iterNum, 28] <- nrow(full) # size of training set
        
        #################
        # Study Labels
        #################
        test_country <- which(countries == test_country) # number corresponding to test country
        full$Study <- as.numeric(full$Study) # turn numeric
        studyVec <- countries <- unique(full$Study)
        num.trainStudy <- length(studyVec)
        indxVec <- vector(length = (nrow(full)))
        
        for(j in 1:K){
            
            indx <- which(full$Study == studyVec[j])
            indxVec[indx] <- j
            
        }
        
        full$Study <- indxVec # replace with new study labels from 1:K
        countries <- unique(full$Study) # number corresponding to it
        
        X <- full[,-c(1,2)] 
        X_test <-  test[,-c(1,2)]

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
                test2[,-c(1,2)] <- (test2[, column + 2] - means[column]) / sds[column]
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
                test2[, column + 2] <- (test2[, column + 2] - means[column]) / sds[column]
            }
        }
        
        ####################################
        # Make Design Matrix with Full
        ####################################

        X <- as.matrix( cbind(1, full[,-c(1,2)]) )
        X_test <- as.matrix( cbind(1, test[,-c(1,2)]) )
        test <- test[,c(1,2)] # remove design matrix because it is just extra memory
        
        if(pcaInd){
            pcs <- prcomp(X)
            # put last PC first so intercept (column of ones) is still first
            pcs$rotation <- cbind(pcs$rotation[, ncol(pcs$rotation)], pcs$rotation[, -ncol(pcs$rotation)])
            X <- cbind(1, pcs$x[,-ncol(pcs$x)] ) # remove last column because just 0s
            X_test <- X_test %*% pcs$rotation
            rm(pcs)
        }
        
        set.seed(iterNum) # added 9 /3/20
        
        #######################
        # Tuning
        #######################
        # Merge Tune
        
        set.seed(cnt) 
        if(stdTn == "zero" & oecTune == "zero"){
            mergedLambda <- data.frame(alpha = 0, lambda = 0)
        }else{
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
        }
        
        # standard stacking tune
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
                      lambda = mergedLambda$lambda,
                      standardize = TRUE,
                      intercept = TRUE,
                      thresh = 1e-10 )
        
        preds <- X_test %*% coef(fit) # predict(fit, X_test)
        mrgBeta <- coef(fit)
        resMat[iterNum, 1] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
        
        # test 2 merging
        preds <- as.matrix( cbind(1, test2[,-c(1,2)] ) ) %*% coef(fit) # predict(fit, X_test)
        mrgBeta <- coef(fit)
        resMat[iterNum, 5] <- sqrt( mean( (test2$Y - preds)^2 ) ) # RMSE
        
        ############
        # stacking 
        ############
        
        predsMat <- matrix(nrow = nrow(full), ncol = K )
        predsMat_testCountry <- matrix(nrow = length(trainIndx), ncol = K )
        betas <- matrix(nrow = ncol(X), ncol = K)
        
        set.seed(cnt) # added 9 /3/20
        
        for(k in 1:K){
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
                rm(fit)
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
                rm(fit)
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
                
                bCountry_noLinear <- b
                preds <- X_test[,-2] %*% b
                resMat[iterNum, 65] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE

                
            }
            
            
        }
        
        predsMat0 <- predsMat # save copy of predsMat
        
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
                             lambda = tune.grid$lambda, 
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
        
        # standard stacking predictions test 2
        preds <- w[1] + as.matrix( cbind(1, test2[,-c(1,2)] ) ) %*% betas %*% w[-1]
        resMat[iterNum, 6] <- sqrt( mean( (test2$Y - preds)^2 ) ) # RMSE
        
        
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
                             lambda = tune.grid$lambda, 
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
        # standard stacking on country
        ######################
        # indices correpsonding to test country in training set for year prior
        indx1 <- which(full$Study == test_country)
        
        modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                             x = as.matrix(predsMat_testCountry),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
        
        stackParam <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)

        if(stackCV == "zero")   stackParam <- 0
        
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
        stackParamSpec <- stackParam
        
        rm(mod)
        
        # standard stacking predictions
        preds <- w[1] + X_test %*% betas %*% w[-1]
        resMat[iterNum, 7] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
        
        #######################
        # standard stacking on country zero out
        ######################
        # indices correpsonding to test country in training set for year prior
        indx1 <- which(full$Study == test_country)
        predsMat_testCountry2 <- predsMat_testCountry
        predsMat_testCountry2[,test_country] <- 0 # zero out country thats being tested
        
        modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                             x = as.matrix(predsMat_testCountry2),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
            
        stackParam <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)    
            
        if(stackCV == "zero")   stackParam <- 0
        
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
        stackParamSpec0 <- stackParam
        
        rm(mod)
        
        # standard stacking predictions
        preds <- w[1] + X_test %*% betas %*% w[-1]
        resMat[iterNum, 8] <- sqrt( mean( (test$Y - preds)^2 ) ) # RMSE
        ############################################################################
        
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
                                     nfolds = K,
                                     nnlsInd = TRUE,
                                     OECeta = OECeta,
                                     sampSzWeight = sampSzWeight,
                                     glmnetOpt = glmnetOpt,
                                     xStandardize = FALSE, 
                                     cv = stackCV,
                                     standardize = TRUE,
                                     weights = Diagonal( nrow(full) ) )
            
            stackParamOEC <- stackParamOEC$lambda
        }
        
        if(stackCVSpec == "cvSpecTS"| stackCVSpec == "hosoSpecTS" | stackCVSpec == "cvSpec"){
            ############################################
            # Tune Stacking Hyperparameter again with tuned eta
            ############################################
            stackParamOECSpec <- oecW_CV(data = full,
                                         tune.grid = tune.grid,
                                         sslLambdas = tuneParamOEC$lambda,
                                         method = "glmnet",
                                         nfolds = K,
                                         nnlsInd = TRUE,
                                         OECeta = OECeta,
                                         sampSzWeight = sampSzWeight,
                                         glmnetOpt = glmnetOpt,
                                         xStandardize = FALSE,
                                         cv = stackCVSpec,
                                         standardize = TRUE,
                                         weights = Diagonal( nrow(full) ),
                                         SpecID = test_country )
            
            stackParamOECSpec <- stackParamOECSpec$lambda
        }else if(stackCVSpec == "standard"){
            
            # use from above (standard stacking) -- tuned above in standard stacking section
            stackParamOECSpec <- stackParamSpec
            
        }else{
            
            #otherwise use from above generalist OEC 
            stackParamOECSpec <- stackParamOEC
        }
        
        # zero out specialist
        if(stackCVSpec0 == "cvSpecTS0"| stackCVSpec0 == "hosoSpecTS0" | stackCVSpec0 == "cvSpec0"){
            ############################################
            # Tune Stacking Hyperparameter again with tuned eta
            ############################################
            stackParamOECSpec0 <- oecW_CV(data = full,
                                         tune.grid = tune.grid,
                                         sslLambdas = tuneParamOEC$lambda,
                                         method = "glmnet",
                                         nfolds = K,
                                         nnlsInd = TRUE,
                                         OECeta = OECeta,
                                         sampSzWeight = sampSzWeight,
                                         glmnetOpt = glmnetOpt,
                                         xStandardize = FALSE, 
                                         cv = stackCVSpec0,
                                         standardize = TRUE,
                                         weights = Diagonal( nrow(full) ),
                                         SpecID = test_country )
            
            stackParamOECSpec0 <- stackParamOECSpec0$lambda
        }else if(stackCVSpec == "standard"){
            
            # use from above (standard stacking) -- tuned above in standard stacking section
            stackParamOECSpec0 <- stackParamSpec0 # zero out
            
        }else{
            
            #otherwise use from above generalist OEC 
            stackParamOECSpec0 <- stackParamOEC
        }
        
        
        if(tuneParamZero)   tuneParam$lambda <- rep(0, K)
        
        
        warmRes <- ridgeWS(data = full, 
                           tuneParam = tuneParamOEC$lambda, 
                           stackParam = stackParamOEC, 
                           nnlsInd = TRUE,
                           stackTune = TRUE, 
                           modelStandardize = TRUE,
                           stackStandardize = TRUE,
                           glmnetOpt = glmnetOpt, # used closed form
                           xStandardize = FALSE,
                           sampSzWeight = 1, 
                           weights = Diagonal( nrow(full) ),
                           lambdaGrid = tune.grid$lambda,
                           AvgW = FALSE) 
        
        sigK <- warmRes$sigK
        sigStack <- warmRes$sigStack
        
        if(etaTuneInd){
            
            OECeta <-   oecEta_CV(data = full,
                                  tune.grid = etas,
                                  sslLambdas = tuneParamOEC$lambda,
                                  stackParam = stackParamOEC,
                                  method = "glmnet",
                                  nfolds = "K",
                                  nnlsInd = TRUE,
                                  sampSzWeight = sampSzWeight,
                                  cv = etaTune,
                                  standardize = TRUE,
                                  weights = NULL,
                                  glmnetOpt = glmnetOpt, 
                                  xStandardize = FALSE,
                                  AvgW = AvgW)
            
        }
        
        if(etaTuneIndSpec){
            
            OECetaSpec <-   oecEta_CV(data = full,
                                      tune.grid = etas,
                                      sslLambdas = tuneParamOEC$lambda,
                                      stackParam = stackParamOECSpec,
                                      method = "glmnet",
                                      nfolds = "K",
                                      nnlsInd = TRUE,
                                      sampSzWeight = sampSzWeight,
                                      cv = "cvSpec",
                                      standardize = TRUE,
                                      weights = NULL,
                                      glmnetOpt = glmnetOpt, 
                                      xStandardize = FALSE,
                                      AvgW = AvgW,
                                      SpecID = test_country)
            
        }else{
            OECetaSpec <- OECeta
        }
        
        if(etaTuneIndSpec){
            
            OECetaSpec0 <-   oecEta_CV(data = full,
                                      tune.grid = etas,
                                      sslLambdas = tuneParamOEC$lambda,
                                      stackParam = stackParamOECSpec0,
                                      method = "glmnet",
                                      nfolds = "K",
                                      nnlsInd = TRUE,
                                      sampSzWeight = sampSzWeight,
                                      cv = "cvSpec0",
                                      standardize = TRUE,
                                      weights = NULL,
                                      glmnetOpt = glmnetOpt,
                                      xStandardize = FALSE,
                                      AvgW = AvgW,
                                      SpecID = test_country)
            
        }else{
            OECetaSpec0 <- OECeta
        }
        
        
        resMat[iterNum, 22] <- OECeta
        
        if(lowLimTune){
            psiL <- oecLow_CV(data = full,
                              tune.grid = psiVec,
                              sslLambdas = tuneParamOEC$lambda,
                              stackParam = stackParamOEC,
                              oecEta = OECeta,
                              method = "glmnet",
                              nfolds = "K",
                              nnlsInd = TRUE,
                              sampSzWeight = sampSzWeight,
                              cv = etaTune,
                              standardize = TRUE,
                              weights = NULL,
                              glmnetOpt = glmnetOpt, 
                              xStandardize = FALSE,
                              AvgW = AvgW)
            
        }
        
        pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
        optObj <- vector(length = itrs) # saves optimal objective
        
        ### OEC
        if(sampSzWeight < 6){

            pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
            for(it in 1:itrs){
                set.seed(it)
                
                oecRidgeWS <- ridgeAltFix(data = full, 
                                          betaStart = warmRes$beta, 
                                          wStart = warmRes$w, 
                                          lambdaVec = tuneParamOEC$lambda, 
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
                                          stackTol = stackTol,
                                          orderRandom = orderRandom)
                
                pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                          mod = oecRidgeWS)
                
                optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
                
            }
            
            
            if(avgOpt){
                predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
            }else{
                bestObjIndx <- which.min(optObj)
                predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
            }
            
            
            resMat[iterNum, 10] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        }
        
        ############################ # update w 's only after all K beta_k have been updated
        pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
        pMat2 <- matrix(ncol = itrs, nrow = nrow(test2) ) # store predictions in each column 
        for(it in 1:itrs){
            set.seed(it)
            
            oecRidgeWS <- ridgeAlt(data = full, 
                                   betaStart = warmRes$beta, 
                                   wStart = warmRes$w, 
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            pMat2[,it] <- oecRidgePred(data = test2[,-c(1,2)],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC1 <- rowMeans( pMat ) # average predictions from the different runs
            predsVecOEC2 <- rowMeans( pMat2 ) # average predictions from the different runs
            resMat[iterNum, 37] <- mean(optObj) # minimum objective value achieved
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC1 <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
            predsVecOEC2 <- pMat2[,bestObjIndx] # choose predictions from the one with the best objective
            resMat[iterNum, 37] <- min(optObj) # minimum objective value achieved
            }
        
        resMat[iterNum, 11] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        
        resMat[iterNum, 12] <- sqrt( mean( (test2$Y - predsVecOEC2)^2 ) )
        
        rm(pMat2)
        
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
                predsVecOEC1 <- oecRidgePred(data = test2[,-c(1,2)],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 56] <- sqrt( mean( (test2$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 57] <- min(oecRidgeWS$objImp) # best objective achieved
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
        predsVecOEC1 <- oecRidgePred(data = test2[,-c(1,2)],
                                     mod = oecRidgeWS)
        
        resMat[iterNum, 54] <- sqrt( mean( (test2$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 55] <- min(oecRidgeWS$objImp) # best objective achieved
        resMat[iterNum, 64] <- OECeta
        
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
                predsVecOEC1 <- oecRidgePred(data = test2[,-c(1,2)],
                                             mod = oecRidgeWS)
                
                resMat[iterNum, 60] <- sqrt( mean( (test2$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 61] <- min(oecRidgeWS$objImp) # best objective achieved
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
        predsVecOEC1 <- oecRidgePred(data = test2[,-c(1,2)],
                                     mod = oecRidgeWS)
        
        resMat[iterNum, 58] <- sqrt( mean( (test2$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 59] <- min(oecRidgeWS$objImp) # best objective achieved
        
        #################################################################################################### 
        
        
        # ----------------------------
        # ridgeAltFixNoIntercept
        oecRidgeWS <- ridgeAlt(data = full, 
                               betaStart = warmRes$beta, 
                               wStart = warmRes$w, 
                               lambdaVec = tuneParamOEC$lambda, 
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
                               weights = Diagonal( nrow(full)),
                               orderRandom = orderRandom)
        
        
        predsMat <- as.vector(oecRidgeWS$w[1]) + X %*% as.matrix( oecRidgeWS$beta )
        
        modTune <- cv.glmnet(y = as.vector(full$Y), 
                             x = as.matrix(predsMat),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
        
        stackParam <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)
        
        if(stackCV == "zero")   stackParam <- 0
        
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
                           tuneParam = tuneParamOEC$lambda, 
                           stackParam = stackParamOECSpec, 
                           nnlsInd = TRUE,
                           stackTune = TRUE, 
                           modelStandardize = TRUE,
                           stackStandardize = TRUE,
                           glmnetOpt = glmnetOpt, # used closed form
                           xStandardize = FALSE,
                           sampSzWeight = 1, 
                           weights = Diagonal( nrow(full) ),
                           lambdaGrid = tune.grid$lambda,
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
        resMat[iterNum, 16] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        
####################################################################################################        
        ######################
        # Specialist OEC -- Country specific OEC-- initialize on country-specific stacking
        ######################
        indx1 <- which(full$Study == test_country)
        pMat2 <- matrix(ncol = itrs, nrow = nrow(test2) ) # store predictions in each column 
        pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
        for(it in 1:itrs){
            set.seed(it)
            oecRidgeWS <- ridgeAltSpec(data = full, 
                                       betaStart = warmRes$beta, 
                                       wStart = wZeroOut,
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            pMat2[,it] <- oecRidgePred(data = test2[,-c(1,2)],
                                       mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC1 <- rowMeans( pMat ) # average predictions from the different runs
            predsVecOEC2 <- rowMeans( pMat2 ) # average predictions from the different runs
            resMat[iterNum, 35] <- mean(optObj) # average objective value achieved
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC1 <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
            predsVecOEC2 <- pMat2[,bestObjIndx] # choose predictions from the one with the best objective
            resMat[iterNum, 35] <- min(optObj) # minimum objective value achieved
        }
        #}
        
        resMat[iterNum, 17] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 13] <- sqrt( mean( (test2$Y - predsVecOEC2)^2 ) )
        
        rm(pMat2)
        
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
                
                resMat[iterNum, 40] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 41] <- min(oecRidgeWS$objImp) # best objective achieved
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
        
        resMat[iterNum, 38] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 39] <- min(oecRidgeWS$objImp) # best objective achieved

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
                
                resMat[iterNum, 44] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 45] <- min(oecRidgeWS$objImp) # best objective achieved
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
        
        resMat[iterNum, 42] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 43] <- min(oecRidgeWS$objImp) # best objective achieved
        resMat[iterNum, 62] <- OECetaSpec
####################################################################################################        

        
        ######################
        # Specialist OEC -- Country specific OEC-- initialize on country-specific stacking
        ######################
        itrs2 <- 5
        
        pMat <- matrix(ncol = itrs, nrow = nrow(X) ) # store predictions in each column 
        pMat2 <- matrix(ncol = itrs, nrow = nrow(test2) ) # store predictions in each column 
        sMat <- matrix(ncol = K, nrow = nrow(X) )
        sMat2 <- matrix(ncol = K, nrow = nrow(test2) )
        
        
        for(j in 1:K){
            
            indx1 <- which(full$Study == j)
                    
            for(it in 1:itrs){
                
                
                # indices correpsonding to test country in training set for year prior
                predsMat_testCountryS <- X[indx1,] %*% betas
                predsMat_testCountryS[,j] <- 0 # zero out country thats being tested
 
                    
                modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                                         x = as.matrix(predsMat_testCountryS),
                                         alpha = 0, 
                                         lambda = tune.grid$lambda, 
                                         standardize = TRUE,
                                         intercept = TRUE,
                                         lower.limits = 0, 
                                         thresh = 1e-10
                    ) 
                    
                stackParamS <- modTune$lambda.min # replace old value with tuned one
                rm(modTune)
                    
                if(stackCV == "zero")   stackParamS <- 0
                    
                
                mod <- glmnet(y = as.vector(full$Y[indx1]), 
                              x = as.matrix(predsMat_testCountryS), 
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
                
                rm(mod)
                

                set.seed(it)
                oecRidgeWS <- ridgeAltSpec(data = full, 
                                           betaStart = warmRes$beta, 
                                           wStart = w,
                                           lambdaVec = tuneParamOEC$lambda, 
                                           mu = stackParamSpec, 
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
                
                
                pMat[,it] <- oecRidgePred(data = X[,-1],
                                          mod = oecRidgeWS)
                
                pMat2[,it] <- oecRidgePred(data = test2[,-c(1,2)],
                                           mod = oecRidgeWS)
                
                optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
                
            }
            
            
            if(avgOpt){
                sMat[,j] <- rowMeans( pMat ) # average predictions from the different runs
                sMat2[,j] <- rowMeans( pMat2 ) # average predictions from the different runs
                
            }else{
                bestObjIndx <- which.min(optObj)
                sMat[,j] <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
                sMat2[,j] <- pMat2[,bestObjIndx] # choose predictions from the one with the best objective
            }
            
        }
        
        
        modTune <- cv.glmnet(y = as.vector(full$Y), 
                             x = as.matrix(sMat),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
        
        stackParamS <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)
        
        if(stackCV == "zero")   stackParamS <- 0
        
        mod <- glmnet(y = as.vector(full$Y), 
                      x = as.matrix(sMat), 
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
           
        
        preds <- cbind(1, sMat2 ) %*% w
        
        resMat[iterNum, 23] <- sqrt( mean( (test2$Y - preds)^2 ) )
        
        
        
        resMat[iterNum, 24] <- sqrt( mean( (test2$Y - rowMeans(sMat2) )^2 ) )
        
        
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
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
                                       lambdaVec = tuneParamOEC$lambda, 
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
                                       stackTol = 1e-9,
                                       orderRandom = orderRandom)
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
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
                                   lambdaVec = tuneParamOEC$lambda, 
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
                                   stackTol = 1e-9,
                                   orderRandom = orderRandom)
        
        
        predsMat <- cbind(1, as.matrix(full[indx1,-c(1,2)] ) ) %*% as.matrix( oecRidgeWS$beta )
        
        modTune <- cv.glmnet(y = as.vector(full$Y[indx1]), 
                             x = as.matrix(predsMat),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
        
        stackParam1 <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)
        
        if(stackCV == "zero")   stackParam1 <- 0
        
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
                               lambdaVec = tuneParamOEC$lambda, 
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
                               weights = Diagonal( nrow(full)),
                               orderRandom = orderRandom)
        
        
        predsMat <-  cbind(1, as.matrix(full[,-c(1,2)] ) ) %*% as.matrix( oecRidgeWS$beta )
        
        modTune <- cv.glmnet(y = as.vector(full$Y), 
                             x = as.matrix(predsMat),
                             alpha = 0, 
                             lambda = tune.grid$lambda, 
                             standardize = TRUE,
                             intercept = TRUE,
                             lower.limits = 0, 
                             thresh = 1e-10
        ) 
        
        stackParam1 <- modTune$lambda.min # replace old value with tuned one
        rm(modTune)
        
        if(stackCV == "zero")   stackParam1 <- 0
        
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
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
            resMat[iterNum, 36] <- mean(optObj) # minimum objective value achieved
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
            resMat[iterNum, 36] <- min(optObj) # minimum objective value achieved
        }
        
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
                
                resMat[iterNum, 48] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 49] <- min(oecRidgeWS$objImp) # best objective achieved
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
        
        resMat[iterNum, 46] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 47] <- min(oecRidgeWS$objImp) # best objective achieved
        
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
                
                resMat[iterNum, 52] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
                resMat[iterNum, 53] <- min(oecRidgeWS$objImp) # best objective achieved
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
        
        resMat[iterNum, 50] <- sqrt( mean( (test$Y - predsVecOEC1)^2 ) )
        resMat[iterNum, 51] <- min(oecRidgeWS$objImp) # best objective achieved
        resMat[iterNum, 63] <- OECetaSpec0
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
                                       up = 1 / K,
                                       stackTol = 1e-9,
                                       orderRandom = orderRandom)
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
        resMat[iterNum, 32] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        
        ######################
        # Specialist OEC -- Country specific Zero Out OEC-- initialize on country-specific stacking
        ######################
        # constrain from above
        indx1 <- which(full$Study == test_country)
        
        pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
        for(it in 1:itrs){
            set.seed(it)
            oecRidgeWS <- ridgeAltSpecWindow(data = full, 
                                             betaStart = warmRes$beta, 
                                             wStart = wZeroOut,
                                             lambdaVec = tuneParamOEC$lambda, 
                                             mu = 0, #stackParamOECSpec, ## 0 because otherwise step size for stacking parameter is not right
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
                                             up = 1 / K,
                                             stackTol = 1e-9,
                                             orderRandom = orderRandom)
            
            
            pMat[,it] <- oecRidgePred(data = X_test[,-1],
                                      mod = oecRidgeWS)
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
        resMat[iterNum, 33] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        
        #######################################################
        ####################
        # PCA Specialist 
        ####################
        
        ncomp <- min( ncol(X) - 1, K ) # components are number of studies or number of features (not including intercept)
        
        warmRes <- ridgeWS(data = full, 
                           tuneParam = tuneParamOEC$lambda, 
                           stackParam = stackParamOEC, 
                           nnlsInd = TRUE,
                           stackTune = TRUE, 
                           modelStandardize = TRUE,
                           stackStandardize = TRUE,
                           glmnetOpt = glmnetOpt, # used closed form
                           xStandardize = FALSE,
                           sampSzWeight = 1, 
                           weights = Diagonal( nrow(full) ),
                           lambdaGrid = tune.grid$lambda,
                           AvgW = FALSE,
                           pcaInd = TRUE,
                           ncomp = ncomp ) # positive results achieved at FALSE 8/18/20 -- do not change
        
        sigK <- warmRes$sigK
        sigStack <- warmRes$sigStack
        
        ####
        
        indx1 <- which(full$Study == test_country)
        
        pMat <- matrix(ncol = itrs, nrow = nrow(X_test) ) # store predictions in each column 
        for(it in 1:itrs){
            set.seed(it)
            oecRidgeWS <- ridgeAltSpecPCA(data = full, 
                                          betaStart = warmRes$beta, 
                                          wStart = warmRes$w,
                                          lambdaVec = tuneParamOEC$lambda, 
                                          mu = warmRes$mu,#stackParamOEC, 
                                          Stackindx = indx1, # rows corresponding to test country
                                          nnlsInd = TRUE,
                                          tol = 0.001,
                                          objCriteria = TRUE,
                                          eta = OECeta,
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
                                          pcaMat = warmRes$pcaMat,
                                          orderRandom = orderRandom)
            
            
            pMat[,it] <- oecRidgeWS$w[1] + X_test %*% oecRidgeWS$beta %*% warmRes$pcaMat %*% oecRidgeWS$w[-1] # oecRidgePred(data = X_test[,-1],
            
            optObj[it] <- min(oecRidgeWS$objImp) # best objective achieved
            
        }
        
        
        if(avgOpt){
            predsVecOEC <- rowMeans( pMat ) # average predictions from the different runs
        }else{
            bestObjIndx <- which.min(optObj)
            predsVecOEC <- pMat[,bestObjIndx] # choose predictions from the one with the best objective
        }
        
        resMat[iterNum, 34] <- sqrt( mean( (test$Y - predsVecOEC)^2 ) )
        
        rm(warmRes, oecRidgeWS, preds, predsMat)
        
        print(paste0("complete ", iterNum))
        print(resMat[iterNum,])
        return(as.vector(resMat[iterNum,]))

    
    } # foreach loop iterating through simulation number


print("For each loop complete: make final results matrix")

resMat <- do.call(rbind, results)

print("save file")
setwd(save.folder)
write.csv(resMat, fileNm) 

colnames(resMat) <- c("merge", "avg", "stacking", "stacking_zeroOut", 
                      "merge2", "stack2", 
                      "stacking_country", "stacking_country_zeroOut",
                      "country", "oec", "oec2", "oec_test2", "oecSpec_test2",
                      "avg", "oecNoInt", "oec_country", "oec_country2", "oec_country3",
                      "oec_country4", "oec_country5",
                      "oec_countryAvg",
                      "eta", "nnet", "nnetAvg", "K", "n_k", "n_test", "N",
                      "oec_countryAvgStack", "oec_AvgStack", "oec_country0", "oec_country6",
                      "oec_countryWindow",
                      "oec_SpecPCA", "objSpec", "objZero", "objGen",
                      "oec_specAnneal_left", "oec_specAnneal_leftObj", "oec_spec_left", "oec_spec_leftObj",
                      "oec_specAnneal_right", "oec_specAnneal_rightObj", "oec_spec_right", "oec_spec_rightObj",
                      "oec_0Anneal_left", "oec_0Anneal_leftObj", "oec_0_left", "oec_0_leftObj",
                      "oec_0Anneal_right", "oec_0Anneal_rightObj", "oec_0_right", "oec_0_rightObj",
                      "oec_genAnneal_left", "oec_genAnneal_leftObj", "oec_gen_left", "oec_gen_leftObj",
                      "oec_genAnneal_right", "oec_genAnneal_rightObj", "oec_gen_right", "oec_gen_rightObj",
                      "etaGen", "etaSpec", "etaSpec0", "linearOnly")


print( resMat )

    # write files
    print("save file")
    setwd(save.folder)
    write.csv(resMat, fileNm) 
    
    
timeEnd <- Sys.time()
print(difftime(timeEnd, timeStart, units='mins'))
