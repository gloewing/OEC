
hosoCV <- function(data,
                   tune.grid,
                   hoso = "merged",
                   method = "glmnet",
                   metric = "RMSE",
                   nfolds = "K",
                   nnlsInd = TRUE,
                   OECeta = 0.5,
                   cvFolds = 10,
                   sampSzWeight = 4,
                   weights = NULL,
                   glmnetOpt = TRUE,
                   xStandardize = FALSE,
                   stackParam = NULL,
                   lambdaVec = NULL, # this is only necessary for cvOEC -- this is the initial values for it
                   specIndx = NULL, # only necessary if using cvOECSpec or cvOECSpec0
                   horizon = 10
                   ){

    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation

    library(caret)
    library(glmnet)

    # rename studies from 1:K
    num.trainStudy <- K <- length(unique(data$Study))
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    rm(indx, studyVec, indxVec)

    if( is.null(weights) )   weights <- Diagonal(nrow(data))#matrix( diag( nrow(data) ))  # Identity weights for stacking if nothing specified

    if(hoso == "merged"){
        message(hoso)
        # hold one study out on the merged dataset--used for Merged tuning

        indxList <- vector("list", length = num.trainStudy)
        HOOList <- vector("list", length = num.trainStudy)

        for(study in 1:num.trainStudy){

            indxList[[study]] <- which(data$Study != study) # indices of studies to train on
            HOOList[[study]] <- which(data$Study == study) # indices of study to hold out

        }

        control <- caret::trainControl(method="cv",
                                       number=num.trainStudy,
                                       index = indxList,
                                       indexOut = HOOList)

        tuneModel <- caret::train(Y ~.,
                                  data = data[,-1],
                                  method = method,
                                  trControl=control,
                                  tuneGrid = tune.grid,
                                  metric = metric)

        testLambda <- tuneModel$bestTune

        return(best = testLambda)

    }else if(hoso == "sse"){
        message(hoso)
        # fit individual studies for stacking
        # fit on actual study and test on all other training studies

        lambdas <- matrix(nrow = num.trainStudy, ncol = ncol(tune.grid) ) # has as many columns as there are parameters

        for(study in 1:num.trainStudy){

            indxList <- vector("list", length = 1)
            HOOList <- vector("list", length = 1)

            indxList[[1]] <- which(data$Study == study) # indices of studies to train on
            HOOList[[1]] <- which(data$Study != study) # indices of study to hold out

            control <- caret::trainControl(method="cv",
                                           number = 1,
                                           index = indxList,
                                           indexOut = HOOList)

            tuneModel <- caret::train(Y ~.,
                                      data=data[,-1],
                                      method = method,
                                      trControl=control,
                                      tuneGrid = tune.grid,
                                      metric = metric)

            lambdas[study,] <- as.numeric( tuneModel$bestTune )

        }

        lambdas <- as.data.frame(lambdas)
        colnames(lambdas) <- colnames(tune.grid)

        return(best = lambdas)

    }else if(hoso == "sseCF"){
        message(hoso)
        # closed form SSE
        # fit individual studies for stacking
        # fit on actual study and test on all other training studies MERGED TOGETHER

        source("OEC Functions.R")

        lambdaMat <- matrix(nrow = num.trainStudy, ncol = ncol(tune.grid) ) # stores best lambdas

        for(study in 1:num.trainStudy){
            rmseVec <- vector(length = nrow(tune.grid) ) # store RMSEs for current study
            indxList <- which(data$Study == study) # indices of studies to train on
            HOOList <- which(data$Study != study) # indices of study to hold out

            nk <- length(indxList) # sample size of kth study

            ## standard deviation of y to standardize the l2 hyperparameter
            ySD <- sqrt( var(data$Y[indxList]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it


            for(j in 1:nrow(tune.grid)){

                # ridge estimator
                betaHat <- ridgeEst(y = as.vector( data$Y[indxList] ),
                                    x = as.matrix( data[indxList, -c(1,2)] ),
                                    intercept = TRUE,
                                    lambda = as.numeric( tune.grid$lambda[j] ) / ySD # standardize like glmnet
                                    )

                # predictions on merged dataset of all held out studies
               # oneVec <- rep(1, length(HOOList))
                preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                rmseVec[j]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse

            }

            # chose tuning parameter with lowest RMSE for this study
            lambdaMat[study,] <- as.numeric( tune.grid[ which.min(rmseVec), ] / ySD )# standardize by study specific standard deviation of y

        }

        lambdaMat <- as.data.frame(lambdaMat)
        colnames(lambdaMat) <- colnames(tune.grid)

        # print(rbind(tune.grid$lambda, rmseMean))
        return(best = lambdaMat)

    }else if(hoso == "sseCF2"){
        message(hoso)
        # same as sseCF but do not standardize the lambdas by sd_y
        # closed form SSE
        # fit individual studies for stacking
        # fit on actual study and test on all other training studies MERGED TOGETHER

        source("OEC Functions.R")

        lambdaMat <- matrix(nrow = num.trainStudy, ncol = ncol(tune.grid) ) # stores best lambdas

        for(study in 1:num.trainStudy){
            rmseVec <- vector(length = nrow(tune.grid) ) # store RMSEs for current study
            indxList <- which(data$Study == study) # indices of studies to train on
            HOOList <- which(data$Study != study) # indices of study to hold out

            nk <- length(indxList) # sample size of kth study

            ## standard deviation of y to standardize the l2 hyperparameter
            ySD <- sqrt( var(data$Y[indxList]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it


            for(j in 1:nrow(tune.grid)){

                # ridge estimator
                betaHat <- ridgeEst(y = as.vector( data$Y[indxList] ),
                                    x = as.matrix( data[indxList, -c(1,2)] ),
                                    intercept = TRUE,
                                    lambda = as.numeric( tune.grid$lambda[j] ) # standardize like glmnet
                )

                # predictions on merged dataset of all held out studies

                preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                rmseVec[j]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse

            }

            # chose tuning parameter with lowest RMSE for this study
            lambdaMat[study,] <- as.numeric( tune.grid[ which.min(rmseVec), ]  )# standardize by study specific standard deviation of y

        }

        lambdaMat <- as.data.frame(lambdaMat)
        colnames(lambdaMat) <- colnames(tune.grid)

        # print(rbind(tune.grid$lambda, rmseMean))
        return(best = lambdaMat)

    }else if (hoso == "stacking"){
        message(hoso)
        # use sse function and all of the same tuning parameters

        # set number of folds to number of training studies unless
        # specified as another number
        if(nfolds == "K" || nfolds > K)       nfolds <- K

        resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets

        for(study in 1:nfolds){
            HOOList <- which(data$Study == study) # indices of study to hold out
            indxList <- which(data$Study != study) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){

                print(paste0("study ", study, "_tune_", tune))

                sse.mod <- sse(formula = Y ~.,
                               data = data[indxList,], # current studies
                               sim.mets = FALSE,
                               model = FALSE,
                               ssl.method = list(method),
                               ssl.tuneGrid = list( tune.grid[tune,] ) # current hyperparameters
                )

                preds <- studyStrap.predict(sse.mod, data[HOOList,-1])
                resMat[study, tune] <- sqrt(mean( (preds[,2] - data$Y[HOOList])^2  )) # stacking predictions

        }

            }

        resMean <- colMeans(resMat)

        testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE

        return(best = testLambda)

    }else if (hoso == "oec"){
        message(hoso)
        source("OEC Functions.R")

        # set number of folds to number of training studies unless
        # specified as another number

        if(nfolds == "K" || nfolds > K)       nfolds <- K

        resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets

        for(study in 1:nfolds){

            # changed 9/11/20 to make indices match
            HOOList <- which(data$Study == study) # indices of study to hold out
            indxList <- which(data$Study != study) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            # if no stack Param provided then tune it, otherwise keep it fixed in ridgeWS
            if(is.null(stackParam)){

                stackTune <- TRUE
                stackParam <- 0

            }else{

                # use the stackParam provided as an argument into the fucntion hosoCV
                stackTune <- FALSE

            }

            for(tune in 1:nrow(tune.grid) ){

                print(paste0("study ", study, "_tune_", tune))
                set.seed(1)
                warmRes <- ridgeWS(data[indxList,], # current studies
                                   tuneParam = tune.grid[tune, 2], # current lambda associated with L2 penalty
                                   stackParam = stackParam,
                                   nnlsInd = nnlsInd,
                                   stackTune = stackTune,
                                   modelStandardize = TRUE,
                                   stackStandardize = TRUE,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack
                stackParam <- warmRes$mu
                tuneParam <- rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAlt(data = data[indxList,], # current studies
                                       betaStart = warmRes$beta,
                                       wStart = warmRes$w,
                                       lambdaVec = tuneParam,
                                       mu = stackParam,
                                       nnlsInd = nnlsInd,
                                       tol = 0.01,
                                       objCriteria = TRUE,
                                       eta = OECeta,
                                       projs = 0,
                                       Avg = FALSE,
                                       wUpdate = "glment",
                                       standardize = TRUE,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       sampSzWeight = sampSzWeight,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        resMean <- colMeans(resMat)

        testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE

        return(best = testLambda)

    }else if (hoso == "cv"){
        message(hoso)
        # 10 fold CV within studies

        source("OEC Functions.R")

        # store optimal parameters for each study
        paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
        colnames(paramMat) <- colnames(tune.grid)

        for(study in 1:K){

            # standard cross validation within each study
            control <- caret::trainControl(method="cv",
                                           number = cvFolds)

            # changed on 9/11/20 to make indices match
            # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out

            indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match

            tuneModel <- caret::train(Y ~.,
                                      data=data[indx, -1],
                                      method = method,
                                      trControl= control,
                                      tuneGrid = tune.grid,
                                      metric = metric)

            paramMat[study,] <- tuneModel$bestTune

        }
        }else if (hoso == "cvCF"){
            message(hoso)
            # 10 fold CV within studies
            # used closed form ridge

            source("OEC Functions.R")
            library(caret)
            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)

            for(study in 1:K){

                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match
                foldsList <- createFolds(indx, k = nfolds) # break up into F fold CV

                rmseMat <- matrix( nrow = nfolds, ncol = nrow(tune.grid) )

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it

                for(j in 1:nrow(tune.grid)){


                    for(f in 1:nfolds){
                        HOOList <- indx[ foldsList[[f]] ] # indices of study to hold out
                        indxList <- indx[ -foldsList[[f]] ]#do.call(c, foldsList[-f]) ] # concatenate indices of studies to train on

                        # ridge estimator
                        betaHat <- ridgeEst(y = as.vector( data$Y[indxList] ),
                                            x = as.matrix( data[indxList, -c(1,2)] ),
                                            intercept = TRUE,
                                            lambda = as.numeric( tune.grid$lambda[j] / ySD) # standardize like glmnet
                        )

                        # predictions on merged dataset of all held out studies
                        # oneVec <- rep(1, length(HOOList))
                        preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                        rmseMat[f, j]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

            }

        }else if (hoso == "zero"){
            message(hoso)
            # matrix of zeros
            paramMat <- as.data.frame( matrix(0, ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

        }else if (hoso == "cvCFTS"){
            message(hoso)
            # TS version
            # 10 fold CV within studies
            # used closed form ridge

            source("OEC Functions.R")
            library(caret)
            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)

            for(study in 1:K){

                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match

                horizon2 <- min(horizon,
                                length(indx) - round(0.5 * length(indx))
                )

                foldsList <- createTimeSlices(indx,
                                              initialWindow = round(0.5 * length(indx)), # minimum training set is half of total observations
                                              horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                              fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
                                              )

                # remove some so we dont do a HO-observation out CV
                len <- length(foldsList$train)
                div <- round(len / nfolds)

                if(div > 0){

                    seqVec <- seq(from = 1, to = len, by = div) # sequence that jumps to ensure we have roughly k folds
                    # always make last index the last one
                    seqVec[length(seqVec)] <- len
                    len <- length(seqVec)
                }else{
                    seqVec <- 1
                    len <- 1
                    foldsList$test <- unique( do.call(c, foldsList$test) ) # collapse into 1
                }


                rmseMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

                # remove some so we do not have so many
                foldsList$test <- foldsList$test[seqVec]
                foldsList$train <- foldsList$train[seqVec]

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1) / nk) #  --- standardized the way glmnet does it

                for(j in 1:nrow(tune.grid)){


                    for(f in 1:len){
                        HOOList <- indx[ foldsList$test[[f]] ]# indices of study to hold out
                        indxList <- indx[ foldsList$train[[f]] ]# concatenate indices of studies to train on

                        # ridge estimator
                        betaHat <- ridgeEst(y = as.vector( data$Y[indxList] ),
                                            x = as.matrix( data[indxList, -c(1,2)] ),
                                            intercept = TRUE,
                                            lambda = as.numeric( tune.grid$lambda[j] / ySD) # standardize like glmnet
                                            )

                        # predictions on merged dataset of all held out studies
                        # oneVec <- rep(1, length(HOOList))
                        preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat )
                        rmseMat[f, j]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

            }

        }else if (hoso == "cvTS"){
            message(hoso)
            # TS version
            # 10 fold CV within studies
            # used glmnet

            source("OEC Functions.R")
            library(caret)
            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)

            for(study in 1:K){

                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match

                horizon2 <- min(horizon,
                                length(indx) - round(0.5 * length(indx))
                )


                foldsList <- createTimeSlices(indx,
                                              initialWindow = round(0.5 * length(indx)), # minimum training set is half of total observations
                                              horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                              fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
                )

                # remove some so we dont do a HO-observation out CV
                len <- length(foldsList$train)
                div <- round(len / nfolds)

                if(div > 0){

                    seqVec <- seq(from = 1, to = len, by = div) # sequence that jumps to ensure we have roughly k folds
                    # always make last index the last one
                    seqVec[length(seqVec)] <- len
                    len <- length(seqVec)
                }else{
                    seqVec <- 1
                    len <- 1
                    foldsList$test <- unique( do.call(c, foldsList$test) ) # collapse into 1
                }


                rmseMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

                # remove some so we do not have so many
                foldsList$test <- foldsList$test[seqVec]
                foldsList$train <- foldsList$train[seqVec]

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1) / nk) #  --- standardized the way glmnet does it

                for(j in 1:nrow(tune.grid)){


                    for(f in 1:len){
                        HOOList <- indx[ foldsList$test[[f]] ]# indices of study to hold out
                        indxList <- indx[ foldsList$train[[f]] ]# concatenate indices of studies to train on

                        # ridge estimator
                        mod <- glmnet(y = as.vector( data$Y[indxList] ),
                                           x = as.matrix( data[indxList, -c(1,2)] ),
                                           intercept = TRUE,
                                           alpha = tune.grid$alpha[j],
                                           lambda = tune.grid$lambda[j],
                                           standardize = TRUE
                                          )

                        betaHat <- as.vector( coef(mod) )

                        # predictions on merged dataset of all held out studies
                        # oneVec <- rep(1, length(HOOList))
                        preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat )
                        rmseMat[f, j]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

            }

        }else if (hoso == "cvOEC"){
            # OEC for each study separately holding the constant
            message(hoso)
            # 10 fold CV within studies
            # used closed form ridge

            source("OEC Functions.R")

            library(caret)

            if(is.null(stackParam)){

                stackTune <- TRUE
                stackParam <- 0

            }else{

                # use the stackParam provided as an argument into the fucntion hosoCV
                stackTune <- FALSE

            }

            sslLambdas <- lambdaVec # initialize lambdas with what was fed to the function

            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)
            allRows <- 1:nrow(data)

            for(study in 1:K){
                message(paste("Tuning Study:", study))
                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match
                foldsList <- createFolds(indx, k = nfolds) # break up into F fold CV

                rmseMat <- matrix( nrow = nfolds, ncol = nrow(tune.grid) )

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it

                for(tune in 1:nrow(tune.grid)){

                    sslLambdas[study] <- tune.grid[tune, 2] # alter study specific

                    for(f in 1:nfolds){
                        HOOList <- indx[ foldsList[[f]] ] # indices of study to hold out
                        indxList <- allRows[ -HOOList ] #indx[ -foldsList[[f]] ]#do.call(c, foldsList[-f]) ] # concatenate indices of studies to train on

                        wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])


                        ######################
                        warmRes <- ridgeWS(data[indxList,], # current studies
                                           tuneParam = sslLambdas, # current lambda associated with L2 penalty of study
                                           stackParam = stackParam,
                                           nnlsInd = nnlsInd,
                                           stackTune = stackTune, # tune this accordingly CONSIDER CHANGING THIS IF YOU RETUNE and loop through 9/18/20
                                           modelStandardize = TRUE,
                                           stackStandardize = TRUE,
                                           glmnetOpt = glmnetOpt,
                                           xStandardize = xStandardize,
                                           weights = wgt)

                        sigK <- warmRes$sigK
                        sigStack <- warmRes$sigStack
                        stackParam <- warmRes$mu
                        tuneParam <- rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies

                        set.seed(1)

                        oecRidgeWS <- ridgeAlt(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = sslLambdas,
                                               mu = stackParam,
                                               nnlsInd = nnlsInd,
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               sampSzWeight = sampSzWeight,
                                               weights = wgt)

                        predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                    mod = oecRidgeWS)

                        rmseMat[f, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions



                        #######################
                        # predictions on merged dataset of all held out studies
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

                # update the current vector of lambdas
                sslLambdas[study] <- tune.grid$lambda[bestLambda] / ySD

            }

            # save the final version of lambdas
            paramMat[,2] <- sslLambdas

        }else if (hoso == "cvOECSpec"){
            # OEC for each study separately holding the constant
            message(hoso)
            # 10 fold CV within studies
            # used closed form ridge

            source("OEC Functions.R")

            library(caret)

            if(is.null(stackParam)){

                stackTune <- TRUE
                stackParam <- 0

            }else{

                # use the stackParam provided as an argument into the fucntion hosoCV
                stackTune <- FALSE

            }

            sslLambdas <- lambdaVec # initialize lambdas with what was fed to the function

            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)
            allRows <- 1:nrow(data)

            for(study in 1:K){

                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match
                foldsList <- createFolds(indx, k = nfolds) # break up into F fold CV

                rmseMat <- matrix( nrow = nfolds, ncol = nrow(tune.grid) )

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it

                for(tune in 1:nrow(tune.grid)){

                    sslLambdas[study] <- tune.grid[tune, 2] # alter study specific

                    for(f in 1:nfolds){
                        HOOList <- indx[ foldsList[[f]] ] # indices of study to hold out
                        indxList <- allRows[ -HOOList ] #indx[ -foldsList[[f]] ]#do.call(c, foldsList[-f]) ] # concatenate indices of studies to train on
                        stackIndx <- which(data$Study[indxList] == specIndx) # rows in current training fold that are specialist study
                        wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

                        ######################
                        warmRes <- ridgeWS(data[indxList,], # current studies
                                           tuneParam = sslLambdas, # current lambda associated with L2 penalty of study
                                           stackParam = stackParam,
                                           nnlsInd = nnlsInd,
                                           stackTune = stackTune, # tune this accordingly CONSIDER CHANGING THIS IF YOU RETUNE and loop through 9/18/20
                                           modelStandardize = TRUE,
                                           stackStandardize = TRUE,
                                           glmnetOpt = glmnetOpt,
                                           xStandardize = xStandardize,
                                           weights = wgt)

                        sigK <- warmRes$sigK
                        sigStack <- warmRes$sigStack
                        stackParam <- warmRes$mu
                        tuneParam <- rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies

                        set.seed(1)

                        oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = sslLambdas,
                                               Stackindx = stackIndx,
                                               mu = stackParam,
                                               nnlsInd = nnlsInd,
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               projs = 0,
                                               Avg = FALSE,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               sampSzWeight = sampSzWeight,
                                               weights = wgt)

                        predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                    mod = oecRidgeWS)

                        rmseMat[f, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions



                        #######################
                        # predictions on merged dataset of all held out studies
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

                # update the current vector of lambdas
                sslLambdas[study] <- tune.grid$lambda[bestLambda] / ySD

            }

            # save the final version of lambdas
            paramMat[,2] <- sslLambdas

        }else if (hoso == "cvOECSpec0"){
            # OEC for each study separately holding the constant
            message(hoso)
            # 10 fold CV within studies
            # used closed form ridge

            source("OEC Functions.R")

            library(caret)

            if(is.null(stackParam)){

                stackTune <- TRUE
                stackParam <- 0

            }else{

                # use the stackParam provided as an argument into the fucntion hosoCV
                stackTune <- FALSE

            }

            sslLambdas <- lambdaVec # initialize lambdas with what was fed to the function

            # store optimal parameters for each study
            paramMat <- as.data.frame( matrix(ncol = ncol(tune.grid), nrow = K) )
            colnames(paramMat) <- colnames(tune.grid)

            # determine number of folds
            nVec <- table(data$Study)
            nfolds <- min( c(10, nVec) ) # do 10 fold CV or do F - fold CV using the minimum number of obs in any study
            trainStudies <- unique(data$Study)
            allRows <- 1:nrow(data)

            for(study in 1:K){

                # standard cross validation within each study

                indx <- which(data$Study == study) # indices of study to hold out -- changed on 9/11/20 to make indices match
                foldsList <- createFolds(indx, k = nfolds) # break up into F fold CV

                rmseMat <- matrix( nrow = nfolds, ncol = nrow(tune.grid) )

                nk <- length(indx)
                ySD <- sqrt( var(data$Y[indx]) * (nk - 1)/ nk) #  --- standardized the way glmnet does it

                for(tune in 1:nrow(tune.grid)){

                    sslLambdas[study] <- tune.grid[tune, 2] # alter study specific

                    for(f in 1:nfolds){
                        HOOList <- indx[ foldsList[[f]] ] # indices of study to hold out
                        indxList <- allRows[ -HOOList ] #indx[ -foldsList[[f]] ]#do.call(c, foldsList[-f]) ] # concatenate indices of studies to train on
                        stackIndx <- which(data$Study[indxList] == specIndx) # rows in current training fold that are specialist study
                        wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

                        ######################
                        warmRes <- ridgeWS(data[indxList,], # current studies
                                           tuneParam = sslLambdas, # current lambda associated with L2 penalty of study
                                           stackParam = stackParam,
                                           nnlsInd = nnlsInd,
                                           stackTune = stackTune, # tune this accordingly CONSIDER CHANGING THIS IF YOU RETUNE and loop through 9/18/20
                                           modelStandardize = TRUE,
                                           stackStandardize = TRUE,
                                           glmnetOpt = glmnetOpt,
                                           xStandardize = xStandardize,
                                           zeroOut = TRUE,
                                           weights = wgt)

                        sigK <- warmRes$sigK
                        sigStack <- warmRes$sigStack
                        stackParam <- warmRes$mu
                        tuneParam <- rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies

                        set.seed(1)

                        oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                                   betaStart = warmRes$beta,
                                                   wStart = warmRes$w,
                                                   lambdaVec = sslLambdas,
                                                   Stackindx = stackIndx,
                                                   mu = stackParam,
                                                   nnlsInd = nnlsInd,
                                                   tol = 0.01,
                                                   objCriteria = TRUE,
                                                   eta = OECeta,
                                                   projs = 0,
                                                   Avg = FALSE,
                                                   wUpdate = "glment",
                                                   standardize = TRUE,
                                                   sigK = sigK,
                                                   sigStack = sigStack,
                                                   sampSzWeight = sampSzWeight,
                                                   weights = wgt)

                        predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                    mod = oecRidgeWS)

                        rmseMat[f, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions



                        #######################
                        # predictions on merged dataset of all held out studies
                    }

                }

                rmseMean <- colMeans(rmseMat)
                bestLambda <- which.min(rmseMean)

                paramMat[study,] <- tune.grid[bestLambda,] / ySD

                # update the current vector of lambdas
                sslLambdas[study] <- tune.grid$lambda[bestLambda] / ySD

            }

            # save the final version of lambdas
            paramMat[,2] <- sslLambdas

        }

        return(best = paramMat)


    }


