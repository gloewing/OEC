# OEC Functions
# Gabe Loewinger

######################
# Constrained Stacking
######################
# stacking with upper and lower limits on coefficient estimates
stackLim <- function(y,
                     x,
                     lambda = 0,
                     low = 0,
                     up = Inf,
                     w0 = NULL,
                     tol = 1e-12,
                     u1 = nrow(x),
                     u2 = 1,
                     weights = NULL,
                     intercept = TRUE){

    if( low <= 0  ){
        # added in 9/2/20 -- criteria to use glmnet

        # if you can use glmnet, use glmnet
        # if low > 0 then cannot use glmnet
        # if lambda > 0 then the below is not converging

        mod <- glmnet(y = y,
                      x = x,
                      lower.limits = low,
                      upper.limits = up,
                      weights = weights,
                      intercept = intercept,
                      lambda = lambda,
                      alpha = 0) # ridge

        w0 <- as.vector(coef(mod))

    }else{

        if(intercept){

            if(is.null(w0)){
                w0 <- coef(lm(y ~ x))
                w0 <- ifelse(is.na(w0), 0, w0)
            }

            x <- cbind(1, x) # intercept

            w0 <- c(w0[1], sapply(w0[-1], function(x) max(low, x) ) )

        }else{
            # if no intercept

            if(is.null(w0)){

                w0 <- coef(  lm(y ~ x - 1)  )
                w0 <- ifelse(is.na(w0), 0, w0)

            }

            w0 <- sapply(w0, function(x) max(low, x) )
        }


        # observation weights
        if(is.null(weights)){
            XX <- t(as.matrix(x)) %*% as.matrix(x)
            XY <- t(as.matrix(x)) %*% as.vector(y)
        }else{
            XX <- t(as.matrix(x)) %*% as.matrix(x * weights)
            XY <- t(as.matrix(x * weights)) %*% as.vector(y)

        }

        step <- 1 / sum(x^2)
        eps <- 1 + tol
        w <- w0 # initialize

        if(intercept){

            if(lambda == 0){

                while(eps > tol){
                    w0 <- w - step * (-XY + XX %*% w)

                    # project
                    w0[-1] <- sapply(w0[-1], function(x) max(low, x) ) # project onto non-neg orthant

                    eps <- sum(  (w - w0)^2 )
                    w <- w0

                }

            }else{

                lambda <- lambda * u2
                ind <- lambda * c(0, rep(1, length(w[-1])) ) # do not penalize intercept

                while(eps > tol){
                    # if zero lambda and intercept
                    w0 <- w - ( step * (-XY + XX %*% w) + ind * w )

                    # project
                    w0[-1] <- sapply(w0[-1], function(x) max(low, x) ) # project onto non-neg orthant


                    eps <- sum(  (w - w0)^2 )
                    w <- w0
                }

            }
        }else{
            # if no intercept

            if(lambda == 0){
                while(eps > tol){
                    w0 <- w - step * (-XY + XX %*% w)

                    # project
                    w0 <- sapply(w0, function(x) max(low, x) ) # project onto non-neg orthant

                    eps <- sum(  (w - w0)^2 )
                    w <- w0

                }

            }else{
                # if nonzero lambda and no intercept
                lambda <- lambda * u2

                while(eps > tol){
                    w0 <- w - ( step * (-XY + XX %*% w) + lambda * w )

                    # project
                    w0 <- sapply(w0, function(x) max(low, x) ) # project onto non-neg orthant


                    eps <- sum(  (w - w0)^2 )
                    w <- w0
                }

            }
        }

    }

    return(w0)

}

######################
# Constrained Stacking
######################
# stacking with upper limits on coefficient estimates for TEST Training set
stackLimUp <- function(y,
                     x,
                     lambda = 0,
                     low = 0,
                     up = Inf,
                     w0 = NULL,
                     tol = 1e-12,
                     u1 = nrow(x),
                     u2 = 1,
                     weights = NULL,
                     intercept = TRUE,
                     testIndx){


        if(intercept){

            if(is.null(w0)){
                w0 <- coef(lm(y ~ x))
                w0 <- ifelse(is.na(w0), 0, w0)
            }

            x <- cbind(1, x) # intercept

            w0 <- c(w0[1], sapply(w0[-1], function(x) max(low, x) ) )
            w0[testIndx + 1] <- min(w0[testIndx + 1], up) # constrain test coefficient index

        }else{
            # if no intercept

            if(is.null(w0)){

                w0 <- coef(  lm(y ~ x - 1)  )
                w0 <- ifelse(is.na(w0), 0, w0)

            }

            w0 <- sapply(w0, function(x) max(low, x) )
            w0[testIndx] <- min(w0[testIndx + 1], up) # constrain test coefficient index
        }


        # observation weights
        if(is.null(weights)){
            XX <- t(as.matrix(x)) %*% as.matrix(x)
            XY <- t(as.matrix(x)) %*% as.vector(y)
        }else{
            XX <- t(as.matrix(x)) %*% as.matrix(x * weights)
            XY <- t(as.matrix(x * weights)) %*% as.vector(y)

        }

        step <- 1 / sum(x^2)
        eps <- 1 + tol
        w <- w0 # initialize

        if(intercept){

            if(lambda == 0){

                while(eps > tol){
                    w0 <- w - step * (-XY + XX %*% w)

                    # project
                    w0[-1] <- sapply(w0[-1], function(x) max(low, x) ) # project onto non-neg orthant
                    w0[testIndx + 1] <- min(w0[testIndx + 1], up) # constrain test coefficient index

                    eps <- sum(  (w - w0)^2 )
                    w <- w0

                }

            }else{

                lambda <- lambda * u2
                ind <- lambda * c(0, rep(1, length(w[-1])) ) # do not penalize intercept

                while(eps > tol){
                    # if zero lambda and intercept
                    w0 <- w - ( step * (-XY + XX %*% w) + ind * w )

                    # project
                    w0[-1] <- sapply(w0[-1], function(x) max(low, x) ) # project onto non-neg orthant
                    w0[testIndx + 1] <- min(w0[testIndx + 1], up) # constrain test coefficient index

                    eps <- sum(  (w - w0)^2 )
                    w <- w0
                }

            }
        }else{
            # if no intercept

            if(lambda == 0){
                while(eps > tol){
                    w0 <- w - step * (-XY + XX %*% w)

                    # project
                    w0 <- sapply(w0, function(x) max(low, x) ) # project onto non-neg orthant
                    w0[testIndx] <- min(w0[testIndx], up) # constrain test coefficient index

                    eps <- sum(  (w - w0)^2 )
                    w <- w0

                }

            }else{
                # if nonzero lambda and no intercept
                lambda <- lambda * u2

                while(eps > tol){
                    w0 <- w - ( step * (-XY + XX %*% w) + lambda * w )

                    # project
                    w0 <- sapply(w0, function(x) max(low, x) ) # project onto non-neg orthant
                    w0[testIndx] <- min(w0[testIndx], up) # constrain test coefficient index

                    eps <- sum(  (w - w0)^2 )
                    w <- w0
                }

            }
        }


    return(w0)

}


########################
# Ridge OEC Functions
########################

ridgeWS <- function(data,
                    tuneParam = 0,
                    stackParam = 0,
                    nnlsInd = TRUE,
                    stackInt = TRUE,
                    stackTune = TRUE,
                    modelStandardize = FALSE,
                    stackStandardize = FALSE,
                    xStandardize = FALSE,
                    sampSzWeight = 4,
                    weights = NULL,
                    lambdaGrid = NULL,
                    glmnetOpt = FALSE,
                    zeroOut = FALSE,
                    AvgW = FALSE,
                    pcaInd = FALSE,
                    ncomp = ncol(data) - 1){


    # tune lamda is a vector of lambdas to tune the stacking regression for cv.glmnet
    # nnlsInd indicates whether stacking weights have nonnegativity constraint
    library(glmnet)

    p <- ncol(data) - 2 # doesnt include intercept
    K <- length(unique(data$Study))
    sigK <- sigStack <- vector(length = K)

    if( is.null(weights) )      weights <- matrix( diag( nrow(data) ) ) # make weights identity if nothing provided

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(length(tuneParam) == 1){
        # if one number is given for ridge penalty hyperparameter
        # then assume it is the same for each study and repeat it K times to make vector of
        # study specific hyperparameters
        tuneParam <- rep(tuneParam, K)
    }

    if(nnlsInd){
        # check whether there is a nonnegativity constraint on stacking weights
        # no constraint on the intercept
        wLimits <- 0

    }else{
        wLimits <- -Inf
    }

    if(!glmnetOpt){

        diagMat <- diag( ncol(data ) - 1 )
        diagMat[1,1] <- 0 # so do not penalize intercept


    }

    ############################################################
    ################################
    # Covaraite Standardization
    ################################
    meanMat <- matrix(nrow = ncol(data) - 2, ncol = K)
    sigmaMat <- matrix(nrow = ncol(data) - 2, ncol = K)

    if(xStandardize){
        y_SD.vec <- vector(length = K)
        # if standardize covariates with 1 / n_k formula


        for(studyNum in 1:K){

            Sindx <- which( data$Study == studyNum    )  # index of kth study observations
            nk <- length( Sindx ) # study sample size

            # sd (MLE) of covaraites
            sd_x <- sqrt(
                        apply( data[Sindx, -c(1,2)], 2, var )  *  (nk - 1) / nk
                        )
            sd_x <- ifelse(sd_x == 0, 1, sd_x) # replace sd = 0 to 1 so it divides by 1

            # mle standard deviation of outcome
            y_SD.vec[studyNum] <- sqrt( var(data$Y[Sindx]) * (nk - 1) / nk) # glmnet formula for study specific MLE std. dev

            cov_ind <- 3:ncol(data) # columns of covariates (first two columns are outcome and study)
            mean_x <- colMeans(data[Sindx, cov_ind]) # study specific covaraite means

            meanMat[,studyNum] <- mean_x
            sigmaMat[,studyNum] <- sd_x

            # scale covaraites -- iterate through each covariate
            for(cov_i in cov_ind){

                data[Sindx, cov_i] <- ( data[Sindx, cov_i] - mean_x[cov_i - 2] ) / sd_x[cov_i - 2]

            }

        }

    }

    ############################################################

    ens_betaMat <- matrix(ncol = K, nrow = ncol(data) - 1)

    # identity matrix for closed form Ridge estimator
    if(!glmnetOpt){
        diagMat <- diag( ncol(data) - 1 ) # p + 1 for intercept
        diagMat[1,1] <- 0 # so do not penalize intercept
    }

    for(studyNum in 1:K){

        # if use glmnet for optimization
        if(glmnetOpt){

            Sindx <- which( data$Study == studyNum    )

            # USE EXACT = TRUE

            cfs <- coef(
                mod <- glmnet(y = as.vector(data[Sindx,2]),
                              x = as.matrix(data[Sindx,-c(1,2)]),
                              alpha = 0,          # no lasso penalty
                              lambda = tuneParam[studyNum], # ridge penalty
                              standardize = modelStandardize,
                              intercept = TRUE
                              ), #,
                              #thresh = 1e-20),
                exact = TRUE,
                y = as.vector(data[Sindx,2]),
                x = as.matrix(data[Sindx,-c(1,2)])
            )

        }else{

            # if use clsoed form expression
            Sindx <- which( data$Study == studyNum    )
            nS <- length(Sindx) # n_k

            # ridge estimator
            cfs <- solve( t( as.matrix(cbind(1, data[Sindx,-c(1,2)]) ) ) %*%
                                as.matrix(cbind(1, data[Sindx,-c(1,2)]) ) + nS *
                                    tuneParam[studyNum] * diagMat ) %*%
                                        (t(as.matrix(cbind(1, data[Sindx,-c(1,2)]) )) %*%
                                            as.vector(data$Y[Sindx]) ) # original y (not scaled)


        }

        # variance of residual
        sigK[studyNum] <- var(  data$Y[Sindx] - as.matrix(cbind(1, data[Sindx,-c(1,2)]) ) %*% as.vector(cfs)   )


        ens_betaMat[, studyNum] <- as.vector(cfs)

    }

    if(pcaInd){
        if(stackInt == TRUE){
            pcs <- prcomp( as.matrix(cbind(1, data[,-c(1,2)])) %*% ens_betaMat, center = FALSE)
        }else{
            pcs <- prcomp( as.matrix(cbind(1, data[,-c(1,2)])) %*% ens_betaMat, center = FALSE )
        }

        pcs <- as.matrix( pcs$rotation[,1:ncomp] )# only keep a certain number of components

    }else{

        pcs <- diag( ncol( ens_betaMat )) # identity matrix
    }




    ########################
    # Stacking - Traditional
    ########################
    # print("stacking")

    # if an intercept is included in the stacking regression
    if(stackInt == TRUE){


        predsMat <- as.matrix( cbind(1, data[,-c(1,2)] ) )  %*% as.matrix( ens_betaMat ) %*% pcs

        if(zeroOut){

            for(studyNum in 1:K){
                sIndex <- which(data$Study == studyNum)
                predsMat[sIndex, studyNum ] <- 0 # zero out
            }

        }

        if(stackTune){
            # tune stack parameter
            modTune <- cv.glmnet(y = as.vector(data[,2]),
                                 x = as.matrix(predsMat),
                                alpha = 0, # no lasso penalty
                                #foldid = data[,1],
                                lambda = lambdaGrid, # ridge penalty
                                standardize = stackStandardize,
                                intercept = TRUE,
                                lower.limits = wLimits,
                                weights = diag(weights) )

            stackParam <- modTune$lambda.min # replace old value with tuned one
        }

        mod <- glmnet(y = as.vector(data[,2]),
                      x = as.matrix(predsMat),
                      alpha = 0,          # no lasso penalty
                      lambda = stackParam, # ridge penalty
                      standardize = stackStandardize,
                      intercept = TRUE,
                      lower.limits = wLimits,
                      weights = diag(weights),
                      thresh = 1e-10 )

        w <- as.vector( mod$beta ) # stacking weights
        w <- c( as.vector( mod$a0 ), w ) # add intecept
        w <- as.vector( coef( mod ,
                        exact = TRUE,
                        y = as.vector(data[,2]),
                        x = as.matrix(predsMat)
                        )
        )

    }else{

        if(stackTune){
            # tune stack parameter
            modTune <- cv.glmnet(y = as.vector(data[,2]),
                                 x = as.matrix(predsMat),
                                 alpha = 0, # no lasso penalty
                                 #foldid = data[,1],
                                 lambda = lambdaGrid, # ridge penalty
                                 standardize = stackStandardize,
                                 intercept = TRUE,
                                 lower.limits = wLimits,
                                 weights = diag(weights) )

            stackParam <- modTune$lambda.min # replace old value with tuned one
        }

        predsMat <- as.matrix( cbind(1, data[,-c(1,2)] ) ) %*% as.matrix( ens_betaMat ) %*% pcs$rotation

        if(zeroOut){

            for(studyNum in 1:K){
                sIndex <- which(data$Study == studyNum)
                predsMat[sIndex, studyNum ] <- 0 # zero out
            }
        }

        mod <- glmnet(y = as.vector(data[,2]),
                      x = as.matrix(predsMat),
                      alpha = 0,          # no lasso penalty
                      lambda = stackParam, # tuneParam, # ridge penalty
                      standardize = stackStandardize,
                      intercept = FALSE,
                      lower.limits = wLimits,
                      weights = diag(weights),
                      thresh = 1e-10 )


        w <- as.vector( coef( mod,
                              exact = TRUE,
                              y = as.vector(data[,2]),
                              x = as.matrix(predsMat))
                        ) # stacking weights

    }

    if(AvgW){
        w <- rep(1 / K, K)
        wT0 <- (
                    sum(data$Y) - sum( as.matrix( cbind(1, data[,-c(1,2)] ) ) %*%
                        as.matrix( ens_betaMat ) %*% w )
                ) / nrow(data)
        w <- c( wT0, w )
    }

    # variance of residuals for stacking
    fitted_values <- predict(mod, as.matrix(predsMat)) # fitted values from stacking regression
    residS <- data$Y - fitted_values # residual from stacking regression

    for(studyNum in 1:K){

        Sindx <- which( data$Study == studyNum    )  # indices of current study

        # variance of residual
        sigStack[studyNum] <- var(  residS[Sindx]  )

    }

    # objective value
    if(pcaInd){
        obj <- objOECPCA(data = data,
                      beta = ens_betaMat,
                      w = w,
                      mu = stackParam,
                      lambdaVec = tuneParam,
                      stackInt = TRUE,
                      eta = 0.5,
                      sampSzWeight = sampSzWeight,
                      weights = weights,
                      sigK = sigK,
                      sigStack = sigStack,
                      pcaMat = pcs)
    }else{
        obj <- objOEC(data = data,
                      beta = ens_betaMat,
                      w = w,
                      mu = stackParam,
                      lambdaVec = tuneParam,
                      stackInt = TRUE,
                      eta = 0.5,
                      sampSzWeight = sampSzWeight,
                      weights = weights,
                      sigK = sigK,
                      sigStack = sigStack)
    }


    # project w onto pcas
    return(list(beta = ens_betaMat, w = w, obj = obj,
                mu = stackParam, sigK = sigK,
                sigStack = sigStack,
                meanMat = meanMat,
                sigmaMat = sigmaMat,
                pcaMat = pcs)
                )
}




#############
# NewWR
#############
ridgeAltFast <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         mu,
                         nnlsInd = TRUE,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL){

    # eta is the parameter that determines convex combination of the losses
    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(nnlsInd){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # No nonegative least squares for glmnet paramter
        lowLim <- -Inf
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            w <- sapply(wStart[-1], function(x) max(0, x) )
        }else{

            w <- wStart[-1]

        }

    }



    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N
        weights <- Diagonal( x = 1 / sigStack[data$Study] )
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)

    SSLindx <- c()
    Stackindx <- c()

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )

    y_sum <- sum(y[Stackindx])
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    objVec <- objImp <- obj0 <- objOEC(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = TRUE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- matrix(ncol = K, nrow = p )

    for(k in 1:K){

        indx <- which(Study == k)
        lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list <- t(X[indx,]) %*% X[indx,]
        Xy_k_list <- t(X[indx,])  %*% y[indx]

        inv_list[[k]] <- as.matrix( solve(
            eta * u1 * w[k]^2 * XX +
                (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )
        ) )


        mat[,k] <- as.vector(inv_list[[k]] %*% (
            ( eta * u1 * w[k] * ( Xy - w0 * X_rowSums ) ) +
                (1 - eta) * u3[k] * Xy_k_list
        ) )

        rm(XX_k_list, Xy_k_list)

    }

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize
    }

    if(dataSplit == 1)    Stackindx <- NULL # for objective evaluation below
    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)
    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){


            betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + mat[,k]

            obj_k <- objOEC(data = data,
                            beta = betaTemp,
                            w = c(w0, w),
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = TRUE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                mod <- glmnet(y = y,
                              x = as.matrix(X %*% beta),
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = standardize,
                              intercept = TRUE,
                              weights = weights,
                              thresh = 1e-10)

                # w coefficients
                wT <-  as.vector(mod$beta)
                wT0 <- as.vector(mod$a0)

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, wT),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, w),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]


    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}
#############
# NewWR
#############

ridgeAlt <- function(data,
                          betaStart,
                          wStart,
                          lambdaVec,
                          mu,
                          nnlsInd = TRUE,
                          tol = 0.001,
                          objCriteria = FALSE,
                          eta = 0.5,
                          dataSplit = 1,
                          projs = 0,
                          Avg = FALSE,
                          wUpdate = "glmnet",
                          standardize = FALSE,
                          xStandardize = FALSE,
                          sampSzWeight = 4,
                          weights = NULL,
                          sigK = NULL,
                          sigStack = NULL,
                          orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(nnlsInd){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # No nonegative least squares for glmnet paramter
        lowLim <- -Inf
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            w <- sapply(wStart[-1], function(x) max(0, x) )
        }else{

            w <- wStart[-1]

        }

    }



    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N
        weights <- Diagonal( x = 1 / sigStack[data$Study] )

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }

        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k] , p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        if(xStandardize){
            # standardize the kth study's design matrix by its sd  only on the RHS

            # add one up front (and take one element of the end) to account for intercept's variance (which is 0)
            sd_xVec <- c(1, sqrt(apply(X[indx_k_SSL,], 2, var) * (nVec[k]-1) / nVec[k] )[-1] )
            dMat <- diag( 1 / sd_xVec )

            XX_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,]) %*% X[indx_k_SSL,] %*% dMat
            Xy_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
            #print(sd_xVec)
            rm(dMat, sd_xVec)

        }else{
            # no standardization of RHS (Study specific portion of loss)

            XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
            Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
        }



    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] ) # X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    X <- X[Stackindx,]
    Xnrow <- N <- nrow(X) # number of rows in all studies
    y_sum <- sum(y[Stackindx])
    #muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOEC(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = TRUE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] ) # * u4[k]   --- above
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC(data = data,
                            beta = betaTemp,
                            w = c(w0, w),
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = TRUE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                mod <- glmnet(y = y,
                              x = as.matrix(X %*% beta),
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights),
                              thresh = 1e-10)

                # w coefficients
                wT <-  as.vector(mod$beta)
                wT0 <- as.vector(mod$a0)

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, wT),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                        w <- rep(1 / K, K) # should be redundant

                        # w0 update
                        wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                        obj_k <- objOEC(data = data,
                                        beta = beta,
                                        w = c(wT0, w),
                                        mu = mu,
                                        lambdaVec = lambdaVec,
                                        stackInt = TRUE,
                                        eta = eta,
                                        StackIndx = Stackindx,
                                        sampSzWeight = sampSzWeight,
                                        sigK = sigK,
                                        sigStack = sigStack,
                                        weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                        )

                        if(obj_k < min(objVec) ){
                            # if the objective improves save it as current best
                            w0 <- wT0

                            cur <- cur + 1
                            objImp <- c(objImp, obj_k)

                        }

                        objVec <- c(objVec, obj_k) # add objective even if no improvement

                    }

        }

        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}

#############
# Specialist
#############
# was ridgeAltSpec until September 14,2020 but updates didnt account for changing w[k] at each iteration
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
ridgeAltSpecFAST <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         Stackindx,
                         mu,
                         nnlsInd = TRUE,
                         low = 0,
                         up = Inf,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL,
                         stackTol = 1e-9){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx)
        weights <- Diagonal( x = 1 / sigStack[data$Study] )
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))

                                           # /
    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
      lambdaList <- vector("list", length = K)

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    # muMat <- diag( rep( mu, K )   )


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )

    for(k in 1:K){

        indx <- which(Study == k)
        lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list <- t(X[indx,]) %*% X[indx,]
        Xy_k_list <- t(X[indx,])  %*% y[indx]

        inv_list[[k]] <- as.matrix( solve(
            eta * u1 * w[k]^2 * XX +
                (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )
                                    ) )


        mat[,k] <- as.vector(inv_list[[k]] %*% (
                                                    ( eta * u1 * w[k] *  Xy  ) +
                                                    (1 - eta) * u3[k] * Xy_k_list
                                                ) )

        mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )

        rm(XX_k_list, Xy_k_list)

    }

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_S(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){

            betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) +
                                    mat[,k] - w0 * mat2[,k]

            obj_k <- objOEC_S(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    # w coefficients
                    wT <-  as.vector(mod$beta)
                    wT0 <- as.vector(mod$a0)

                }else{
                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                    mod <- stackLim(y = y,
                                          x = as.matrix(X %*% beta),
                                          low = low,
                                          # up = up,  doesnt currently work for upper lims
                                          lambda = mu,
                                          w0 = c(w0, w), # warm start
                                          tol = stackTol,  # better results from 1e-12
                                          u1 = u1,
                                          u2 = u2,
                                          weights = weights)


                    wT <- mod[-1]
                    wT0 <- mod[1]
                }


                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}



#############
# Specialist
#############
# stopped using Setpember 14, 2020 because doesnt account for w[k] changing at each iteration
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
# ridge alt specialist but zero out specialist country in w (stacking weights)
ridgeAltSpec0FAST <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         Stackindx,
                         mu,
                         nnlsInd = TRUE,
                         low = 0,
                         up = Inf,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL,
                         stackTol = 1e-9){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )
    nTest <- length(Stackindx )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

        if(length( indx ) == nTest){

            if(all(indx == Stackindx)){
                # test country
                testIndx <- j
            }

        }


    }

    b0 <- rep(0, nrow(betaStart)) # to 0 out
    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    # if(nnlsInd & low <= 0){
    #     # nonegative least squares for glmnet paramter
    #     lowLim <- 0
    # }else{
    #     # No nonegative least squares for glmnet paramter
    #     lowLim <- -Inf
    # }

    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level

    # zero out
    betaStart[,testIndx] <- 0
    wStart[1 + testIndx] <- 0

    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }



    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) ## LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx) #  # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******

    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    # muMat <- diag( rep( mu, K )   )


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )

    for(k in 1:K){

        indx <- which(Study == k)
        lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list <- t(X[indx,]) %*% X[indx,]
        Xy_k_list <- t(X[indx,])  %*% y[indx]

        inv_list[[k]] <- as.matrix( solve(
            eta * u1 * w[k]^2 * XX +
                (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )
        ) )


        mat[,k] <- as.vector(inv_list[[k]] %*% (
            ( eta * u1 * w[k] *  Xy  ) +
                (1 - eta) * u3[k] * Xy_k_list
        ) )

        mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )

        rm(XX_k_list, Xy_k_list)

    }

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_S(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){

            if(k == testIndx){
                betaTemp[,k] <- b0 # 0 out the coefficient corresponding to the test study
            }else{
                betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) +
                    mat[,k] - w0 * mat2[,k]
            }

            obj_k <- objOEC_S(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    # w coefficients
                    wT <-  as.vector(mod$beta)
                    wT0 <- as.vector(mod$a0)

                    wT[testIndx] <- 0 # zero out the coefficent associated with the test set

                }else{
                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                    mod <- stackLim(y = y,
                                    x = as.matrix(X %*% beta),
                                    low = low,
                                    # up = up,  doesnt currently work for upper lims
                                    lambda = mu,
                                    w0 = c(w0, w), # warm start
                                    tol = stackTol,  # better results from 1e-12
                                    u1 = u1,
                                    u2 = u2,
                                    weights = weights)


                    wT <- mod[-1]
                    wT0 <- mod[1]

                    wT[testIndx] <- 0 # zero out the coefficent associated with the test set
                }


                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}

#############
# Specialist
#############
# Started using Setpember 14, 2020 because old one (RidgeAltSpecFast) didnt account for w[k] changing at each update
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
ridgeAltSpec0 <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         Stackindx,
                         mu,
                         nnlsInd = TRUE,
                         low = 0,
                         up = Inf,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         xStandardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL,
                         stackTol = 1e-9,
                         pcaMat = NULL,
                         orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset
    # pcaMat is a matrix of loading vectors

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) #
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )
    nTest <- length(Stackindx )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

        if(length( indx ) == nTest){

            if(all(indx == Stackindx)){
                # test country
                testIndx <- j
            }

        }


    }

    wStart[testIndx + 1] <- 0 # zero out coefficient corresponding to test training set (+1 because of intercept)
    betaStart[,testIndx] <- 0 # zero out coefficients corresponding to test training study
    b0 <- rep(0, nrow(betaStart)) # to 0 out
    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)


    if(low <= 0 || up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf && low != -Inf){

            w <- wStart[-1]

        }

    }

    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] /
        #     sum(diag(weights[Stackindx, Stackindx]))

    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx) # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning


    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }


        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k] , p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        if(xStandardize){
            # standardize the kth study's design matrix by its sd  only on the RHS

            # add one up front (and take one element of the end) to account for intercept's variance (which is 0)
            sd_xVec <- c(1, sqrt(apply(X[indx_k_SSL,], 2, var) * (nVec[k]-1) / nVec[k] )[-1] )
            dMat <- diag( 1 / sd_xVec )

            XX_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,]) %*% X[indx_k_SSL,] %*% dMat
            Xy_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
            #print(sd_xVec)
            rm(dMat, sd_xVec)

        }else{
            # no standardization of RHS (Study specific portion of loss)

            XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
            Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
        }

    }


    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    # muMat <- diag( rep( mu, K )   )


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_S(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){

            if(k == testIndx){
                betaTemp[,k] <- b0 # 0 out the coefficient corresponding to the test study
            }else{
                inv <- solve(
                    eta * u1 * w[k]^2 * XX +
                        (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  ) # * u4[k]
                )

                betaTemp[,k] <- inv %*% (
                    eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                        (1 - eta) * u3[k] * Xy_k_list[[k]]
                )
            }

            obj_k <- objOEC_S(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    # w coefficients
                    wT <-  as.vector(mod$beta)
                    wT0 <- as.vector(mod$a0)

                    wT[testIndx] <- 0 # zero out stacking

                }else{
                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                    mod <- stackLim(y = y,
                                    x = as.matrix(X %*% beta),
                                    low = low,
                                    # up = up,  doesnt currently work for upper lims
                                    lambda = mu,
                                    w0 = c(w0, w), # warm start
                                    tol = stackTol,  # better results from 1e-12
                                    u1 = u1,
                                    u2 = u2,
                                    weights = weights)


                    wT <- mod[-1]
                    wT0 <- mod[1]

                    wT[testIndx] <- 0 # zero out stacking
                }


                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}

#############
# Specialist
#############

ridgeAltSpec <- function(data,
                            betaStart,
                            wStart,
                            lambdaVec,
                            Stackindx,
                            mu,
                            nnlsInd = TRUE,
                            low = 0,
                            up = Inf,
                            tol = 0.001,
                            objCriteria = FALSE,
                            eta = 0.5,
                            dataSplit = 1,
                            projs = 0,
                            Avg = FALSE,
                            wUpdate = "glmnet",
                            standardize = FALSE,
                            xStandardize = FALSE,
                            sampSzWeight = 4,
                            weights = NULL,
                            sigK = NULL,
                            sigStack = NULL,
                            stackTol = 1e-9,
                            pcaMat = NULL,
                            orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset
    # pcaMat is a matrix of loading vectors

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx)
        weights <- Diagonal( x = 1 / sigStack[data$Study] )
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }

        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k] , p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        if(xStandardize){
            # standardize the kth study's design matrix by its sd  only on the RHS

            # add one up front (and take one element of the end) to account for intercept's variance (which is 0)
            sd_xVec <- c(1, sqrt(apply(X[indx_k_SSL,], 2, var) * (nVec[k]-1) / nVec[k] )[-1] )
            dMat <- diag( 1 / sd_xVec )

            XX_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,]) %*% X[indx_k_SSL,] %*% dMat
            Xy_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
            #print(sd_xVec)
            rm(dMat, sd_xVec)

        }else{
            # no standardization of RHS (Study specific portion of loss)

            XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
            Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
        }

    }


    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_S(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  )
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC_S(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    # w coefficients
                    wT <-  as.vector(mod$beta)
                    wT0 <- as.vector(mod$a0)

                }else{
                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                    mod <- stackLim(y = y,
                                    x = as.matrix(X %*% beta),
                                    low = low,
                                    lambda = mu,
                                    w0 = c(w0, w), # warm start
                                    tol = stackTol,
                                    u1 = u1,
                                    u2 = u2,
                                    weights = weights)


                    wT <- mod[-1]
                    wT0 <- mod[1]
                }


                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}
#############
# Specialist
#############

ridgeAltSpecPCA <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         Stackindx,
                         mu,
                         nnlsInd = TRUE,
                         low = 0,
                         up = Inf,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL,
                         stackTol = 1e-9,
                         pcaMat = NULL,
                         orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset
    # pcaMat is a matrix of loading vectors

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    # removed for speed

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level

    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] )  # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx) # # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******

    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }



        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k], p - 1)   ) )#Diagonal( x = c(0, rep( lambdaVec[k] * u4[k], p - 1)   )  )  # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }


    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    # muMat <- diag( rep( mu, K )   )


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_SPCA(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL, # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                                         pcaMat = pcaMat
                                         )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    while(eps > tol){

        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){

            inv <- solve(
                eta * u1 * ((pcaMat %*% w)[k])^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  ) # * u4[k]
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * (pcaMat %*% w)[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% (pcaMat %*% w)[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC_SPCA(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL, # weights
                              pcaMat = pcaMat
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta %*% pcaMat),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    wT <-  as.vector( coef(mod)[-1] ) #pcaMat %*%

                    # w coefficients
                    wT0 <- as.vector( coef(mod) )[1] # as.vector(mod$a0)

                }else{
                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                    mod <- stackLim(y = y,
                                    x = as.matrix(X %*% beta %*% pcaMat),
                                    low = low,
                                    # up = up,  doesnt currently work for upper lims
                                    lambda = mu,
                                    w0 = c(w0, w), # warm start
                                    tol = stackTol,  # better results from 1e-12
                                    u1 = u1,
                                    u2 = u2,
                                    weights = weights)

                    # w coefficients
                    wT0 <- (mod)[1] # as.vector(mod$a0)
                    wT <-  mod[-1] #as.vector(mod$beta) pcaMat %*%


                }


                obj_k <- objOEC_SPCA(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL, # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                                  pcaMat = pcaMat
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% pcaMat %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_SPCA(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL, # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                                  pcaMat = pcaMat
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}


#############
# Specialist
#############
# make the test country training set be constrained to have stacking weight <= 1 / K
# Started using Setpember 14, 2020 because old one (RidgeAltSpecFast) didnt account for w[k] changing at each update
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
ridgeAltSpecWindow <- function(data,
                         betaStart,
                         wStart,
                         lambdaVec,
                         Stackindx,
                         mu,
                         nnlsInd = TRUE,
                         low = 0,
                         up = 1 / K,
                         tol = 0.001,
                         objCriteria = FALSE,
                         eta = 0.5,
                         dataSplit = 1,
                         projs = 0,
                         Avg = FALSE,
                         wUpdate = "glmnet",
                         standardize = FALSE,
                         sampSzWeight = 4,
                         weights = NULL,
                         sigK = NULL,
                         sigStack = NULL,
                         stackTol = 1e-9,
                         pcaMat = NULL,
                         orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses
    # Stackindx is a vector of row indices corresponding to the target dataset
    # pcaMat is a matrix of loading vectors

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) #
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )
    nTest <- length(Stackindx )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

        if(length( indx ) == nTest){

            if(all(indx == Stackindx)){
                # test country
                testIndx <- j
            }

        }


    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)


    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(Stackindx) # # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(Stackindx) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

    }


    itr <- cur <- 1

    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }



        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k], p - 1)   ) )#Diagonal( x = c(0, rep( lambdaVec[k] * u4[k], p - 1)   )  )  # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }


    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL

    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL # do not need these weights anymore
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize this
    }
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )

    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking

    objVec <- objImp <- obj0 <- objOEC_S(data = data,
                                         beta = betaStart,
                                         w = wStart,
                                         mu = mu,
                                         lambdaVec = lambdaVec,
                                         stackInt = TRUE,
                                         eta = eta,
                                         StackIndx = Stackindx,
                                         sampSzWeight = sampSzWeight,
                                         sigK = sigK,
                                         sigStack = sigStack,
                                         weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){


            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  ) # * u4[k]
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC_S(data = data,
                              beta = betaTemp,
                              w = c(w0, w),
                              mu = mu,
                              lambdaVec = lambdaVec,
                              stackInt = TRUE,
                              eta = eta,
                              StackIndx = Stackindx,
                              sampSzWeight = sampSzWeight,
                              sigK = sigK,
                              sigStack = sigStack,
                              weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                    # use custom function if lower limits > 0
                    # note upper limits cannot be used here
                # only constrains upper limits for test index
                    mod <- stackLimUp(y = y,
                                    x = as.matrix(X %*% beta),
                                    low = low,
                                    up = up,  #doesnt currently work for upper lims
                                    lambda = mu,
                                    w0 = c(w0, w), # warm start
                                    tol = stackTol,  # better results from 1e-12
                                    u1 = u1,
                                    u2 = u2,
                                    weights = weights,
                                    testIndx = testIndx)


                    wT <- mod[-1]
                    wT0 <- mod[1]
                #}


                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, wT),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC_S(data = data,
                                  beta = beta,
                                  w = c(wT0, w),
                                  mu = mu,
                                  lambdaVec = lambdaVec,
                                  stackInt = TRUE,
                                  eta = eta,
                                  StackIndx = Stackindx,
                                  sampSzWeight = sampSzWeight,
                                  sigK = sigK,
                                  sigStack = sigStack,
                                  weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}
#############
# NewWR
#############
# *** same as ridgeAltFix but no intercept

# fix stacking weights
ridgeAltFixNoIntercept <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 0,
                     Avg = FALSE,
                     wUpdate = "glmnet",
                     standardize = FALSE,
                     sampSzWeight = 4,
                     weights = NULL,
                     sigK = NULL,
                     sigStack = NULL,
                     up = Inf,
                     low = -Inf,
                     simplex = FALSE,
                     stackTol = 1e-6){

    # eta is the parameter that determines convex combination of the losses
    # stackTol is the tolerance for projected gradient descent for stacking

    # solver for w update
    library(CVXR)
    if( is.element("OSQP", installed_solvers()) ){
        solverNm <- "OSQP"
    }else{
        solverNm <- "OSQP"
    }

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if( is.na(nnlsInd)){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # specify value
        lowLim <- nnlsInd
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart

    if(Avg){
        len <- length(wStart)
        w <- wStart <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart

        }

    }

    if(length(wStart) > K)     wStart <- wStart[-1] # chop off intercept if there still is one

    w_Star <- w # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }


        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] # X^TX calculate so do not have to at each iteration
    Xy <- t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    X <- X[Stackindx,]
    Xnrow <- N <- nrow(X) # number of rows in all studies
    y_sum <- sum(y[Stackindx])
    muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOEC(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = FALSE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################
    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC(data = data,
                            beta = betaTemp,
                            w = w,
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = FALSE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                mod <-       stackLim(y = y,
                                     x = as.matrix(X %*% beta),
                                     low = low,
                                     # up = up,  doesnt currently work for upper lims
                                     lambda = mu,
                                     intercept = FALSE,
                                     w0 = w, # warm start
                                     tol = stackTol,  # better results from 1e-12
                                     u1 = u1,
                                     u2 = u2,
                                     weights = weights)


                wT <- mod

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = wT,
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = FALSE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c( wT )
                    w <- wT
                    # w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = w,
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = FALSE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }

        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]


    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}

#############
# NewWR
#############
# fix stacking weights
ridgeAltFix <- function(data,
                        betaStart,
                        wStart,
                        lambdaVec,
                        mu,
                        nnlsInd = TRUE,
                        tol = 0.001,
                        objCriteria = FALSE,
                        eta = 0.5,
                        dataSplit = 1,
                        projs = 0,
                        Avg = FALSE,
                        wUpdate = "glmnet",
                        standardize = FALSE,
                        sampSzWeight = 4,
                        weights = NULL,
                        sigK = NULL,
                        sigStack = NULL,
                        up = Inf,
                        low = -Inf,
                        simplex = FALSE,
                        stackTol = 1e-6,
                        orderRandom = TRUE){

    # eta is the parameter that determines convex combination of the losses
    # stackTol is the tolerance for projected gradient descent for stacking

    # solver for w update
    library(CVXR)
    if( is.element("OSQP", installed_solvers()) ){
        solverNm <- "OSQP"
    }else{
        solverNm <- "OSQP"
    }

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if( is.na(nnlsInd)){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # specify value
        lowLim <- nnlsInd
    }

    if(low <= 0 | up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }



        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k] , p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }

    #
    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop
    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize
    }

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOEC(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = TRUE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )


    if(dataSplit == 1)    Stackindx <- NULL

    #############################
    # Iteratively update
    #############################

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        if(!orderRandom)      stud <- 1:K   # do not randomize order

        for(k in stud){

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  )
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )



            obj_k <- objOEC(data = data,
                            beta = betaTemp,
                            w = c(w0, w),
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = TRUE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                if(glmnetInd){
                    # use glmnet for optimization
                    mod <- glmnet(y = y,
                                  x = as.matrix(X %*% beta),
                                  alpha = 0,
                                  lambda = mu,
                                  lower.limits = low,
                                  upper.limits = up,
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = weights,
                                  thresh = 1e-10)

                    wT <-  as.vector( coef(mod)[-1] ) #pcaMat %*%

                    # w coefficients
                    wT0 <- as.vector( coef(mod) )[1] # as.vector(mod$a0)
                    #wT <-  wT[-1] #as.vector(mod$beta)

                }else{

                mod <-       stackLim(y = y,
                                      x = as.matrix(X %*% beta),
                                      low = low,
                                      lambda = mu,
                                      w0 = c(w0, w), # warm start
                                      tol = stackTol,  # better results from 1e-12
                                      u1 = u1,
                                      u2 = u2,
                                      weights = weights)

                    wT <- mod[-1]
                    wT0 <- mod[1]
                }

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, wT),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, w),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}


#############
# NewWR
#############

ridgeAltFixFAST <- function(data,
                        betaStart,
                        wStart,
                        lambdaVec,
                        mu,
                        nnlsInd = TRUE,
                        tol = 0.001,
                        objCriteria = FALSE,
                        eta = 0.5,
                        dataSplit = 1,
                        projs = 0,
                        Avg = FALSE,
                        wUpdate = "glmnet",
                        standardize = FALSE,
                        sampSzWeight = 4,
                        weights = NULL,
                        sigK = NULL,
                        sigStack = NULL,
                        up = Inf,
                        low = -Inf,
                        simplex = FALSE,
                        stackTol = 1e-6){

    # eta is the parameter that determines convex combination of the losses
    # stackTol is the tolerance for projected gradient descent for stacking

    # solver for w update
    library(CVXR)
    if( is.element("OSQP", installed_solvers()) ){
        solverNm <- "OSQP"
    }else{
        solverNm <- "OSQP"
    }

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if( is.na(nnlsInd)){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # specify value
        lowLim <- nnlsInd
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(low != -Inf){
            # truncate from below

            w <- wStart[-1] <- sapply(wStart[-1], function(x) max(low, x) )
        }

        if(up != Inf){
            # truncate from above

            w <- wStart[-1] <- sapply(wStart[-1], function(x) min(up, x) )
        }

        # if neither, then do not alter
        if(up != Inf & low != -Inf){

            w <- wStart[-1]

        }

    }


    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N #  LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] )
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)

    SSLindx <- c()
    Stackindx <- c()

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] )# X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    y_sum <- sum(y[Stackindx])
    X <- X[Stackindx,]
    y <- y[Stackindx]
    Xnrow <- nrow(X) # number of rows in stacking


    if(sampSzWeight < 6){
        weights <- rep(1, length( Stackindx) ) # NULL
    }else{
        weights <- diag( weights[Stackindx, Stackindx] ) # vectorize
    }

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOEC(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = TRUE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #########################
    # Inverse Matrix
    #########################
    inv_list <- vector("list", length = K)
    mat <- matrix(ncol = K, nrow = p )

    for(k in 1:K){

        indx <- which(Study == k)
        lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list <- t(X[indx,]) %*% X[indx,]
        Xy_k_list <- t(X[indx,])  %*% y[indx]

        inv_list[[k]] <- as.matrix( solve(
            eta * u1 * w[k]^2 * XX +
                (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )
        ) )


        mat[,k] <- as.vector(inv_list[[k]] %*% (
            ( eta * u1 * w[k] * ( Xy - w0 * X_rowSums ) ) +
                (1 - eta) * u3[k] * Xy_k_list
        ) )

        rm(XX_k_list, Xy_k_list)

    }

    if(dataSplit == 1)    Stackindx <- NULL

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight
    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){


            betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + mat[,k]


            obj_k <- objOEC(data = data,
                            beta = betaTemp,
                            w = c(w0, w),
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = TRUE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                mod <-       stackLim(y = y,
                                      x = as.matrix(X %*% beta),
                                      low = low,
                                      lambda = mu,
                                      w0 = c(w0, w), # warm start
                                      tol = stackTol,  # better results from 1e-12
                                      u1 = u1,
                                      u2 = u2,
                                      weights = weights)


                wT <- mod[-1]
                wT0 <- mod[1]

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, wT),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC(data = data,
                                beta = beta,
                                w = c(wT0, w),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }


        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}

#############
# NewWR
#############
# standardize each design matrix based upon kth studies mean and variance

ridgeAltX <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 0,
                     Avg = FALSE,
                     wUpdate = "glmnet",
                     standardize = FALSE,
                     sampSzWeight = 4,
                     weights = NULL,
                     sigK = NULL,
                     sigStack = NULL){

    # eta is the parameter that determines convex combination of the losses

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(nnlsInd){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # No nonegative least squares for glmnet paramter
        lowLim <- -Inf
    }

    # removed for speed
    # if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            w <- sapply(wStart[-1], function(x) max(0, x) )
        }else{

            w <- wStart[-1]

        }

    }



    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    X_list <- vector("list", length = K)
    meanMat <- matrix(ncol = K, nrow = ncol(data) - 2) # matrix of covariate means for each study
    sigmaMat <- matrix(ncol = K, nrow = ncol(data) - 2) # matrix of covariate sds for each study

    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }


        nFull <- length(indx_k_SSL) # sample size of merged

        # scale Covaraites
        means <- colMeans(as.matrix(X[indx_k_SSL,-1]))
        sds <- sqrt( apply(as.matrix(X[indx_k_SSL,-1]), 2, var) *  (nFull - 1) / nFull ) # use mle formula to match with GLMNET
        sds <- ifelse(sds == 0, 1, sds) # replace sd = 0 to 1 so it divides by 1
        meanMat[,k] <- means
        sigmaMat[,k] <- sds
        #
        X_list[[k]] <- matrix(nrow = nrow(X), ncol = ncol(X) )
        X_list[[k]][,1] <- 1

        for(column in 2:(ncol(X) )){
            # starts from 2 to skip intercet
            # center scale entire design matrix
            X_list[[k]][, column ] <- (X[, column] - means[column - 1]) / sds[column - 1]
        }

        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        XX_k_list[[k]] <- t(X_list[[k]][indx_k_SSL,]) %*% X_list[[k]][indx_k_SSL,]
        Xy_k_list[[k]] <- t(X_list[[k]][indx_k_SSL,])  %*% y[indx_k_SSL]
    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    # XX <- t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] # X^TX calculate so do not have to at each iteration
    # Xy <- t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] # X^T y
    # X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    X <- X[Stackindx,]
    Xnrow <- N <- nrow(X) # number of rows in all studies
    y_sum <- sum(y[Stackindx])
    muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    # if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOECX(data = data,
                                       beta = betaStart,
                                       w = wStart,
                                       mu = mu,
                                       lambdaVec = lambdaVec,
                                       stackInt = TRUE,
                                       eta = eta,
                                       StackIndx = Stackindx,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
    )

    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){
            XXw <- 0

            for(l in seq(1,K)[-k]){
                XXw <- XXw + w[l] * X_list[[l]][Stackindx,] %*% beta[,l]
            }

            XXw <- t(X_list[[k]][Stackindx,]) %*% weights[Stackindx, Stackindx] %*% XXw

            XX <- t(X_list[[k]][Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X_list[[k]][Stackindx,]

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( t(X_list[[k]][Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] -
                    w0 * rowSums( t(X_list[[k]][Stackindx,]) %*% weights[Stackindx, Stackindx] ) - XXw ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOECX(data = data,
                            beta = betaTemp,
                            w = c(w0, w),
                            mu = mu,
                            lambdaVec = lambdaVec,
                            stackInt = TRUE,
                            eta = eta,
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            sigK = sigK,
                            sigStack = sigStack,
                            weights = NULL # weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){
                predsMat <- predsX(meanMat = meanMat,
                                   sigmaMat = sigmaMat,
                                   data = data[,-c(1,2)],
                                   beta = beta)

                mod <- glmnet(y = y,
                              x = as.matrix(predsMat),
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights),
                              thresh = 1e-10)

                # w coefficients
                wT <-  as.vector(mod$beta)
                wT0 <- as.vector(mod$a0)

                obj_k <- objOECX(data = data,
                                beta = beta,
                                w = c(wT0, wT),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20


                obj_k <- objOECX(data = data,
                                beta = beta,
                                w = c(wT0, w),
                                mu = mu,
                                lambdaVec = lambdaVec,
                                stackInt = TRUE,
                                eta = eta,
                                StackIndx = Stackindx,
                                sampSzWeight = sampSzWeight,
                                sigK = sigK,
                                sigStack = sigStack,
                                weights = NULL # weights  -- objOEC doesnt take weights anymore, and stacking variance adjustment is done another way
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }

        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]

    }

    return(list(beta = beta_Star,
                w = w_Star,
                obj = objVec,
                itrs = itr,
                objImp = objImp,
                sigmaMat = sigmaMat,
                meanMat = meanMat))

}

##########################
# Zero Out OEC
##########################
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using 10/4/20
ridgeAlt0 <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 0,
                     Avg = FALSE,
                     wUpdate = "glmnet",
                     standardize = FALSE,
                     xStandardize = FALSE,
                     sampSzWeight = 4,
                     weights = NULL,
                     sigK = NULL,
                     sigStack = NULL){

    # eta is the parameter that determines convex combination of the losses

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    if(nnlsInd){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # No nonegative least squares for glmnet paramter
        lowLim <- -Inf
    }

    # removed for speed

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(Avg){
        len <- length(wStart[-1])
        w <- wStart[-1] <- rep(1 / K, K)
    }else{

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            w <- sapply(wStart[-1], function(x) max(0, x) )
        }else{

            w <- wStart[-1]

        }

    }



    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    nVec <- as.vector( table(data$Study) ) # number of rows in each study
    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)
    N <- nrow(data)

    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) # Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK # aacount for tuning

        if(dataSplit != 1){
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
            weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        }
    }


    itr <- cur <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }

        lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k] , p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty

        if(xStandardize){
            # standardize the kth study's design matrix by its sd  only on the RHS

            # add one up front (and take one element of the end) to account for intercept's variance (which is 0)
            sd_xVec <- c(1, sqrt(apply(X[indx_k_SSL,], 2, var) * (nVec[k]-1) / nVec[k] )[-1] )
            dMat <- diag( 1 / sd_xVec )

            XX_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,]) %*% X[indx_k_SSL,] %*% dMat
            Xy_k_list[[k]] <- dMat %*% t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
            #print(sd_xVec)
            rm(dMat, sd_xVec)

        }else{
            # no standardization of RHS (Study specific portion of loss)

            XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
            Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
        }



    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- as.matrix( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% X[Stackindx,] ) # X^TX calculate so do not have to at each iteration
    Xy <- as.vector( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] %*% y[Stackindx] ) # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) %*% weights[Stackindx, Stackindx] )
    X <- X[Stackindx,]
    Xnrow <- N <- nrow(X) # number of rows in all studies
    y_sum <- sum(y[Stackindx])
    #muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objImp <- obj0 <- objOEC0(data = data,
                                        beta = betaStart,
                                        w = wStart,
                                        mu = mu,
                                        lambdaVec = lambdaVec,
                                        stackInt = TRUE,
                                        eta = eta,
                                        StackIndx = Stackindx,
                                        sampSzWeight = sampSzWeight,
                                        weights = weights
                                        )


    #############################
    # Iteratively update
    #############################

    # check to see how objective is weighted
    # sampSzWeight

    # nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)

    # if scale LHS of objective by 1 / N to account for variable sample sizes
    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1

        betaTemp <- beta # beta used as temp

        # beta update
        stud <- sample.int(K, K, replace = FALSE)
        for(k in stud){

            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] ) #
            )

            betaTemp[,k] <- inv %*% (
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
            )

            obj_k <- objOEC0(data = data,
                             beta = betaTemp,
                             w = c(w0, w),
                             mu = mu,
                             lambdaVec = lambdaVec,
                             stackInt = TRUE,
                             eta = eta,
                             StackIndx = Stackindx,
                             sampSzWeight = sampSzWeight,
                             weights = weights
            )

            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

                cur <- cur + 1
                objImp <- c(objImp, obj_k)

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }

            objVec <- c(objVec, obj_k) # add objective even if no improvement


            # w update
            if(!Avg){

                mod <- glmnet(y = y,
                              x = as.matrix(X %*% beta),
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights),
                              thresh = 1e-10)

                # w coefficients
                wT <-  as.vector(mod$beta)
                wT0 <- as.vector(mod$a0)

                obj_k <- objOEC0(data = data,
                                 beta = beta,
                                 w = c(wT0, wT),
                                 mu = mu,
                                 lambdaVec = lambdaVec,
                                 stackInt = TRUE,
                                 eta = eta,
                                 StackIndx = Stackindx,
                                 sampSzWeight = sampSzWeight,
                                 weights = weights
                                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w_Star <- c(wT0, wT)
                    w <- wT
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }else{
                ####################
                # if average weights
                ####################

                w <- rep(1 / K, K) # should be redundant

                # w0 update
                wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20

                obj_k <- objOEC0(data = data,
                                 beta = beta,
                                 w = c(wT0, w),
                                 mu = mu,
                                 lambdaVec = lambdaVec,
                                 stackInt = TRUE,
                                 eta = eta,
                                 StackIndx = Stackindx,
                                 sampSzWeight = sampSzWeight,
                                 weights = weights
                )

                if(obj_k < min(objVec) ){
                    # if the objective improves save it as current best
                    w0 <- wT0

                    cur <- cur + 1
                    objImp <- c(objImp, obj_k)

                }

                objVec <- c(objVec, obj_k) # add objective even if no improvement

            }

        }

        eps <- obj0 - objImp[cur]

        obj0 <- objImp[cur]


    }

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))

}



ridgeAltAvg <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 0,
                     wUpdate = "glmnet"){

    # eta is the parameter that determines convex combination of the losses

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    if(nnlsInd){
        # if there is a nonnegativity constraint on the stacking weights
        # do not apply to intercept

        len <- length(wStart[-1])
        w <- rep(1 / len, len)

    }else{

        len <- length(wStart[-1])
        w <- rep(1 / len, len)

    }

    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    wOld <- wStart # for update difference calcs
    betaOld <- betaStart

    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)

    itr <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- list(length = K)
    XX_k_list <- list(length = K)
    Xy_k_list <- list(length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }



        lambdaList[[k]] <- diag( rep( lambdaVec[k], p )   )
        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- t(X[Stackindx,]) %*% X[Stackindx,] # X^TX calculate so do not have to at each iteration
    Xy <- t(X[Stackindx,]) %*% y[Stackindx] # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) )
    X <- X[Stackindx,]
    Xnrow <- nrow(X)
    y_sum <- sum(y[Stackindx])
    muMat <- diag( rep( mu, K )   )

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objOEC(data = data,
                     beta = betaStart,
                     w = wStart,
                     mu = mu,
                     lambdaVec = lambdaVec,
                     stackInt = TRUE,
                     eta = eta,
                     StackIndx = Stackindx)

    #############################
    # Iteratively update
    #############################

    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1
        # beta update
        for(k in 1:K){

            inv <- solve(
                eta * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] + lambdaList[[k]] )
            )

            beta[,k] <- inv %*% (
                eta * w[k] * ( Xy - w0 * X_rowSums  - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta)  *  Xy_k_list[[k]]
            )
        }

        if(wUpdate == "glmnet"){
            # if use glmnet to update parameters
            mod <- glmnet(y = y,
                          x = as.matrix(X %*% beta),
                          alpha = 0,
                          lambda = mu,
                          lower.limits = lowLim,
                          standardize = TRUE,
                          intercept = TRUE)

            # coefficients
            w <-  as.vector(mod$beta)
            w0 <- as.vector(mod$a0)

        }else{

            # w0 update
            w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow

            # w update
            bInv <- solve(    eta * ( t(beta) %*% XX %*% beta + muMat )    )
            w <- bInv %*% (
                eta * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) ) )

            )

        }

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            if(projs > 0){

                # if project onto unit simplex

                for(pr in 1:projs){
                    # Projection onto Non Negative Orthant
                    w <- sapply(w, function(x) max(0, x) )
                    # Projection onto Sum of 1
                    w <- w - ((sum(w) - 1) / K)

                    # Projection onto Non Negative Orthant
                    w <- sapply(w, function(x) max(0, x) )
                }

            }


        }

        obj_r <- objOEC(data = data,
                        beta = beta,
                        w = c(w0, w),
                        mu = mu,
                        lambdaVec = lambdaVec,
                        stackInt = TRUE,
                        eta = eta,
                        StackIndx = Stackindx)



        if(obj_r < objVec[itr - 1]){
            # if the objective improves save it as current best
            beta_Star <- beta
            w_Star <- c(w0, w)
        }else{
            # if it is not better, then do not update the best

            if(objCriteria == TRUE){
                # if obj criteria is true, use the best iteration as this iterations value
                w <- w_Star[-1]
                w0 <- w_Star[1]
                beta <- beta_Star
            }
        }

        # calculate change from last iteration
        eps <- sum( (betaOld - beta)^2  ) + sum( (wOld - c(w0, w))^2  )

        objVec <- c(objVec, obj_r )
        # update previous iteration values
        betaOld <- beta
        wOld <- c(w0, w)
    }

    #print(paste("Final Iterations ", itr))
    # wStar <- c(w0, w) # concatenate stacking weights and intercept

    # objVal <- objOEC(data = data,
    #                  beta = beta_Star,
    #                  w = w_Star,
    #                  mu = mu,
    #                  lambdaVec = lambdaVec,
    #                  stackInt = TRUE,
    #                  eta = eta)

    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))

}


objLinear <- function(X,
                      Y,
                      beta,
                      weights = 1){

    # returns RSS from OLS
    # as.vector( diag( weights) ) is there in case a diagonal matrix of weights is passed
    # this can handle either just the scalar 1 (default) or a vector of weights or a diagonal matrix of weights

    return( sum(  as.vector( diag( weights) ) * (Y - X %*% beta)^2  ) )
}



#  objective of Ridge OEC
objOEC <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL,
                   sigK = NULL,
                   sigStack = NULL){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS
    # sigK is a K x 1 vector of study specific residual variances from study specific regressions
    # sigStack is a K x 1 vector of study specific residual variances from the stacking regression

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )

    if(stackInt){

        s <- -1  # remove first index so do not penalize intercept in stacking term
    }else{
        s <- 1:length(w) # do not remove first index because no intercept
    }

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N ## LHS Objective weight   - N x 1 vector
         # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

    }


    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2] - w[1] - cbind(1, data[,-c(1,2)]) %*% as.matrix( beta ) %*% w[-1] )^2 )

    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2]  -  cbind(1, data[,-c(1,2)]) %*% beta %*% w)^2     )


    }

    for(j in 1:K){

        indx <- which(Study == j) # rows from this study

        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% as.matrix( beta[-1, j] ) )^2     )


    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[s]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}

#  objective of Ridge OEC
objOECPCA <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL,
                   sigK = NULL,
                   sigStack = NULL,
                   pcaMat){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS
    # sigK is a K x 1 vector of study specific residual variances from study specific regressions
    # sigStack is a K x 1 vector of study specific residual variances from the stacking regression

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )

    if(stackInt){

        s <- -1  # remove first index so do not penalize intercept in stacking term
    }else{
        s <- 1:length(w) # do not remove first index because no intercept
    }

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # # LHS Objective weight   - N x 1 vector
    # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

    }


    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2] - w[1] - cbind(1, data[,-c(1,2)]) %*% beta %*% pcaMat %*% w[-1] )^2 )
      }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2]  -  cbind(1, data[,-c(1,2)]) %*% beta %*% pcaMat %*% w)^2     )

    }

    for(j in 1:K){

        indx <- which(Study == j) # rows from this study

        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% beta[-1, j] )^2     )


    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[s]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}



#  objective of Ridge OEC
objOEC_S <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL,
                   sigK = NULL,
                   sigStack = NULL){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS
    # sigK is a K x 1 vector of study specific residual variances from study specific regressions
    # sigStack is a K x 1 vector of study specific residual variances from the stacking regression

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )

    if(stackInt){

        s <- -1  # remove first index so do not penalize intercept in stacking term
    }else{
        s <- 1:length(w) # do not remove first index because no intercept
    }

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(StackIndx) #  LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight   - N x 1 vector
        # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

    }


    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2] - w[1] -
                                         cbind(1, data[StackIndx, -c(1,2)]) %*% as.matrix( beta ) %*% w[-1] )^2 )
    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2]  -  cbind(1, data[StackIndx,-c(1,2)]) %*% as.matrix( beta ) %*% w)^2     )


    }

    for(j in 1:K){

        indx <- which(Study == j) # rows from this study

        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% as.matrix( beta[-1, j] ) )^2     )



    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[s]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}



#  objective of Ridge OEC
objOEC_SPCA <- function(data,
                     beta,
                     w,
                     mu = NULL,
                     lambdaVec = NULL,
                     stackInt = TRUE,
                     eta = 0.5,
                     StackIndx = NULL,
                     sampSzWeight = 4,
                     weights = NULL,
                     sigK = NULL,
                     sigStack = NULL,
                     pcaMat){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS
    # sigK is a K x 1 vector of study specific residual variances from study specific regressions
    # sigStack is a K x 1 vector of study specific residual variances from the stacking regression

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )

    if(stackInt){

        s <- -1  # remove first index so do not penalize intercept in stacking term
    }else{
        s <- 1:length(w) # do not remove first index because no intercept
    }

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / length(StackIndx) # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / length(StackIndx)  # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(StackIndx)   # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

    }


    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2] - w[1] -
                                         cbind(1, data[StackIndx, -c(1,2)]) %*% beta %*% pcaMat %*% w[-1] )^2 )
    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2]  -  cbind(1, data[StackIndx,-c(1,2)]) %*% beta %*% pcaMat %*% w)^2     )


    }

    for(j in 1:K){

        indx <- which(Study == j) # rows from this study

        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% beta[-1, j] )^2     )


    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[s]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}

#  objective of Ridge OEC -- with standardized (study specific) design matrices
objOECX <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL,
                   sigK = NULL,
                   sigStack = NULL){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS
    # sigK is a K x 1 vector of study specific residual variances from study specific regressions
    # sigStack is a K x 1 vector of study specific residual variances from the stacking regression

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )
    obj <- preds <- 0

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)


    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK

         }

    X_list <- vector("list", length = K)


    for(j in 1:K){

        ###############
        # pre processing
        ###############


        indx <- which(Study == j) # rows from this study

        nFull <- length(indx) # sample size of merged

        # scale Covaraites
        means <- colMeans(as.matrix(data[indx,-c(1,2)]))
        sds <- sqrt( apply(as.matrix(data[indx,-c(1,2)]), 2, var) *  (nFull - 1) / nFull )# use mle formula to match with GLMNET
        sds <- ifelse(sds == 0, 1, sds)
        #

        X <- matrix( ncol = ncol(data) - 1, nrow = nrow(data) )
        X[,1] <- 1 # for intercept

        for(column in 2:ncol(X)){

            # center scale entire design matrix
            X[, column] <- ( data[, column + 1] - means[column - 1] )  / sds[column - 1]
        }



        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] - beta[1, j] - X[indx,-1] %*% beta[-1, j] )^2    )

        #############
        ### stacking
        ############
        # if intercept in stacking regression add column of ones to beta
        if(stackInt){

            preds <- preds + X %*% beta[,j] * w[j + 1]

        }else{

            # stacking function -- divide by N to make objective consistent with glmnet
            preds <- preds + X %*% beta[,j] * w[j]

        }

    }

    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- obj + eta / 2 * sum( u1 * ( data[,2] - w[1] - preds )^2 )

    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- obj + eta / 2 * sum( u1 * ( data[,2] - preds)^2     )

    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[-1]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}

#########################
# Objective of Ridge OEC Zeroed Out
#########################
objOEC0 <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL,
                   sigK = NULL,
                   sigStack = NULL){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
    beta <- as.matrix( beta )
    w <- as.vector( w )

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }else if(sampSzWeight == 5){
        # just variance adjustment on RHS
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / ( rep(sigStack, nVec) * N )  # LHS Objective weight  - N x 1 vector
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
    }


    # # zero out
    predsMat <- cbind(1, data[,-c(1,2)]) %*% beta

    for(studyNum in 1:K){
        sIndex <- which(Study == studyNum)
        predsMat[sIndex, studyNum ] <- 0 # zero out
    }

    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- u1 * eta / 2 * sum( ( data[,2] - w[1] - predsMat %*% w[-1] )^2 )

    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- u1 * eta / 2 * sum( ( data[,2]  -  predsMat %*% beta %*% w)^2     )

    }

    for(j in 1:K){

        indx <- which(Study == j) # rows from this study

        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% as.matrix( beta[-1, j] ) )^2     )


    }

    lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
    # add penalties
    obj <- obj + (1 - eta) * sum( beta[-1,]^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w[-1]^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}

# allows for obsWeights and StackIndx (when some observations are used for both)
objOECweighted <- function(data,
                   beta,
                   w,
                   mu = NULL,
                   lambdaVec = NULL,
                   stackInt = TRUE,
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = 4,
                   weights = NULL){

    # returns objective of OEC Linear
    # sampSzWeight determines how obejective weights sample sizes on RHS

    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))

    if( is.null(weights) ){

        weights <- matrix( diag( N ) ) # make weights identity if nothing provided

    }else{
        weights <- N * weights / sum(weights) # standardize weights
    }


    if(is.null(StackIndx)){
        # if null then use all data for both stacking regression and SSL regression
        SIndx <- SSLIndx <- 1:nrow(data)

    }else{
        # if not then use the rows given above for stacking and the rest for the SSL
        SIndx <- StackIndx
        SSLIndx <- setdiff( 1:nrow(data), StackIndx )
    }

    if(is.null(mu)){
        # no stacking penalty
        mu <- 0
    }

    if(is.null(lambdaVec)){
        # no stacking penalty
        lambdaVec <- rep(0, K)
    }

    # different weights for different loss functions
    # different weights for different loss functions
    if(sampSzWeight == 1){

        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / nVec
        u4 <- rep(1, K)
    }else if(sampSzWeight == 2){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- 1 / (nVec * K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 3){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- rep( 1 / K,  K )
    }else if(sampSzWeight == 4){
        # LHS weights
        u1 <- 1 / N # LHS Objective weight
        u2 <- 1 # LHS Objective weight

        # RHS weights
        u3 <- rep(1 / N, K)
        u4 <- nVec / N
    }


    # if intercept in stacking regression add column of ones to beta
    if(stackInt){

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- u1 * eta / 2 * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ),
                                         Y = as.matrix( data[SIndx,2] ),
                                         beta = as.matrix( cbind(1, beta) ) %*% as.vector(w),
                                         weights = weights[SIndx, SIndx])
    }else{

        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- u1 * eta / 2 * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ),
                                         Y = as.matrix( data[SIndx,2] ),
                                         beta = as.matrix( beta ) %*% as.vector(w),
                                         weights = weights)
    }

        for(j in 1:K){

            indx <- which(Study == j) # rows from this study
            sslRows <- intersect( SSLIndx, indx ) # rows from this study that are used for SSL training

            #divide by weights to make objective consider sample sizes
            obj <- obj + (1 - eta) * u3[j] * objLinear(X = as.matrix( cbind(1, data[sslRows,-c(1,2)] ) ),
                                                        Y = as.vector( data[sslRows, 2] ),
                                                        beta = as.matrix(beta[,j],
                                                        weights = 1
                                                        )
            )

        }

        lambdaVec <- lambdaVec * u4 # scale L2 penalty on study specific betas by sample sizes
        # add penalties
        obj <- obj + (1 - eta) * sum( beta^2 %*% lambdaVec ) / 2 + eta * mu * u2 * sum( w^2 ) / 2 # the 2s are to make objective consistent with glmnet

    return(obj)

}

oecRidgePred <- function(mod, data){
    # data is just design matrix--no X or Y
    # mod is OEC ridge model

    # add 1 for intercept
    data <- as.matrix( cbind(1, data) )

    if(length(mod$w) == ncol(mod$beta) + 1 ){
        # if stacking intercept
        preds <- data %*% mod$beta %*% mod$w[-1]# prediction matrix
        preds <- preds + mod$w[1] # add intercept

    }else{

        preds <- data %*% mod$beta %*% mod$w

    }

    return(preds)
}



oecRidgePredX <- function(mod, data){
    # data is just design matrix--no X or Y
    # mod is OEC ridge model

    # add 1 for intercept
    data <- X <- as.matrix( cbind(1, data) )

    K <- ncol(mod$meanMat) # num studies
    preds <- 0

    for(k in 1:K){

        # standardize with study specific means and variances
        for(j in 2:ncol(data)){
            X[,j] <- (data[,j] - mod$meanMat[j - 1,k]) / mod$sigmaMat[j - 1,k]
        }

        if(length(mod$w) == ncol(mod$beta) + 1 ){
            # if stacking intercept add intercept below
            preds <- preds + X %*% mod$beta[,k] * mod$w[k + 1]# prediction matrix


        }else{

            preds <- preds + X %*% mod$beta[,k] * mod$w[k]

        }

    }

    if(length(mod$w) == ncol(mod$beta) + 1 )     preds <- preds + mod$w[1] # add intercept


    return(preds)
}


predsX <- function(meanMat, sigmaMat, data, beta){
    # make design matrix for stacking with standardizing X
    # meanMat and sigmaMat are p x K matrices with varainces and means of covaraites
    # data is just design matrix not with outcome out study labels

    # add 1 for intercept
    data <- X <- as.matrix( cbind(1, data ) )

    K <- ncol(meanMat) # num studies
    preds <- 0
    predsMat <- matrix(nrow = nrow(data), ncol = K )

    for(k in 1:K){

        # standardize with study specific means and variances
        for(j in 2:ncol(data)){
            X[,j] <- (data[,j] - meanMat[j - 1, k]) / sigmaMat[j - 1, k]
        }

        predsMat[, k ] <- preds + X %*% beta[,k]


    }


    return(predsMat)
}



#############
# Ridge Alt -- same as Ridge Alt glmnet but saves stacking weights and betas at each step of optimization
#############

ridgeAltPlot <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 0,
                     Avg = FALSE,
                     wUpdate = "glmnet"){

    # eta is the parameter that determines convex combination of the losses

    if(nnlsInd){
        # nonegative least squares for glmnet paramter
        lowLim <- 0
    }else{
        # No nonegative least squares for glmnet paramter
        lowLim <- -Inf
    }

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

    wMat <- matrix(ncol = length(wStart), nrow = 200) # assumes 200 maximum iterations of optimization
    betaMat <- array(NA, c(200, nrow(beta), ncol(beta) )    )

    wMat[1,] <- wStart
    betaMat[1,,] <- beta

    if(Avg){
        len <- length(wStart[-1])
        w <- rep(1 / len, len)
    }

    if(nnlsInd){
        # if there is a nonnegativity constraint on the stacking weights
        # do not apply to intercept

        w <- sapply(wStart[-1], function(x) max(0, x) )
    }else{

        w <- wStart[-1]

    }

    w_Star <- c(w0, w) # start best w as the initial value (with no negative numbers if needed)

    wOld <- wStart # for update difference calcs
    betaOld <- betaStart

    X <- as.matrix( cbind(1, data[,-c(1,2)] ) )
    y <- as.matrix( data[,2] )
    Study <- data[,1]
    K <- length(unique(Study))
    p <- ncol(X)

    itr <- 1

    #############################
    # Calculate study specific value for below
    #############################
    lambdaList <- list(length = K)
    XX_k_list <- list(length = K)
    Xy_k_list <- list(length = K)
    SSLindx <- c()
    Stackindx <- c()

    for(k in 1:K){
        indx <- which(Study == k)

        if(dataSplit < 1){

            n_k <- length(indx) # n_k
            nSplit <- round(n_k * dataSplit) # n_k split
            indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
            indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking

            # save for below
            Stackindx <- c(Stackindx, indx_k_Stack)

        }else{
            # use all data for both
            indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both

        }



        lambdaList[[k]] <- diag( rep( lambdaVec[k], p )   )
        XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    }

    ##########################################################
    # Calculate stacking values for below
    ##########################################################
    # calculate quantities that will be used below so
    # do not have to repeatedly calculate in while loop

    # if dataSplit is 1 then use all the data for stacking. First set to all the rows then below reset to NULL
    if(dataSplit == 1)    Stackindx <- 1:nrow(X)


    XX <- t(X[Stackindx,]) %*% X[Stackindx,] # X^TX calculate so do not have to at each iteration
    Xy <- t(X[Stackindx,]) %*% y[Stackindx] # X^T y
    X_rowSums <- rowSums( t(X[Stackindx,]) )
    X <- X[Stackindx,]
    Xnrow <- nrow(X)
    y_sum <- sum(y[Stackindx])
    muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]

    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    if(dataSplit == 1)    Stackindx <- NULL

    objVec <- objOEC(data = data,
                     beta = betaStart,
                     w = wStart,
                     mu = mu,
                     lambdaVec = lambdaVec,
                     stackInt = TRUE,
                     eta = eta,
                     StackIndx = Stackindx)

    #############################
    # Iteratively update
    #############################

    while(eps > tol){

        #print(paste("iter", itr))
        itr <- itr + 1
        # beta update
        for(k in 1:K){

            inv <- solve(
                eta * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] + lambdaList[[k]] )
            )

            beta[,k] <- inv %*% (
                eta * w[k] * ( Xy - w0 * X_rowSums  - XX %*% beta[,-k] %*% w[-k] ) +
                    (1 - eta)  *  Xy_k_list[[k]]
            )
        }

        # w update
        if(!Avg){
            # if not using average weights

            if(wUpdate == "glmnet"){
                # if use glmnet to update parameters
                mod <- glmnet(y = y,
                              x = as.matrix(X %*% beta),
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = TRUE,
                              intercept = TRUE)

                # coefficients
                w <-  as.vector(mod$beta)
                w0 <- as.vector(mod$a0)

            }else{

                # w0 update
                w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow

                # w update
                bInv <- solve(    eta * ( t(beta) %*% XX %*% beta + muMat )    )
                w <- bInv %*% (
                    eta * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) ) )

                )



            }

        }else{
            # if average just update the intercept with closed form expression
            # w0 update
            w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow

        }

        if(nnlsInd){
            # if there is a nonnegativity constraint on the stacking weights
            # do not apply to intercept

            if(projs > 0){

                # if project onto unit simplex

                for(pr in 1:projs){
                    # project onto unit simplex

                    # Projection onto Non Negative Orthant
                    w <- sapply(w, function(x) max(0, x) )
                    # Projection onto Sum of 1
                    w <- w - ((sum(w) - 1) / K)

                    # Projection onto Non Negative Orthant
                    w <- sapply(w, function(x) max(0, x) )
                }

            }else{
                # just standard projection onto non negative orthant -- not onto unit simplex

                # Projection onto Non Negative Orthant
                w <- sapply(w, function(x) max(0, x) )
            }


        }

        obj_r <- objOEC(data = data,
                        beta = beta,
                        w = c(w0, w),
                        mu = mu,
                        lambdaVec = lambdaVec,
                        stackInt = TRUE,
                        eta = eta,
                        StackIndx = Stackindx)



        if(obj_r < objVec[itr - 1]){
            # if the objective improves save it as current best
            beta_Star <- beta
            w_Star <- c(w0, w)
        }else{
            # if it is not better, then do not update the best

            if(objCriteria == TRUE){
                # if obj criteria is true, use the best iteration as this iterations value
                w <- w_Star[-1]
                w0 <- w_Star[1]
                beta <- beta_Star
            }
        }

        # calculate change from last iteration
        eps <- sum( (betaOld - beta)^2  ) + sum( (wOld - c(w0, w))^2  )

        objVec <- c(objVec, obj_r )
        # update previous iteration values
        betaOld <- beta
        wOld <- c(w0, w)

        # save values
        wMat[itr,] <- wOld
        betaMat[itr,,] <- betaOld
    }



    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, betaMat = betaMat, wMat = wMat))

}

##################
# Ridge OEC w Tuner
##################
# tunes l2 hyperparameter associated with stacking regression
oecW_HOSO <- function(data,
                   tune.grid,
                   sslLambdas,
                   method = "glmnet",
                   nfolds = "K",
                   nnlsInd = TRUE,
                   OECeta = 0.5,
                   sampSzWeight = 4,
                   standardize = TRUE,
                   weights = NULL,
                   glmnetOpt = TRUE,
                   xStandardize = FALSE){

    # sslLambdas is the hyperparameters for the SSL ridge terms
    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation

    library(caret)
    library(glmnet)

    # rename studies from 1:K
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    rm(indx, studyVec, indxVec)

    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems

        source("OEC Functions.R")

        # set number of folds to number of training studies unless
        # specified as another number

        if(nfolds == "K" || nfolds > K)       nfolds <- K
        if( is.null(weights) )   weights <- diag( nrow(data) ) # identity

        resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets

        for(study in 1:nfolds){

            HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){

                print(paste0("study ", study, "_tune_", tune))

                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas[-study], # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack
                stackParam <- warmRes$mu
                tuneParam <- sslLambdas[-study] #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies

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



        resMean <- colMeans(resMat)

        testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE

        }

        return(best = testLambda)


}


##################
# ridge estimator
##################
ridgeEst <- function(y,
                     x,
                     lambda,
                     intercept = TRUE,
                     coefWeights = NULL#,
                     # xStandardize = FALSE # standardize covariates
                    ){

    if(intercept){
        #oneVec <- rep(1, nrow(x)) # for intercept
        x <- as.matrix( cbind(1, x) )
        if(is.null(coefWeights)){
            diagMat <- as.matrix( diag( ncol(x) ) )
            diagMat[1,1] <- 0 # do not penalize intercept
        }else{
            diagMat <- as.matrix( diag(coefWeights ) )
        }

    }else{
        x <- as.matrix(x)
        if(is.null(coefWeights)){
            diagMat <- as.matrix( diag( ncol(x) ) )
        }else{
            diagMat <- as.matrix( diag(coefWeights ) )
        }
    }

    y <- as.vector( y )
    n <- length(y)
    lambda <- as.numeric(lambda)

    # ridge estimator
    cfs <- solve( t( x ) %*% x + n * lambda * diagMat ) %*% (t(x) %*% y ) # original y (not scaled)

    return(cfs)

}

##################
# Weighted ridge estimator
# penalize betas based on matrix of weights
##################
ridgeEstW <- function(y,
                     x,
                     lambda,
                     intercept = TRUE,
                     coefWeights = NULL#,
                     # xStandardize = FALSE # standardize covariates
){

    if(intercept){
        #oneVec <- rep(1, nrow(x)) # for intercept
        x <- as.matrix( cbind(1, x) )


    }else{
        x <- as.matrix(x)

    }

    y <- as.vector( y )
    n <- length(y)
    lambda <- as.numeric(lambda)

    # ridge estimator
    cfs <- solve( t( x ) %*% x + n * lambda * coefWeights ) %*% (t(x) %*% y ) # original y (not scaled)

    return(cfs)

}


#################
# multi Penalty Tuner
#################

multiPenTune <- function(data,
                         tune.grid,
                         penalty = 1,
                         nfolds = 10,
                         sslLambdas,
                         weights = NULL){

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
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems

    source("OEC Functions.R")

    # set number of folds to number of training studies unless
    # specified as another number

    if(nfolds == "K" || nfolds > K)       nfolds <- K
    if( is.null(weights) )   weights <- Diagonal( nrow(data) ) # identity

    resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets

    foldsList <- createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

    for(study in 1:nfolds){

        HOOList <- foldsList[[study]] # indices of study to hold out
        indxList <- do.call(c, foldsList[-study]) # concatenate indices of studies to train on

        #standrize weights

        for(tune in 1:nrow(tune.grid) ){


            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = tune.grid$lambda[tune],
                               nnlsInd = TRUE,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = TRUE,
                               stackStandardize = TRUE,
                               glmnetOpt = TRUE,
                               xStandardize = FALSE,
                               sampSzWeight = 1,
                               weights = weights[indxList,indxList])

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack
            stackParam <- tune.grid$lambda[tune]
            tuneParam <- sslLambdas
            betas <- warmRes$beta

            if(penalty > 4){
                pen <- penalty
                penalty <- penalty - 4
            }else{
                pen <- penalty
            }

            if(penalty == 1){

                betaCov <- cov(t(betas))

            }else if(penalty == 2){

                betaCov <- cov(t(betas))

                # remove intercept terms
                betaCov[1,] <- 0
                betaCov[,1] <- 0

            }else if(penalty == 3){

                betaCov <- cov(t(betas))
                betaCovDiag <- diag( diag(betaCov) )

            }else if(penalty == 4){

                betaCov <- cov(t(betas))

                # remove intercept terms
                betaCov[1,] <- 0
                betaCov[,1] <- 0
                betaCov <- diag( diag(betaCov) )

            }

            if(pen > 4){
                betaCov <- solve(betaCov)
            }

            oecRidgeWS <- ridgeEstW(data$Y[indxList],
                                                x = data[indxList, -c(1,2)],
                                                lambda = tune.grid$lambda[tune],
                                                intercept = TRUE,
                                                coefWeights = betaCov#,
                                                # xStandardize = FALSE # standardize covariates
                                                )

            predsVecOEC <- cbind(1, as.matrix( data[HOOList, -c(1,2)]) ) %*% oecRidgeWS

            resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

        }

    }

    # moved outside for loop 9/7/20
    resMean <- colMeans(resMat)

    testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE

    return(list(lambdaOpt = testLambda, res = resMat))
}

##################
# Ridge OEC w Tuner
##################
# tunes l2 hyperparameter associated with stacking regression but with balanced groups with respect to study
# NOT a hold one study out -- use this so it has the same number of studies
oecW_CV <- function(data,
                    tune.grid,
                    sslLambdas,
                    method = "glmnet",
                    nfolds = "K",
                    nnlsInd = TRUE,
                    OECeta = 0.5,
                    sampSzWeight = 4,
                    cv = "cv",
                    standardize = TRUE,
                    weights = NULL,
                    glmnetOpt = TRUE,
                    xStandardize = FALSE,
                    SpecID = NULL,
                    horizon = 10){

    # sslLambdas is the hyperparameters for the SSL ridge terms
    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation
    # cv is either "cv" for balanced folds, or "hoso" for hold one study out

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
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)

    source("OEC Functions.R")

    # set number of folds to number of training studies unless
    # specified as another number

    if(nfolds == "K" || nfolds > K)       nfolds <- K
    if( is.null(weights) )   weights <- Diagonal( nrow(data) ) # identity

    resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets
    allRows <- 1:nrow(data)

    if(cv == "cv"){
        message(cv)
        foldsList <- createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

        for(study in 1:nfolds){

            HOOList <- foldsList[[study]] # indices of study to hold out
            indxList <- allRows[-HOOList] #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack
                stackParam <- tune.grid$lambda[tune]
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
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
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        # moved outside for loop 9/7/20
        resMean <- colMeans(resMat)

        testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE


    }else if(cv == "hoso"){
        message(cv)
        for(study in 1:nfolds){

            HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas[-study], # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack
                stackParam <- warmRes$mu
                tuneParam <- sslLambdas[-study] #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
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
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }


        }

        resMean <- colMeans(resMat)
        testLambda <- tune.grid[ which.min(resMean),   ] # choose hyperparameters with best RMSE

    }else if(cv == "cvSpecTS"){
        message(cv)

        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # in case there is not enough in the horizon, choose the minimum
        horizon2 <- min(horizon,
                       length(indx1) - round(0.5 * length(indx1))
                       )

        # time series folds generator
        foldsList <- createTimeSlices(indx1,
                                      initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                      horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                      fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
        ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

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


        resMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

        # remove some so we do not have so many
        foldsList$test <- foldsList$test[seqVec]
        foldsList$train <- foldsList$train[seqVec]

        allRows <- seq(1, nrow(data))

        for(study in 1:len){

            HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
            specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

            remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
            specIndx <- indx1[specIndx] # make these in terms of actual row numbers

            indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = tune.grid$lambda[tune],
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = OECeta,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }



        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }else if(cv == "cvSpecTS0"){
        message(cv)

        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # in case there is not enough in the horizon, choose the minimum
        horizon2 <- min(horizon,
                        length(indx1) - round(0.5 * length(indx1))
        )

        # time series folds generator
        foldsList <- createTimeSlices(indx1,
                                      initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                      horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                      fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
        )  # make folds with balanced studies

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



        resMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

        # remove some so we do not have so many
        foldsList$test <- foldsList$test[seqVec]
        foldsList$train <- foldsList$train[seqVec]

        allRows <- seq(1, nrow(data))

        for(study in 1:len){

            HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
            specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

            remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
            specIndx <- indx1[specIndx] # make these in terms of actual row numbers

            indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt # ,
                                   # zeroOut = TRUE # decided not to include since its not a zero out Specialist ridgeWS
                                   )

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = tune.grid$lambda[tune],
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = OECeta,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }



        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }else if(cv == "cvSpecTSHOO"){
        # iterates through all specialists
        message(cv)

        resMat2 <- matrix(nrow = K, ncol = nrow(tune.grid)) # store RMSE of held out sets

        for(specID in 1:K){

            # specialist CV
            indx1 <- which(data$Study == specID)

            # in case there is not enough in the horizon, choose the minimum
            horizon2 <- min(horizon,
                            length(indx1) - round(0.5 * length(indx1))
            )

            # time series folds generator
            foldsList <- createTimeSlices(indx1,
                                          initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                          horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                          fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
            ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

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


            resMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

            # remove some so we do not have so many
            foldsList$test <- foldsList$test[seqVec]
            foldsList$train <- foldsList$train[seqVec]

            allRows <- seq(1, nrow(data))

            for(study in 1:len){

                HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
                specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

                remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
                specIndx <- indx1[specIndx] # make these in terms of actual row numbers

                indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

                #standrize weights
                wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

                for(tune in 1:nrow(tune.grid) ){
                    set.seed(1)
                    warmRes <- ridgeWS(data = data[indxList,], # current studies
                                       tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                       stackParam = tune.grid$lambda[tune],
                                       nnlsInd = nnlsInd,
                                       stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                       modelStandardize = standardize,
                                       stackStandardize = standardize,
                                       glmnetOpt = glmnetOpt,
                                       xStandardize = xStandardize,
                                       sampSzWeight = sampSzWeight,
                                       weights = wgt)

                    sigK <- warmRes$sigK
                    sigStack <- warmRes$sigStack


                    stackParam <- warmRes$mu
                    tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam,
                                               mu = tune.grid$lambda[tune],
                                               nnlsInd = nnlsInd,
                                               Stackindx = specIndx, # rows corresponding to test country
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               weights = wgt)

                    predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                mod = oecRidgeWS)

                    resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

                }


            }

            resMat2[specID, ] <- colMeans(resMat, na.rm = TRUE)

            }

            resMean <- colMeans(resMat2, na.rm = TRUE)
            testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }else if(cv == "cvSpec0TSHOO"){
        # iterates through all specialists
        message(cv)

        resMat2 <- matrix(nrow = K, ncol = nrow(tune.grid)) # store RMSE of held out sets

        for(specID in 1:K){

            # specialist CV
            indx1 <- which(data$Study == specID)

            # in case there is not enough in the horizon, choose the minimum
            horizon2 <- min(horizon,
                            length(indx1) - round(0.5 * length(indx1))
            )

            # time series folds generator
            foldsList <- createTimeSlices(indx1,
                                          initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                          horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                          fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
            ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

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


            resMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

            # remove some so we do not have so many
            foldsList$test <- foldsList$test[seqVec]
            foldsList$train <- foldsList$train[seqVec]

            allRows <- seq(1, nrow(data))

            for(study in 1:len){

                HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
                specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

                remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
                specIndx <- indx1[specIndx] # make these in terms of actual row numbers

                indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

                #standrize weights
                wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

                for(tune in 1:nrow(tune.grid) ){
                    set.seed(1)
                    warmRes <- ridgeWS(data = data[indxList,], # current studies
                                       tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                       stackParam = tune.grid$lambda[tune],
                                       nnlsInd = nnlsInd,
                                       stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                       modelStandardize = standardize,
                                       stackStandardize = standardize,
                                       glmnetOpt = glmnetOpt,
                                       xStandardize = xStandardize,
                                       sampSzWeight = sampSzWeight,
                                       weights = wgt)

                    sigK <- warmRes$sigK
                    sigStack <- warmRes$sigStack

                    # print(paste0("study ", study, "_tune_", tune))

                    stackParam <- warmRes$mu
                    tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam,
                                               mu = tune.grid$lambda[tune],
                                               nnlsInd = nnlsInd,
                                               Stackindx = specIndx, # rows corresponding to test country
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = OECeta,
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               weights = wgt)

                    predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                mod = oecRidgeWS)

                    resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

                }


            }

            resMat2[specID, ] <- colMeans(resMat, na.rm = TRUE)

        }

        resMean <- colMeans(resMat2, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }else if(cv == "hosoSpecTS"){
        message(cv)

        for(study in 1:nfolds){

            HOOList <- which(data$Study == study)
            len <- length(HOOList)

            endIndx <- seq(1, round(0.9 * len) ) # take the first 90% as training for specialist country
            specIndx <- HOOList[endIndx] # indices of portion to train on for specialist specifically
            HOOList <- HOOList[-endIndx] # test set
            indxList <- seq(1, nrow(data))[-HOOList]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = stackParam,
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = tune.grid$lambda[tune],
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = OECeta,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }




        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE

    }else if(cv == "cvSpec"){
        message(cv)
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # time series folds generator
        foldsList <- createFolds(factor(data$Study)[indx1], k = nfolds) # make folds with balanced studies

        # resMat <- matrix(nrow = nfolds, ncol = nrow(tune.grid)) # store RMSE of held out sets

        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            specIndx <- which(data$Study[indxList] == SpecID)
            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = tune.grid$lambda[tune],
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = OECeta,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }



        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }else if(cv == "cvSpec0"){
        message(cv)
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # time series folds generator
        foldsList <- createFolds(factor(data$Study)[indx1], k = nfolds) # make folds with balanced studies

        # resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets

        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on
            specIndx <- which(data$Study[indxList] == SpecID)

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])

            for(tune in 1:nrow(tune.grid) ){
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = tune.grid$lambda[tune],
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt,
                                   zeroOut = TRUE)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = tune.grid$lambda[tune],
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = OECeta,
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }



        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE



    }


    return(best = testLambda)


}

stackCVTS <- function(data,
                      tune.grid,
                      nnls = TRUE,
                      nfolds = 10,
                      horizon= 10
                      ){

        if(nnls)      low <- 0  # nnls for stacking


        # specialist CV
        indx1 <- 1:nrow(data)

        # in case there is not enough in the horizon, choose the minimum
        horizon2 <- min(horizon,
                        length(indx1) - round(0.5 * length(indx1))
        )

        # time series folds generator
        foldsList <- createTimeSlices(indx1,
                                      initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                      horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                      fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
        )

        # remove some so we dont do a HO-observation out CV

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


        resMat <- matrix(nrow = len, ncol = nrow(tune.grid)) # store RMSE of held out sets

        # remove some so we do not have so many
        foldsList$test <- foldsList$test[seqVec]
        foldsList$train <- foldsList$train[seqVec]

        for(study in 1:len){

            HOOList <- foldsList$test[[study]] # indices of portion to hold out
            specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

            #standrize weights

            for(tune in 1:nrow(tune.grid) ){

                mod <- glmnet(y = as.vector(data$Y[specIndx]),
                              x = as.matrix(data[specIndx,-1]),
                              alpha = tune.grid$alpha[tune],
                              lambda = tune.grid$lambda[tune],
                              standardize = TRUE,
                              intercept = TRUE,
                              lower.limits = low )

                predsVecOEC <- predict(mod, as.matrix( data[HOOList, -1]) )

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }



        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        testLambda <- tune.grid[ which.min(resMean),  ] # choose hyperparameters with best RMSE

    return(testLambda)

}


#################
# Ridge OEC w Tuner
##################
# tunes l2 hyperparameter associated with stacking regression but with balanced groups with respect to study
# NOT a hold one study out -- use this so it has the same number of studies
oecEta_CV <- function(data,
                    tune.grid,
                    sslLambdas,
                    stackParam = 0,
                    method = "glmnet",
                    nfolds = "K",
                    nnlsInd = TRUE,
                    sampSzWeight = 4,
                    cv = "cv",
                    standardize = TRUE,
                    weights = NULL,
                    glmnetOpt = TRUE,
                    xStandardize = FALSE,
                    AvgW = FALSE,
                    SpecID = NULL,
                    horizon = 10
                    ){

    #glmnetOpt is for ridgeWS
    # xStandardize is for ridgeWS
    # SpecID is study number of test study

    # sslLambdas is the hyperparameters for the SSL ridge terms
    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation
    # cv is either "cv" for balanced folds, or "hoso" for hold one study out

    library(caret)
    library(glmnet)

    # rename studies from 1:K
    num.trainStudy <- K <- length(unique(data$Study))
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    rm(indx, studyVec, indxVec)

    source("OEC Functions.R")

    # set number of folds to number of training studies unless
    # specified as another number

    if(nfolds == "K" || nfolds > K)       nfolds <- K
    if( is.null(weights) )   weights <- Diagonal( nrow(data) ) # identity

    resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets
    allRows <- 1:nrow(data)

    if(cv == "cv"){

        foldsList <- createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies

        for(study in 1:nfolds){

            HOOList <- foldsList[[study]] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas  # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAlt(data = data[indxList,], # current studies
                                       betaStart = warmRes$beta,
                                       wStart = warmRes$w,
                                       lambdaVec = tuneParam,
                                       mu = stackParam,
                                       nnlsInd = nnlsInd,
                                       tol = 0.01,
                                       objCriteria = TRUE,
                                       eta = tune.grid[tune],
                                       projs = 0,
                                       Avg = AvgW,
                                       wUpdate = "glment",
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        # moved outside 9/7/20
        resMean <- colMeans(resMat)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "hoso"){

        for(study in 1:nfolds){

            HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas[-study], # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            stackParam <- warmRes$mu
            tuneParam <- sslLambdas[-study] # use the same parameter for all studies

            for(tune in 1:length(tune.grid) ){

                set.seed(1)
                oecRidgeWS <- ridgeAlt(data = data[indxList,], # current studies
                                       betaStart = warmRes$beta,
                                       wStart = warmRes$w,
                                       lambdaVec = tuneParam,
                                       mu = stackParam,
                                       nnlsInd = nnlsInd,
                                       tol = 0.01,
                                       objCriteria = TRUE,
                                       eta = tune.grid[tune],
                                       projs = 0,
                                       Avg = AvgW,
                                       wUpdate = "glment",
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

            resMean <- colMeans(resMat)
            etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvSpecTS"){
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # in case there is not enough in the horizon, choose the minimum
        horizon2 <- min(horizon,
                        length(indx1) - round(0.5 * length(indx1))
        )

        # time series folds generator
        foldsList <- createTimeSlices(indx1,
                                      initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
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



        resMat <- matrix(nrow = len, ncol = length(tune.grid)) # store RMSE of held out sets

        # remove some so we do not have so many
        foldsList$test <- foldsList$test[seqVec]
        foldsList$train <- foldsList$train[seqVec]
        allRows <- 1:nrow(data)

        for(study in 1:len){

            HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
            specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

            remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
            specIndx <- indx1[specIndx] # indices in terms of actual rows of dataset
            indxList <- allRows[-remInd]    # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas  # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                       betaStart = warmRes$beta,
                                       wStart = warmRes$w,
                                       lambdaVec = tuneParam,
                                       mu = stackParam,
                                       nnlsInd = nnlsInd,
                                       Stackindx = specIndx, # rows corresponding to test country
                                       tol = 0.01,
                                       objCriteria = TRUE,
                                       eta = tune.grid[tune],
                                       projs = 0,
                                       Avg = AvgW,
                                       wUpdate = "glment",
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvSpecTSHOO"){

        resMat2 <- matrix(nrow = K, ncol = length(tune.grid)) # store RMSE of held out sets

        for(SpecID in 1:K){
            # specialist CV
            indx1 <- which(data$Study == SpecID)

            # in case there is not enough in the horizon, choose the minimum
            horizon2 <- min(horizon,
                            length(indx1) - round(0.5 * length(indx1))
            )

            # time series folds generator
            foldsList <- createTimeSlices(indx1,
                                          initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                          horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                          fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
            ) # make folds with balanced studies

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



            resMat <- matrix(nrow = len, ncol = length(tune.grid)) # store RMSE of held out sets

            # remove some so we do not have so many
            foldsList$test <- foldsList$test[seqVec]
            foldsList$train <- foldsList$train[seqVec]
            allRows <- 1:nrow(data)

            for(study in 1:len){

                HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
                specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

                remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
                specIndx <- indx1[specIndx] # indices in terms of actual rows of dataset
                indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

                #standrize weights
                wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = stackParam,
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                for(tune in 1:length(tune.grid) ){

                    stackParam <- warmRes$mu
                    tuneParam <- sslLambdas # use the same parameter for all studies
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam,
                                               mu = stackParam,
                                               nnlsInd = nnlsInd,
                                               Stackindx = specIndx, # rows corresponding to test country
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = tune.grid[tune],
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               weights = wgt)

                    predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                mod = oecRidgeWS)

                    resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

                }
            }

            resMat2[SpecID, ] <- colMeans(resMat, na.rm = TRUE)

        }

        resMean <- colMeans(resMat2, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvSpecTS0"){
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # in case there is not enough in the horizon, choose the minimum
        horizon2 <- min(horizon,
                        length(indx1) - round(0.5 * length(indx1))
        )

        # time series folds generator
        foldsList <- createTimeSlices(indx1,
                                      initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                      horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                      fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
        ) # make folds with balanced studies

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



        resMat <- matrix(nrow = len, ncol = length(tune.grid)) # store RMSE of held out sets

        # remove some so we do not have so many
        foldsList$test <- foldsList$test[seqVec]
        foldsList$train <- foldsList$train[seqVec]
        allRows <- 1:nrow(data)

        for(study in 1:len){

            HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
            specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

            remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
            specIndx <- indx1[specIndx] # indices in terms of actual rows of dataset
            indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt #,
                               # zeroOut = FALSE, # decided to do false since it is not a specialist tuning
                               )

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){


                stackParam <- warmRes$mu
                tuneParam <- sslLambdas # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = stackParam,
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = tune.grid[tune],
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvSpec0TSHOO"){

        resMat2 <- matrix(nrow = K, ncol = length(tune.grid)) # store RMSE of held out sets

        for(SpecID in 1:K){
            # specialist CV
            indx1 <- which(data$Study == SpecID)

            # in case there is not enough in the horizon, choose the minimum
            horizon2 <- min(horizon,
                            length(indx1) - round(0.5 * length(indx1))
            )

            # time series folds generator
            foldsList <- createTimeSlices(indx1,
                                          initialWindow = round(0.5 * length(indx1)), # minimum training set is half of total observations
                                          horizon = horizon2, # number of test-points in each fold #round(0.5 * length(indx1) / nfolds ), #nfolds, #round(0.5 * length(indx1) / nfolds ), # makes k = nfolds functionally
                                          fixedWindow = FALSE # if FALSE, the training set always start at the first sample and the training set size will vary over data splits.
            ) # make folds with balanced studies

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



            resMat <- matrix(nrow = len, ncol = length(tune.grid)) # store RMSE of held out sets

            # remove some so we do not have so many
            foldsList$test <- foldsList$test[seqVec]
            foldsList$train <- foldsList$train[seqVec]
            allRows <- 1:nrow(data)

            for(study in 1:len){

                HOOList <- indx1[ foldsList$test[[study]] ] # indices of portion to hold out
                specIndx <- foldsList$train[[study]] # indices of portion to train on for specialist specifically

                remInd <- indx1[-specIndx] # indices to remove from training set (because we are removing everything from before this point in test country's trainign set)
                specIndx <- indx1[specIndx] # indices in terms of actual rows of dataset
                indxList <- allRows[-remInd]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

                #standrize weights
                wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
                set.seed(1)
                warmRes <- ridgeWS(data = data[indxList,], # current studies
                                   tuneParam = sslLambdas, # current lambda associated with L2 penalty
                                   stackParam = stackParam,
                                   nnlsInd = nnlsInd,
                                   stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                                   modelStandardize = standardize,
                                   stackStandardize = standardize,
                                   glmnetOpt = glmnetOpt,
                                   xStandardize = xStandardize,
                                   sampSzWeight = sampSzWeight,
                                   weights = wgt)

                sigK <- warmRes$sigK
                sigStack <- warmRes$sigStack

                for(tune in 1:length(tune.grid) ){

                    stackParam <- warmRes$mu
                    tuneParam <- sslLambdas  # use the same parameter for all studies
                    set.seed(1)
                    oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                               betaStart = warmRes$beta,
                                               wStart = warmRes$w,
                                               lambdaVec = tuneParam,
                                               mu = stackParam,
                                               nnlsInd = nnlsInd,
                                               Stackindx = specIndx, # rows corresponding to test country
                                               tol = 0.01,
                                               objCriteria = TRUE,
                                               eta = tune.grid[tune],
                                               projs = 0,
                                               Avg = AvgW,
                                               wUpdate = "glment",
                                               standardize = TRUE,
                                               sampSzWeight = sampSzWeight,
                                               sigK = sigK,
                                               sigStack = sigStack,
                                               weights = wgt)

                    predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                                mod = oecRidgeWS)

                    resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

                }
            }

            resMat2[SpecID, ] <- colMeans(resMat, na.rm = TRUE)

        }

        resMean <- colMeans(resMat2, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "hosoSpecTS"){


        for(study in 1:nfolds){

            HOOList <- which(data$Study == study)
            len <- length(HOOList)

            endIndx <- seq(1, round(0.9 * len) ) # take the first 90% as training for specialist country
            specIndx <- HOOList[endIndx] # indices of portion to train on for specialist specifically
            HOOList <- HOOList[-endIndx] # test set
            indxList <- seq(1, nrow(data))[-HOOList]   #do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = stackParam,
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = tune.grid[tune],
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE

    }else if(cv == "cvSpec"){
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # time series folds generator
        foldsList <- createFolds(factor(data$Study)[indx1], k = nfolds) # make folds with balanced studies

        resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets

        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            specIndx <- which(data$Study[indxList] == SpecID)
            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = stackParam,
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = tune.grid[tune],
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvSpec0"){
        # specialist CV
        indx1 <- which(data$Study == SpecID)

        # time series folds generator
        foldsList <- createFolds(factor(data$Study)[indx1], k = nfolds) # make folds with balanced studies

        resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets

        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on
            specIndx <- which(data$Study[indxList] == SpecID)

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt,
                               zeroOut = TRUE)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                           betaStart = warmRes$beta,
                                           wStart = warmRes$w,
                                           lambdaVec = tuneParam,
                                           mu = stackParam,
                                           nnlsInd = nnlsInd,
                                           Stackindx = specIndx, # rows corresponding to test country
                                           tol = 0.01,
                                           objCriteria = TRUE,
                                           eta = tune.grid[tune],
                                           projs = 0,
                                           Avg = AvgW,
                                           wUpdate = "glment",
                                           standardize = TRUE,
                                           sampSzWeight = sampSzWeight,
                                           sigK = sigK,
                                           sigStack = sigStack,
                                           weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat, na.rm = TRUE)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }

    return(best = etaStar)


}

#################
# Lower Limits OEC w Tuner
##################
# tunes lower limits
# NOT a hold one study out -- use this so it has the same number of studies
oecLow_CV <- function(data,
                      tune.grid,
                      sslLambdas,
                      stackParam = 0,
                      oecEta = 0.99,
                      method = "glmnet",
                      nfolds = "K",
                      nnlsInd = TRUE,
                      sampSzWeight = 4,
                      cv = "cv",
                      standardize = TRUE,
                      weights = NULL,
                      glmnetOpt = TRUE,
                      xStandardize = FALSE,
                      AvgW = FALSE){

    #glmnetOpt is for ridgeWS
    # xStandardize is for ridgeWS

    # sslLambdas is the hyperparameters for the SSL ridge terms
    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation
    # cv is either "cv" for balanced folds, or "hoso" for hold one study out

    library(caret)
    library(glmnet)

    # rename studies from 1:K
    num.trainStudy <- K <- length(unique(data$Study))
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    rm(indx, studyVec, indxVec)

    source("OEC Functions.R")

    # set number of folds to number of training studies unless
    # specified as another number

    if(nfolds == "K" || nfolds > K)       nfolds <- K
    if( is.null(weights) )   weights <- Diagonal( nrow(data) ) # identity

    resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets

    if(cv == "cv"){

        foldsList <- createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies
        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- foldsList[[study]] # indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltFix(data = data[indxList,], # current studies
                                       betaStart = warmRes$beta,
                                       wStart = warmRes$w,
                                       lambdaVec = tuneParam,
                                       mu = stackParam,
                                       tol = 0.01,
                                       nnlsInd = tune.grid[tune] / K,
                                       low = tune.grid[tune] / K,
                                       up = Inf,
                                       objCriteria = TRUE,
                                       eta = oecEta,
                                       projs = 0,
                                       Avg = AvgW,
                                       wUpdate = "glment",
                                       standardize = TRUE,
                                       sampSzWeight = sampSzWeight,
                                       sigK = sigK,
                                       sigStack = sigStack,
                                       weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        resMean <- colMeans(resMat)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "hoso"){

        for(study in 1:nfolds){

            HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas[-study], # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            stackParam <- warmRes$mu
            tuneParam <- sslLambdas[-study] #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies

            for(tune in 1:length(tune.grid) ){

                set.seed(1)
                oecRidgeWS <- ridgeAltFix(data = data[indxList,], # current studies
                                          betaStart = warmRes$beta,
                                          wStart = warmRes$w,
                                          lambdaVec = tuneParam,
                                          mu = stackParam,
                                          tol = 0.01,
                                          nnlsInd = tune.grid[tune] / K,
                                          low = tune.grid[tune] / K,
                                          up = Inf,
                                          objCriteria = TRUE,
                                          eta = oecEta,
                                          projs = 0,
                                          Avg = AvgW,
                                          wUpdate = "glment",
                                          standardize = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          sigK = sigK,
                                          sigStack = sigStack,
                                          weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat)

        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }


    return(best = etaStar)


}



#################
# Upper Limits OEC w Tuner
##################
# tunes lower limits
# NOT a hold one study out -- use this so it has the same number of studies
oecHigh_CV <- function(data,
                      tune.grid,
                      sslLambdas,
                      stackParam = 0,
                      oecEta = 0.99,
                      method = "glmnet",
                      nfolds = "K",
                      nnlsInd = TRUE,
                      sampSzWeight = 4,
                      cv = "cv",
                      standardize = TRUE,
                      weights = NULL,
                      glmnetOpt = TRUE,
                      xStandardize = FALSE,
                      SpecID = NULL,
                      AvgW = FALSE){

    #glmnetOpt is for ridgeWS
    # xStandardize is for ridgeWS

    # sslLambdas is the hyperparameters for the SSL ridge terms
    # nnlsind only relevant if oec is chosen as hoso
    # OEC eta is only relevant if oec is chosen as hoso -- trade off in OEC loss
    # cvFolds is standard cross validation
    # cv is either "cv" for balanced folds, or "hoso" for hold one study out

    library(caret)
    library(glmnet)

    # rename studies from 1:K
    num.trainStudy <- K <- length(unique(data$Study))
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    num.trainStudy <- K <- length(studyVec)
    indxVec <- vector( length = nrow(data) )

    for(j in 1:K){

        indx <- which(data$Study == studyVec[j])
        indxVec[indx] <- j

    }

    data$Study <- indxVec # replace with new study labels from 1:K
    trainStudies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
    rm(indx, studyVec, indxVec)

    source("OEC Functions.R")

    # set number of folds to number of training studies unless
    # specified as another number

    if(nfolds == "K" || nfolds > K)       nfolds <- K
    if( is.null(weights) )   weights <- Diagonal( nrow(data) ) # identity

    resMat <- matrix(nrow = nfolds, ncol = length(tune.grid)) # store RMSE of held out sets

    if(cv == "cv"){

        indx1 <- which(data$Study == SpecID)
        foldsList <- createFolds(indx1, k = nfolds) # make folds with balanced studies
        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ]# indices of study to hold out
            indxList <- allRows[-HOOList] #  # concatenate indices of studies to train on
            indx2 <- which(data$Study[indxList] == SpecID)#

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas # # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpec0(data = data[indxList,], # current studies
                                          betaStart = warmRes$beta,
                                          wStart = warmRes$w,
                                          lambdaVec = tuneParam,
                                          mu = stackParam,
                                          Stackindx = indx2, # rows corresponding to test country
                                          tol = 0.01,
                                          nnlsInd = TRUE,
                                          low = 0,
                                          up = tune.grid[tune],
                                          objCriteria = TRUE,
                                          eta = oecEta,
                                          projs = 0,
                                          Avg = AvgW,
                                          wUpdate = "glment",
                                          standardize = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          sigK = sigK,
                                          sigStack = sigStack,
                                          weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        resMean <- colMeans(resMat)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "hoso"){

        for(study in 1:nfolds){

            HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas[-study], # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            stackParam <- warmRes$mu
            tuneParam <- sslLambdas[-study] ## use the same parameter for all studies

            for(tune in 1:length(tune.grid) ){

                set.seed(1)
                oecRidgeWS <- ridgeAltFix(data = data[indxList,], # current studies
                                          betaStart = warmRes$beta,
                                          wStart = warmRes$w,
                                          lambdaVec = tuneParam,
                                          mu = stackParam,
                                          tol = 0.01,
                                          nnlsInd = tune.grid[tune] / K,
                                          low = tune.grid[tune] / K,
                                          up = Inf,
                                          objCriteria = TRUE,
                                          eta = oecEta,
                                          projs = 0,
                                          Avg = AvgW,
                                          wUpdate = "glment",
                                          standardize = TRUE,
                                          sampSzWeight = sampSzWeight,
                                          sigK = sigK,
                                          sigStack = sigStack,
                                          weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }
        }

        resMean <- colMeans(resMat)

        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }else if(cv == "cvWindow"){

        indx1 <- which(data$Study == SpecID)
        foldsList <- createFolds(indx1, k = nfolds) # make folds with balanced studies
        allRows <- 1:nrow(data)

        for(study in 1:nfolds){

            HOOList <- indx1[ foldsList[[study]] ]# indices of study to hold out
            indxList <- allRows[-HOOList] # do.call(c, foldsList[-study]) # concatenate indices of studies to train on
            indx2 <- which(data$Study[indxList] == SpecID)#indx1[ -foldsList[[study]] ]

            #standrize weights
            wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            set.seed(1)
            warmRes <- ridgeWS(data = data[indxList,], # current studies
                               tuneParam = sslLambdas, # current lambda associated with L2 penalty
                               stackParam = stackParam,
                               nnlsInd = nnlsInd,
                               stackTune = FALSE, # dont tune because we fix the parameter in each iteration
                               modelStandardize = standardize,
                               stackStandardize = standardize,
                               glmnetOpt = glmnetOpt,
                               xStandardize = xStandardize,
                               sampSzWeight = sampSzWeight,
                               weights = wgt)

            sigK <- warmRes$sigK
            sigStack <- warmRes$sigStack

            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))

                stackParam <- warmRes$mu
                tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
                set.seed(1)
                oecRidgeWS <- ridgeAltSpecWindow(data = data[indxList,], # current studies
                                            betaStart = warmRes$beta,
                                            wStart = warmRes$w,
                                            lambdaVec = tuneParam,
                                            mu = 0,#stackParam,
                                            Stackindx = indx2, # rows corresponding to test country
                                            tol = 0.01,
                                            nnlsInd = TRUE,
                                            low = 0,
                                            up = tune.grid[tune],
                                            objCriteria = TRUE,
                                            eta = oecEta,
                                            projs = 0,
                                            Avg = AvgW,
                                            wUpdate = "glment",
                                            standardize = TRUE,
                                            sampSzWeight = sampSzWeight,
                                            sigK = sigK,
                                            sigStack = sigStack,
                                            weights = wgt)

                predsVecOEC <- oecRidgePred(data = data[HOOList, -c(1,2)],
                                            mod = oecRidgeWS)

                resMat[study, tune] <- sqrt(mean( (predsVecOEC - data$Y[HOOList])^2  )) # stacking predictions

            }

        }

        resMean <- colMeans(resMat)
        etaStar <- tune.grid[ which.min(resMean)  ] # choose hyperparameters with best RMSE


    }


    return(best = etaStar)


}


