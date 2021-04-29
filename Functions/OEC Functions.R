
logisticNNLS <- function(dat, w_start = NULL, nnls = TRUE, stackInt = TRUE, 
                         psi = 0.1, simplexInd = TRUE, solverName = "ECOS", 
                         accelInd = FALSE, tau = 0.01, alpha = 0, itrMax = 10,
                         stackTune = 0, stepDec = 3){
    # tau is error tolerance for GD convergence
    # alpha is regularization parameter
    # accelInd is indicator of whether to use regular GD or ACD
    # itr max is max iterations for gradient descent
    # stack Int indicator of whetehr to use intercept
    # stepDec how many times to cut step size in half before stopping
    
    # first column of dat is Y
    stackPenal <- stackTune
    #print(paste("dim dat:", dim(dat)))
    if(is.null(solverName)){
        # if no solver specified then use Gradient Descent
        
        if(stackInt){
            #print("stack int")
            dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        }
        
        p <- ncol(dat) - 1 # dim including intercept vector of ones
        
        if(!is.null(w_start)){
            # starting value if none specified
            w_opt <- w_start # starts out just returning the starting value
            obj0 <- objValLogisticIndiv(X = dat[,-1], 
                                        Y = dat[,1], 
                                        beta = w_start, 
                                        stackTune = stackPenal) # objective value at starting point
            #print(paste("starting objective wNNLS", obj0))
            
        }else{
            
            w_opt <- w_start <- rep(1 / p, p)
            obj0 <- objValLogisticIndiv( X = dat[,-1], 
                                         Y = dat[,1], 
                                         beta = w_start, 
                                         stackTune = stackPenal) # objective value at starting point set arbitrarily small since we are maximizing and no w_start given
            #print("no starting value give so no starting objective")
            
        }
        
        w_r <- w_start # arbitrarily set to begin with; rth step
        w_prev <- w_start # r-1 th step is initialized as input
        gamma <- 0 # initialize at 0 so first step of AGD is just normal GD step
        obj1 <- obj0 # initialize for below to ensure loss function is improving
        
        # gradient
        if(stackInt){
            # if there is an intercept
            wGrad <- t( as.matrix(dat[,-1])) %*% 
                (as.vector(dat[,1]) - (1 / (1 + exp(as.matrix(-dat[,-1]) %*% 
                   as.vector(w_r) )  )  ) ) - 2 * stackTune * c(0, w_r[-1]) # stackTune is regularization parameter
        }else{
            # if there is no intercept
            wGrad <- t( as.matrix(dat[,-1])) %*% 
                (as.vector(dat[,1]) - (1 / (1 + exp(as.matrix(-dat[,-1]) %*% 
                     as.vector(w_r) )  )  ) ) - 2 * stackTune * w_r # stackTune is regularization parameter
        }
        
        # determine step size
        if(accelInd){
            # if accelerated gradient descent
            delta <- eta * norm(wGrad, "2") / norm(w_start, "2") # starting step size for accelerated grad descent is ratio of norms
        }else{
            # if not accelerated gradient descent
            etaSVD <- svd(t(as.matrix(dat[,-c(1,2)])) %*% as.matrix(dat[,-c(1,2)]) ) # take out Y and intercept of 1s
            L <- max( etaSVD$d ) / 4 + stackTune # max singular value + lambda (penalty hyperparameter)
            delta <- 1 / L
            
        }
        
        eps <- tau + 1
        itrs <- 1 # initialize step at 1
        
        # gradient descent
        while (abs(eps) > tau && itrs <= itrMax){
            
            deltaWork <- delta # working delta reverts back to delta between iterations but is updated within to ensure improving loss
            
            
            # gradient 
            
            if(stackInt){
                # if there is an intercept
                objGradient <- t(as.matrix(dat[, -1])) %*% 
                    ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                    2 * stackTune * c(0, w_r[-1])# stackTune is regularization parameter
            }else{
                # if there is no intercept
                objGradient <- t(as.matrix(dat[, -1])) %*% 
                    ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                    2 * stackTune * w_r # stackTune is regularization parameter
            }
            
            # accelerated  gradient descent step for jth study
            if(accelInd){
                # accelerated gradient descent
                gamma <- (j - 1) / (j - 1 + 3) * gamma + objGradient # accelerated grad descent
                
                w_r <- w_prev + delta * gamma # update for accelerated gradient descent
                
                if(nnls){
                    # projected GD
                    if(stackInt){
                        w_r <- c(w_r[1], w_r[-1] * I(w_r[-1] >= 0) ) # do not impose non-negativity on intercept
                    }else{
                        w_r <- w_r * I(w_r >= 0)
                        
                    }
                    
                    
                    
                }
                
                w_prev <- w_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                itrs <- itrs + 1 # move step forward
                
                #print(w_r)
                
                obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                             Y = dat[,1], 
                                             beta = w_r, 
                                             stackTune = stackPenal ) # current objective value
                eps <- obj1 - obj0
                #print(eps)
                obj0 <- obj1 # update so previous objective is last iteration's
                
            }else{ 
                # regular gradient descent
                  
                optFlag <- TRUE
                
                # ensure increasing gradient
                while(optFlag){
                    w_work <- w_prev + deltaWork * objGradient # update for gradient descent
                    
                    if(nnls){
                        # projected GD for NNLS
                        if(stackInt){
                            w_r <- c(w_r[1], w_r[-1] * I(w_r[-1] >= 0) ) # do not impose non-negativity on intercept
                        }else{
                            w_r <- w_r * I(w_r >= 0)
                            
                        }
                    }
                    
                    # check to see if loss improved
                    obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                                 Y = dat[,1], 
                                                 beta = w_work, 
                                                 stackTune = stackPenal) # current objective value
                    
                    if( round(obj1, 4) <= round(obj0, 4)  ){
                        #print("update deltawork")
                        # if no improvement cute step size in half
                        deltaWork <- deltaWork / 2
                        
                        if(deltaWork < delta * (1/2)^stepDec ){
                            #print("Break reach delta limit")
                            eps <- 0 # interrupt while loop
                            itrs <- itrMax + 1
                            optFlag <- FALSE
                        }
                        
                    }else{
                        
                        # if improvement then set w_r to working w
                        #print("improve")
                        w_r <- w_work
                        eps <- obj1 - obj0
                        #print(paste0("eps: ", eps))
                        obj0 <- obj1 # update so previous objective is last iteration's
                        optFlag <- FALSE # break while loop
                        
                        # update gradient
                        #####################
                        # added 2/25/20 while experimenting with armijos rule
                        #####################
                        if(stackInt){
                            # if there is an intercept
                            objGradient <- t(as.matrix(dat[, -1])) %*% 
                                ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                                2 * stackTune * c(0, w_r[-1])# stackTune is regularization parameter
                        }else{
                            # if there is no intercept
                            objGradient <- t(as.matrix(dat[, -1])) %*% 
                                ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                                2 * stackTune * w_r # stackTune is regularization parameter
                        }
                       ############################# 
                    }
                }
                
                #print(w_r)
                w_prev <- w_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                itrs <- itrs + 1 # move step forward
                
                #obj1 <- objValLogisticIndiv(dat[,-1], dat[,1], w_r) # current objective value
                
            }
        } # added here
        
        w_opt <- w_r # final is w_opt
         # print(paste("final obj:", obj1))
         # print(paste("final eps:", eps))
         # print(paste("final iterations:", itrs))
        if( abs(eps) <= tau ){
           # print("abs(eps) <= tau")
        }
        #print(paste0("itrs: ", itrs))
        
    }else if(solverName == "glmnet"){
        
        if(stackInt){
            #print("stack int")
            dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        }
        
        p <- ncol(dat) - 1 # dim including intercept vector of ones
        
        # use glmnet
        
        # if(stackInt){
        #     print("stack int")
        #     dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        # }
        
        if(is.null(w_start)){
            # starting value if none specified

            w_start <- rep(1 / p, p)
            
        }
        # print("dim beta")
        # print(dim(dat[,-1]))
        # print("length w")
        # print(length(w_start))
        # print(head(dat))
        
        obj0 <- objValLogisticIndiv(X = dat[,-1], 
                                    Y = dat[,1], 
                                    beta = w_start, 
                                    stackTune = stackPenal) # objective value at starting point
        #print(paste("starting objective", obj0))
        
        if(nnls){
            # with nonnegativity constraint
            if(stackInt){
                loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                family = "binomial", lower.limits = 0, alpha = 0,
                                lambda = stackPenal, 
                                standardize = TRUE,
                                intercept = TRUE)
                
                w_r <- as.vector(loReg$beta)
                w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                
            }else{
                loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-1]), 
                                family = "binomial", lower.limits = 0, alpha = 0,
                                lambda = stackPenal, 
                                standardize = TRUE,
                                intercept = FALSE)
                
                w_r <- as.vector(loReg$beta)
            }
            
            
        }else{
            # without nonnegativity constraint
            if(stackInt){
                loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                family = "binomial", alpha = 0,
                                lambda = stackPenal, 
                                standardize = TRUE,
                                intercept = TRUE)
                
                w_r <- as.vector(loReg$beta)
                w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                
            }else{
                loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-1]), 
                                family = "binomial", alpha = 0,
                                lambda = stackPenal, 
                                standardize = TRUE,
                                intercept = FALSE)
                
                w_r <- as.vector(loReg$beta)
            }
            
        }
        
        # print("dim beta")
        # print(dim(beta))
        # print("length w")
        # print(length(w_r))
        
       # w_opt <- as.vector(loReg$beta)
        obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                     Y = dat[,1], 
                                     beta = w_r, 
                                     stackTune = stackPenal ) # current objective value
        
        if(obj1 > obj0){
            # only use new w if it produces improvement in objective value
            w_opt <- w_r
            #print(paste0("glmnet w Opt improvement: ", obj1 - obj0))
        }else{
            w_opt <- w_start
        }
        
    }else{
        # solve with CVXR
        
        library(CVXR)
        if(stackInt){
            
            w <- Variable(ncol(dat) ) # minus 1 columns for Y plus 1 for intercept
            dat <- cbind( dat[, 1], 1, dat[, -1] )
        }else{
            w <- Variable(ncol(dat) - 1) # minus 1 columns for Y
        }
        
        
        if(nnls && simplexInd){
            # with nonnegativity constraint and no penalization
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w))
            prob <- Problem(Maximize(obj), 
                            constraints = list(w >= 0, sum(w) == 1)) # nnls constraint
        }else if(nnls && !simplexInd){
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w))
            prob <- Problem(Maximize(obj), 
                            constraints = list(w >= 0)) # nnls constraint
            
        }else{
            # without nonnegativity constraint but with penalization term on outside
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w)) - stackPenal * norm2(w[-1])^2
            prob <- Problem(Maximize(obj)) # nnls constraint
        }
        
        
        
        
        result <- solve(prob)
        w_opt <- result$getValue(w)
    }
    
    
    return(w_opt)
    
    
}

logisticNNLS2 <- function(dat, w_start = NULL, nnls = TRUE, stackInt = TRUE, 
                         psi = 0.1, simplexInd = TRUE, solverName = "ECOS", 
                         accelInd = FALSE, tau = 0.01, alpha = 0, itrMax = 10,
                         stackTune = 0, stepDec = 3){
    
    # same as other LogisticNNLS but uses GLM when less than 1 covariate 
    # tau is error tolerance for GD convergence
    # alpha is regularization parameter
    # accelInd is indicator of whether to use regular GD or ACD
    # itr max is max iterations for gradient descent
    # stack Int indicator of whetehr to use intercept
    # stepDec how many times to cut step size in half before stopping
    
    p <- ncol(dat) - 1
    
    # first column of dat is Y
    stackPenal <- stackTune
    #print(paste("dim dat:", dim(dat)))
    if(is.null(solverName)){
        # if no solver specified then use Gradient Descent
        
        if(stackInt){
            #print("stack int")
            dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        }
        
        p <- ncol(dat) - 1 # dim including intercept vector of ones
        
        if(!is.null(w_start)){
            # starting value if none specified
            w_opt <- w_start # starts out just returning the starting value
            obj0 <- objValLogisticIndiv(X = dat[,-1], 
                                        Y = dat[,1], 
                                        beta = w_start, 
                                        stackTune = stackPenal) # objective value at starting point
            #print(paste("starting objective wNNLS", obj0))
            
        }else{
            
            w_opt <- w_start <- rep(1 / p, p)
            obj0 <- objValLogisticIndiv( X = dat[,-1], 
                                         Y = dat[,1], 
                                         beta = w_start, 
                                         stackTune = stackPenal) # objective value at starting point set arbitrarily small since we are maximizing and no w_start given
            #print("no starting value give so no starting objective")
            
        }
        
        w_r <- w_start # arbitrarily set to begin with; rth step
        w_prev <- w_start # r-1 th step is initialized as input
        gamma <- 0 # initialize at 0 so first step of AGD is just normal GD step
        obj1 <- obj0 # initialize for below to ensure loss function is improving
        
        # gradient
        if(stackInt){
            # if there is an intercept
            wGrad <- t( as.matrix(dat[,-1])) %*% 
                (as.vector(dat[,1]) - (1 / (1 + exp(as.matrix(-dat[,-1]) %*% 
                                                        as.vector(w_r) )  )  ) ) - 2 * stackTune * c(0, w_r[-1]) # stackTune is regularization parameter
        }else{
            # if there is no intercept
            wGrad <- t( as.matrix(dat[,-1])) %*% 
                (as.vector(dat[,1]) - (1 / (1 + exp(as.matrix(-dat[,-1]) %*% 
                                                        as.vector(w_r) )  )  ) ) - 2 * stackTune * w_r # stackTune is regularization parameter
        }
        
        # determine step size
        if(accelInd){
            # if accelerated gradient descent
            delta <- eta * norm(wGrad, "2") / norm(w_start, "2") # starting step size for accelerated grad descent is ratio of norms
        }else{
            # if not accelerated gradient descent
            etaSVD <- svd(t(as.matrix(dat[,-c(1,2)])) %*% as.matrix(dat[,-c(1,2)]) ) # take out Y and intercept of 1s
            L <- max( etaSVD$d ) / 4 + stackTune # max singular value + lambda (penalty hyperparameter)
            delta <- 1 / L
            
        }
        
        eps <- tau + 1
        itrs <- 1 # initialize step at 1
        
        # gradient descent
        while (abs(eps) > tau && itrs <= itrMax){
            
            deltaWork <- delta # working delta reverts back to delta between iterations but is updated within to ensure improving loss
            
            
            # gradient 
            
            if(stackInt){
                # if there is an intercept
                objGradient <- t(as.matrix(dat[, -1])) %*% 
                    ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                    2 * stackTune * c(0, w_r[-1])# stackTune is regularization parameter
            }else{
                # if there is no intercept
                objGradient <- t(as.matrix(dat[, -1])) %*% 
                    ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                    2 * stackTune * w_r # stackTune is regularization parameter
            }
            
            # accelerated  gradient descent step for jth study
            if(accelInd){
                # accelerated gradient descent
                gamma <- (j - 1) / (j - 1 + 3) * gamma + objGradient # accelerated grad descent
                
                w_r <- w_prev + delta * gamma # update for accelerated gradient descent
                
                if(nnls){
                    # projected GD
                    w_r <- w_r * I(w_r >= 0)
                }
                
                w_prev <- w_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                itrs <- itrs + 1 # move step forward
                
                #print(w_r)
                
                obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                             Y = dat[,1], 
                                             beta = w_r, 
                                             stackTune = stackPenal ) # current objective value
                eps <- obj1 - obj0
                #print(eps)
                obj0 <- obj1 # update so previous objective is last iteration's
                
            }else{ 
                # regular gradient descent
                
                optFlag <- TRUE
                
                # ensure increasing gradient
                while(optFlag){
                    w_work <- w_prev + deltaWork * objGradient # update for gradient descent
                    
                    if(nnls){
                        # projected GD for NNLS
                        w_work <- w_work * I(w_work >= 0)
                    }
                    
                    # check to see if loss improved
                    obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                                 Y = dat[,1], 
                                                 beta = w_work, 
                                                 stackTune = stackPenal) # current objective value
                    
                    if( round(obj1, 4) <= round(obj0, 4)  ){
                        #print("update deltawork")
                        # if no improvement cute step size in half
                        deltaWork <- deltaWork / 2
                        
                        if(deltaWork < delta * (1/2)^stepDec ){
                            #print("Break reach delta limit")
                            eps <- 0 # interrupt while loop
                            itrs <- itrMax + 1
                            optFlag <- FALSE
                        }
                        
                    }else{
                        
                        # if improvement then set w_r to working w
                        #print("improve")
                        w_r <- w_work
                        eps <- obj1 - obj0
                        #print(paste0("eps: ", eps))
                        obj0 <- obj1 # update so previous objective is last iteration's
                        optFlag <- FALSE # break while loop
                        
                        # update gradient
                        #####################
                        # added 2/25/20 while experimenting with armijos rule
                        #####################
                        if(stackInt){
                            # if there is an intercept
                            objGradient <- t(as.matrix(dat[, -1])) %*% 
                                ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                                2 * stackTune * c(0, w_r[-1])# stackTune is regularization parameter
                        }else{
                            # if there is no intercept
                            objGradient <- t(as.matrix(dat[, -1])) %*% 
                                ( as.vector(dat[,1]) - 1 / (1 + exp(as.matrix(-dat[, -1]) %*% as.vector(w_r))) ) - 
                                2 * stackTune * w_r # stackTune is regularization parameter
                        }
                        ############################# 
                    }
                }
                
                #print(w_r)
                w_prev <- w_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                itrs <- itrs + 1 # move step forward
                
                #obj1 <- objValLogisticIndiv(dat[,-1], dat[,1], w_r) # current objective value
                
            }
        } # added here
        
        w_opt <- w_r # final is w_opt
        # print(paste("final obj:", obj1))
        # print(paste("final eps:", eps))
        # print(paste("final iterations:", itrs))
        if( abs(eps) <= tau ){
            # print("abs(eps) <= tau")
        }
        #print(paste0("itrs: ", itrs))
        
    }else if(solverName == "glmnet"){
        
        if(stackInt){
            #print("stack int")
            dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        }
        
        p <- ncol(dat) - 1 # dim including intercept vector of ones
        
        # use glmnet
        
        # if(stackInt){
        #     print("stack int")
        #     dat <- as.data.frame( cbind( dat[, 1], 1, dat[, -1] ) )
        # }
        
        if(is.null(w_start)){
            # starting value if none specified
            
            w_start <- rep(1 / p, p)
            
        }
        # print("dim beta")
        # print(dim(dat[,-1]))
        # print("length w")
        # print(length(w_start))
        # print(head(dat))
        
        obj0 <- objValLogisticIndiv(X = dat[,-1], 
                                    Y = dat[,1], 
                                    beta = w_start, 
                                    stackTune = stackPenal) # objective value at starting point
        #print(paste("starting objective", obj0))
        
        if(nnls){
            # with nonnegativity constraint
            if(stackInt){
                
                if(p > 1){
                    # if more than one covaraite
                    loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                    family = "binomial", lower.limits = 0, alpha = 0,
                                    lambda = stackPenal, 
                                    standardize = TRUE,
                                    intercept = TRUE)
                    
                    w_r <- as.vector(loReg$beta)
                    w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                }else{
                    # if one covariate
                    library(CVXR)
                    
                    beta <- Variable(2)
                    obj <- -sum(logistic(- cbind(1, dat[y <= 0, ]) %*% beta)) - 
                                    sum(logistic(cbind(1, dat[y == 1, ]) %*% beta))
                    constraints <- list(beta >= 0)
                    prob <- Problem(Maximize(obj), constraints)
                    result <- solve(prob)
                    w_r <- as.vector( prob$getValue(beta) )
                    
                    rm(prob)

                }
                
                
            }else{
                
                if(p > 1){
                    # if more than one covaraite
                    loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                    family = "binomial", lower.limits = 0, alpha = 0,
                                    lambda = stackPenal, 
                                    standardize = TRUE,
                                    intercept = FALSE)
                    
                    w_r <- as.vector(loReg$beta)
                    w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                }else{
                    # if one covariate
                    library(CVXR)
                    
                    beta <- Variable(1)
                    obj <- -sum(logistic(-dat[y <= 0, ] %*% beta)) - 
                        sum(logistic(dat[y == 1, ] %*% beta))
                    constraints <- list(beta >= 0)
                    prob <- Problem(Maximize(obj), constraints)
                    result <- solve(prob)
                    w_r <- as.vector( prob$getValue(beta) )
                    
                    rm(prob)
                    
                }
                

            }
            
            
        }else{
            # without nonnegativity constraint
            if(stackInt){
                if(p > 1){
                    # if more than one covaraite
                    loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                    family = "binomial", 
                                    alpha = 0,
                                    lambda = stackPenal, 
                                    standardize = TRUE,
                                    intercept = TRUE)
                    
                    w_r <- as.vector(loReg$beta)
                    w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                }else{
                    # if one covariate
                    colnames(dat) <- c("y", paste0("x_", 1:(ncol(dat) -1 ) ))
                    loReg <- glm(y ~., data = dat, 
                                 family = "binomial", model = FALSE)
                    
                    w_r <- as.vector(loReg$coefficients)
                    
                }
                
            }else{
                if(p > 1){
                    # if more than one covaraite
                    loReg <- glmnet(y = as.vector(dat[,1]), x = as.matrix(dat[,-c(1,2)]), # dat includes column of ones so remove that so do not penalize intercept coefficient 
                                    family = "binomial", 
                                    alpha = 0,
                                    lambda = stackPenal, 
                                    standardize = TRUE,
                                    intercept = FALSE)
                    
                    w_r <- as.vector(loReg$beta)
                    w_r <- c( as.vector( loReg$a0 ), w_r) # add intecept
                }else{
                    # if one covariate
                    colnames(dat) <- c("y", paste0("x_", 1:(ncol(dat) -1 ) ))
                    loReg <- glm(y ~. -1, data = dat, 
                                 family = "binomial", model = FALSE)
                    
                    w_r <- as.vector(loReg$coefficients)
                    
                }
            }
            
        }
        
        # print("dim beta")
        # print(dim(beta))
        # print("length w")
        # print(length(w_r))
        
        # w_opt <- as.vector(loReg$beta)
        obj1 <- objValLogisticIndiv( X = dat[,-1], 
                                     Y = dat[,1], 
                                     beta = w_r, 
                                     stackTune = stackPenal ) # current objective value
        
        if(obj1 > obj0){
            # only use new w if it produces improvement in objective value
            w_opt <- w_r
            #print(paste0("glmnet w Opt improvement: ", obj1 - obj0))
        }else{
            w_opt <- w_start
        }
        
    }else{
        # solve with CVXR
        
        library(CVXR)
        if(stackInt){
            
            w <- Variable(ncol(dat) ) # minus 1 columns for Y plus 1 for intercept
            dat <- cbind( dat[, 1], 1, dat[, -1] )
        }else{
            w <- Variable(ncol(dat) - 1) # minus 1 columns for Y
        }
        
        
        if(nnls && simplexInd){
            # with nonnegativity constraint and no penalization
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w))
            prob <- Problem(Maximize(obj), 
                            constraints = list(w >= 0, sum(w) == 1)) # nnls constraint
        }else if(nnls && !simplexInd){
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w))
            prob <- Problem(Maximize(obj), 
                            constraints = list(w >= 0)) # nnls constraint
            
        }else{
            # without nonnegativity constraint but with penalization term on outside
            obj <- -sum(dat[dat$Y <= 0, -1] %*% w) - sum(logistic(-dat[,-1] %*% w)) - stackPenal * norm2(w[-1])^2
            prob <- Problem(Maximize(obj)) # nnls constraint
        }
        
        
        
        
        result <- solve(prob)
        w_opt <- result$getValue(w)
    }
    
    
    return(w_opt)
    
    
}

logBetaOpt <- function(dat, 
                       beta_start, 
                       w, 
                       lambda, 
                       eta = 0.01, # was 0.001 2-29
                       itrs = 5, 
                       accelInd = FALSE, 
                       stackTune = 0, 
                       penVec = 0,
                       studyBeta_itrMax = 1000,
                       stepDec = 3){
    
    
    # penaltyVec is a K x 1 vector where each element is the penalization hyperparamerter for the corresponding L2 norm
    # accelInd indicates whether to do accelerated gradient descent
    # beta is K x p matrix
    # dat is dataframe with following columns: Study   Y   X_1   X_2  ....  X_p
    # lambda is tradeoff between stacking and individual loss functions
    # stepDec number of times to cut beta in half before stopping
    
    penaltyVec <- penVec
    stackPenalty <- stackTune
    K <- length(unique(dat$Study))
    indxList <- vector("list", length = K)
    
    beta_r <- beta_start # arbitrarily set to begin with; rth step
    beta_prev <- beta_start # r-1 th step is initialized as input
    
    deltaVec <- vector(length = K) # each element is the starting constant for accelerated gradient descent
    
    if(length(w) == K + 1){
        # if there is a stack intercept
        
        # break up into different elements for dimensions below
        
        w_int <- w[1]
        w <- w[-1] 
    }else{
        w_int <- 0
    }
    
    studyOrder <- sample.int(K, K) # randomize order
    for(j in studyOrder){
        
        indxList[[j]] <- which(dat$Study == j) # list of indices of row numbers for each study
        
        betaGrad <- lambda * w[j] * t( as.matrix(dat[,-c(1,2)])) %*% 
            (as.vector(dat$Y) - ( 1 / (1 + exp(-as.matrix(dat[,-c(1,2)]) %*% 
                                                 t(as.matrix(beta_prev)) %*% as.vector(w) + w_int ))) ) +  # was - before η
            (1 - lambda) * ( t( as.matrix(dat[ indxList[[j]], -c(1,2)])) %*%
            (as.vector(dat$Y[ indxList[[j]] ]) -
                ( 1 / (1 + exp(-as.matrix(dat[ indxList[[j]], -c(1,2)]) %*% 
                                  as.vector(beta_prev[j, ]))))  ) -
                           2 * penaltyVec[j] * beta_prev[j, ]       # regularization term
            )
        if(accelInd){
            # if accelerated gradient descent
            deltaVec[j] <- eta * norm(betaGrad, "2") / norm(beta_start[j, ], "2") # starting step size for accelerated grad descent is ratio of norms
        }else{
            # if not accelerated gradient descent
            # associated with individal study loss
            etaSVD <- svd(t(as.matrix(dat[dat$Study == j,-c(1,2,3)])) %*% 
                              as.matrix(dat[dat$Study == j,-c(1,2,3)]) ) # remove intercept of ones
            
            # associated with stacking loss
            stackSVD <- svd(t(as.matrix( w[j] * dat[,-c(1,2,3)])) %*% 
                                as.matrix( w[j] * dat[,-c(1,2,3)]) ) # remove intercept of ones
            
            # lipschitz constant
            L <-  max( etaSVD$d ) / 4 + max( stackSVD$d ) / 4   # max singular value
            L <- L + penaltyVec[j] + stackTune  # add in regularization terms
            
            # step size
            deltaVec[j] <- 1 / L
        }
    }
    
    
    len <- length(w)
    
    # this is (t(beta) of how it is in the rest of the code)
    
    # K x p matrix -- start with 0 so first step in accelerate grad descent is normal grad descent step
    gammaMat <- matrix(0, ncol = ncol(beta_start), nrow = nrow(beta_start)) # each row is the previous step from accelerate gradient descent for each study
    
    #print("accelerated gradient descent")
    it <- 1 # iterations
    
    obj0 <- obj0_k <- objLogisticEns(data = dat, beta = beta_start, wVal = w, 
                                     lambdaVal = lambda, penVec = penaltyVec, 
                                     stackInt = TRUE, 
                                     stackTune = stackPenalty ) # arbitrary starting objective value
    eps <- eta + 1 # arbitrarily set to be larger than error tolerance
    
    flag <- TRUE
    
    while(flag == TRUE){
        # iterate through a set number of descent updates
        
        #print(paste0("itrs: ", i))
        
        betaWork <- beta_r
        
        studyOrder <- sample.int(K, K)
        
        for (j in studyOrder){
            optFlag <- TRUE
            
            deltaWork <- deltaVec
            betaWork <- beta_r # just added 1-10-20 may need to take out --- pretty sure this should not be included. updates appropriately below
            
            # iterate through the updates for each study's vector of beta coefficients
            # beta is a K x p matrix
            
            # gradient of objective for j^th study
            objGradient <- lambda * w[j] * t(as.matrix(dat[, -c(1,2)])) %*% ( as.vector(dat$Y) - 
                                    1 / (1 + exp(-as.matrix(dat[, -c(1,2)]) %*% 
                                    t(as.matrix(beta_prev)) %*% as.vector(w) + w_int)) ) +  # was - before η
                                    (1 - lambda) * (t( as.matrix(dat[ indxList[[j]], -c(1,2)])) %*%
                                    (  as.vector(dat$Y[ indxList[[j]] ]) -
                                    1 / (1 + exp( -as.matrix(dat[ indxList[[j]], -c(1,2)]) %*% 
                                                      as.vector(beta_prev[j, ]) )   ) ) - 
                                    2 * penaltyVec[j] * beta_prev[j, ]  )     # regularization term
            
            # accelerated  gradient descent step for jth study
            if(accelInd){
                # accelerated gradient descent
                gammaMat[j, ] <- (j - 1) / (j - 1 + 3) * gammaMat[j, ] + objGradient # accelerated grad descent
                
                beta_r[j, ] <- beta_prev[j, ] + deltaVec[j] * gammaMat[j, ] # update for accelerated gradient descent
                
                beta_prev <- beta_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                
            }else{
                
                # regular gradient descent
                #gammaMat[j, ] <- (j - 1) / (j - 1 + 3) * gammaMat[j, ] + objGradient # accelerated grad descent
                itrCntr <- 1
                while(optFlag){
                    
                    betaWork[j,] <- beta_prev[j, ] + deltaWork[j] * objGradient # update for accelerated gradient descent
                    
                    obj1_k <- objLogisticEns(data = dat, beta = betaWork, wVal = w, 
                                             lambdaVal = lambda, penVec = penaltyVec, 
                                             stackInt = TRUE, 
                                             stackTune = stackPenalty)
                    
                    if( round(obj1_k, 3) <= round(obj0_k, 3)  ){
                        # if no improvement cut step size in half
                        
                        deltaWork[j] <- deltaWork[j] * 0.5
                        #print(deltaWork[j])
                        #print(paste0("deltawork: ", deltaWork))
                        
                        if(deltaWork[j] < ( deltaVec[j] * (1/2)^stepDec )  ){
                            # if no improvement in this coordinate (i.e.., study's beta) for a long time then 
                            # interrupt while loop
                            #eps <- 0 
                            #i <- itrs + 1
                            optFlag <- FALSE
                            
                        }
                        
                    }else{
                        # if there is improvement then update beta and objective value

                        obj0_k <- obj1_k
                        itrCntr <- itrCntr + 1
                        # update if there is an improvement
                        beta_r[j, ] <- betaWork[j,] # update for accelerated gradient descent
                        beta_prev <- beta_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                        
                        ##############################
                        # Added in on  2/25/20 when exploring armijos rule
                        ##############################
                        # update gradient
                        
                        # gradient of objective for j^th study
                        objGradient <- lambda * w[j] * t(as.matrix(dat[, -c(1,2)])) %*%
                                        ( as.vector(dat$Y) -  1 / (1 + exp(-as.matrix(dat[, -c(1,2)]) %*%
                                            t(as.matrix(beta_prev)) %*% as.vector(w) + w_int)) ) +  # was - before η
                            (1 - lambda) * (t( as.matrix(dat[ indxList[[j]], -c(1,2)])) %*%
                                                (  as.vector(dat$Y[ indxList[[j]] ]) -
                                                       1 / (1 + exp( -as.matrix(dat[ indxList[[j]], -c(1,2)]) %*%
                                                                         as.vector(beta_prev[j, ]) )   ) ) -
                                                2 * penaltyVec[j] * beta_prev[j, ]  )     # regularization term
                        #####################################
                        
                        if(itrCntr > studyBeta_itrMax){
                            # break if too many updates
                            optFlag <- FALSE
                            
                        }
                    }
                    
                }
            } 
            
        } # maybe move down to after this chunk of code
        
        obj1 <- objLogisticEns(data = dat, beta = beta_r, 
                               wVal = w, 
                               lambdaVal = lambda, 
                               penVec = penaltyVec,
                               stackInt = TRUE, stackTune = stackPenalty)
        
        if(round(obj0, 4) < round(obj1, 4) ){
            # if improvement in objective
            eps <- obj1 - obj0 # obj difference
            obj0 <- obj1 # update previous step's obj with this step's obj
        }else{
            #print("no improvement")
            # otherwise set difference in objectives between subsequent steps to 0
            # do not update obj0
            eps <- 0
        }
        
        
        it <- it + 1 # move step forward
        #print(paste0("itrs limit:", itrs))
        #print(paste0("i: ", i))
        
        if(it > itrs){
            # if limit of iterations is reached stop while loop
            #print("i > itrs")
            flag <- FALSE
            
        } 
        
        if(abs(eps) <= eta){
            # if error is small enough stop while loop
            #print("abs(eps) <= eta")
            flag <- FALSE
        }
        
    }
    
    return(beta_prev)
    
}

logBetaOptA <- function(dat, 
                       beta_start, 
                       w, 
                       lambda, 
                       eta = 0.001, 
                       itrs = 5, 
                       accelInd = FALSE, 
                       stackTune = 0, 
                       penVec = 0,
                       studyBeta_itrMax = 1000){
    
    # uses Armijo's rule
    # penaltyVec is a K x 1 vector where each element is the penalization hyperparamerter for the corresponding L2 norm
    # accelInd indicates whether to do accelerated gradient descent
    # beta is K x p matrix
    # dat is dataframe with following columns: Study   Y   X_1   X_2  ....  X_p
    # lambda is tradeoff between stacking and individual loss functions
    
    penaltyVec <- penVec
    stackPenalty <- stackTune
    K <- length(unique(dat$Study))
    indxList <- vector("list", length = K)
    
    beta_r <- beta_start # arbitrarily set to begin with; rth step
    beta_prev <- beta_start # r-1 th step is initialized as input
    
    deltaVec <- vector(length = K) # each element is the starting constant for accelerated gradient descent
    
    if(length(w) == K + 1){
        # if there is a stack intercept
        
        # break up into different elements for dimensions below
        
        w_int <- w[1]
        w <- w[-1] 
    }else{
        w_int <- 0
    }
    
    studyOrder <- sample.int(K, K) # randomize order
    for(j in studyOrder){
        
        indxList[[j]] <- which(dat$Study == j) # list of indices of row numbers for each study
        
        betaGrad <- lambda * w[j] * t( as.matrix(dat[,-c(1,2)])) %*% 
            (as.vector(dat$Y) - ( 1 / (1 + exp(-as.matrix(dat[,-c(1,2)]) %*% 
                                                   t(as.matrix(beta_prev)) %*% as.vector(w) + w_int ))) ) +  # was - before η
            (1 - lambda) * ( t( as.matrix(dat[ indxList[[j]], -c(1,2)])) %*%
                                 (as.vector(dat$Y[ indxList[[j]] ]) -
                                      ( 1 / (1 + exp(-as.matrix(dat[ indxList[[j]], -c(1,2)]) %*% 
                                                         as.vector(beta_prev[j, ]))))  ) -
                                 2 * penaltyVec[j] * beta_prev[j, ]       # regularization term
            )
        if(accelInd){
            # if accelerated gradient descent
            deltaVec[j] <- eta * norm(betaGrad, "2") / norm(beta_start[j, ], "2") # starting step size for accelerated grad descent is ratio of norms
        }else{
            # if not accelerated gradient descent
            etaSVD <- svd(t(as.matrix(dat[dat$Study == j,-c(1,2,3)])) %*% as.matrix(dat[dat$Study == j,-c(1,2,3)]) ) # remove intercept of ones
            deltaVec[j] <- 4 / max( etaSVD$d ) # max singular value
        }
    }
    
    
    len <- length(w)
    
    # this is (t(beta) of how it is in the rest of the code)
    
    # K x p matrix -- start with 0 so first step in accelerate grad descent is normal grad descent step
    gammaMat <- matrix(0, ncol = ncol(beta_start), nrow = nrow(beta_start)) # each row is the previous step from accelerate gradient descent for each study
    
    #print("accelerated gradient descent")
    it <- 1 # iterations
    
    obj0 <- obj0_k <- objLogisticEns(data = dat, beta = beta_start, wVal = w, 
                                     lambdaVal = lambda, penVec = penaltyVec, 
                                     stackInt = TRUE, 
                                     stackTune = stackPenalty ) # arbitrary starting objective value
    eps <- eta + 1 # arbitrarily set to be larger than error tolerance
    
    flag <- TRUE
    
    while(flag == TRUE){
        # iterate through a set number of descent updates
        
        #print(paste0("itrs: ", i))
        
        betaWork <- beta_r
        
        studyOrder <- sample.int(K, K)
        
        for (j in studyOrder){
            optFlag <- TRUE
            
            deltaWork <- deltaVec
            betaWork <- beta_r # just added 1-10-20 may need to take out --- pretty sure this should not be included. updates appropriately below
            
            # iterate through the updates for each study's vector of beta coefficients
            # beta is a K x p matrix
            
            # gradient of objective for j^th study
            # armijo's rule
            objGradientFn <- function(beta_prev){
                lambda * w[j] * t(as.matrix(dat[, -c(1,2)])) %*% ( as.vector(dat$Y) - 
                                    1 / (1 + exp(-as.matrix(dat[, -c(1,2)]) %*% 
                                          t(as.matrix(beta_prev)) %*% as.vector(w) + w_int)) ) +  # was - before η
                (1 - lambda) * (t( as.matrix(dat[ indxList[[j]], -c(1,2)])) %*%
                            (  as.vector(dat$Y[ indxList[[j]] ]) -
                                    1 / (1 + exp( -as.matrix(dat[ indxList[[j]], -c(1,2)]) %*% 
                                        as.vector(beta_prev[j, ]) )   ) ) - 
                                             2 * penaltyVec[j] * beta_prev[j, ]  )     # regularization term
                }
            # accelerated  gradient descent step for jth study
            if(accelInd){
                # accelerated gradient descent
                gammaMat[j, ] <- (j - 1) / (j - 1 + 3) * gammaMat[j, ] + objGradient # accelerated grad descent
                
                beta_r[j, ] <- beta_prev[j, ] + deltaVec[j] * gammaMat[j, ] # update for accelerated gradient descent
                
                beta_prev <- beta_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                
            }else{
                
                # regular gradient descent
                #gammaMat[j, ] <- (j - 1) / (j - 1 + 3) * gammaMat[j, ] + objGradient # accelerated grad descent
                itrCntr <- 1
                while(optFlag){
                    
                    betaWork[j,] <- armijo(f, x, sigma=0.1, maximise = TRUE)     # beta_prev[j, ] + deltaWork[j] * objGradient # update for accelerated gradient descent
                    
                    obj1_k <- objLogisticEns(data = dat, beta = betaWork, wVal = w, 
                                             lambdaVal = lambda, penVec = penaltyVec, 
                                             stackInt = TRUE, 
                                             stackTune = stackPenalty)
                    
                    if( round(obj1_k, 3) <= round(obj0_k, 3)  ){
                        # if no improvement cut step size in half
                        
                        deltaWork[j] <- deltaWork[j] * 0.5
                        #print(deltaWork[j])
                        #print(paste0("deltawork: ", deltaWork))
                        
                        if(deltaWork[j] < ( deltaVec[j] * (1/2)^stepDec )  ){
                            # if no improvement in this coordinate (i.e.., study's beta) for a long time then 
                            # interrupt while loop
                            #eps <- 0 
                            #i <- itrs + 1
                            optFlag <- FALSE
                            
                        }
                        
                    }else{
                        # if there is improvement then update beta and objective value
                        
                        obj0_k <- obj1_k
                        itrCntr <- itrCntr + 1
                        # update if there is an improvement
                        beta_r[j, ] <- betaWork[j,] # update for accelerated gradient descent
                        beta_prev <- beta_r # update previous beta to be current beta   # should this be here or above inside previous for loop
                        
                        if(itrCntr > studyBeta_itrMax){
                            # break if too many updates
                            optFlag <- FALSE
                            
                        }
                    }
                    
                }
            } 
            
        } # maybe move down to after this chunk of code
        
        obj1 <- objLogisticEns(data = dat, beta = beta_r, 
                               wVal = w, 
                               lambdaVal = lambda, 
                               penVec = penaltyVec,
                               stackInt = TRUE, stackTune = stackPenalty)
        
        if(round(obj0, 4) < round(obj1, 4) ){
            # if improvement in objective
            eps <- obj1 - obj0 # obj difference
            obj0 <- obj1 # update previous step's obj with this step's obj
        }else{
            #print("no improvement")
            # otherwise set difference in objectives between subsequent steps to 0
            # do not update obj0
            eps <- 0
        }
        
        
        it <- it + 1 # move step forward
        #print(paste0("itrs limit:", itrs))
        #print(paste0("i: ", i))
        
        if(it > itrs){
            # if limit of iterations is reached stop while loop
            #print("i > itrs")
            flag <- FALSE
            
        } 
        
        if(abs(eps) <= eta){
            # if error is small enough stop while loop
            #print("abs(eps) <= eta")
            flag <- FALSE
        }
        
    }
    
    return(beta_prev)
    
}




logisticIndivGrad <- function(w_r, dat, alpha = 0){
    return(
        t( as.matrix(dat[,-1])) %*% 
            (as.vector(dat[,1]) - (1 / (1 + exp(as.matrix(-dat[,-1]) %*% 
                                                    as.vector(w_r) )  )  ) ) - 2 * alpha * c(0, w_r[-1])
    )
}



logisticEnsGrad <- function(beta_prev, dat, w, rowIndx, studyNum, stackTune = 0, penVec = 0){
    
    return(
        lambda * w[studyNum] * t(as.matrix(dat[, -c(1,2)])) %*% ( as.vector(dat$Y) - 1 / (1 + exp(-as.matrix(dat[, -c(1,2)]) %*% 
                                                                                                      t(as.matrix(beta_prev)) %*% as.vector(w[-1]) + w[1])) ) +  # was - before η
            (1 - lambda) * t( as.matrix(dat[ rowIndx, -c(1,2)])) %*%
            (  as.vector(dat$Y[ rowIndx ]) -
                   1 / (1 + exp( -as.matrix(dat[ rowIndx, -c(1,2)]) %*% as.vector(beta_prev) ))) - 2 * c(0, stackTune * w[-1] )
        
    )
}


EnsOptLogistic <- function(dat, beta_start, w_start, alpha_vec = NULL, lambda = 0.5, tau = 0.001, itrMax = 100, betaItrs = 1000, wItrs = 1000,
                           n_k = 1000, stackPenal = 0.1, optW = TRUE, stackInt = TRUE, modelInt = TRUE, itersbeta = 3, itersW = 3,
                           eta_beta = 0.001, eta_w = 2^(-6), nnlsInd = TRUE, stackPsi = 0.25, solverName = "glmnet",
                           stackTune = 0, penaltyVec = 0, avgInd = FALSE){
    
    K <- length(unique(dat$Study))
    # beta_start <- as.matrix(beta_start)
    # dat <- as.matrix(dat)
    # w_start <- as.vector(w_start)
    
    #library(CVXR)
    #library(rje)
    # itersbeta and itersW are the number of gradient descent updates to do for each alterating optimization
    # η_beta is the step size for the gradient descent update w.r.t. beta
    # η_w is the step size for the gradient descent update w.r.t. w (stacking weights)
    # alpha_vec is K x 1 vector with regularization parameters for K studies,
    # ϵ is error
    # τ is error tolerance
    if(modelInt == TRUE){
        dat <- cbind(dat[,1:2], 1, dat[,-c(1,2)]) # add column of ones for intercept to design matrix
    }
    
    if(avgInd){
        
        if(length(w_start) == K){
            w_start <- rep(1/ K, K)
            
        }else{
            w_start <- c(0, rep(1/ K, K ) )
        }
        
        w_work <- w_opt <- w_start
        
    }
    
    # number of observations
    n = nrow(dat)
    # number of estimators
    d = ncol(dat) - 2 # first two columns are Study and Y
    #n_k = 5000; # would need to alter to generalize
    itrs = 1 # iteration starting point
    # while loop
    eps = 1 + tau # initial error term to ensure above tolerance
    stackPenalty <- stackTune
    
    beta_opt <- beta_start
    
    objVec1 <- vector(length = itrMax + 1) # vector of objective values at each step
    objVec2 <- vector(length = itrMax + 1) # vector of objective values at each step
    
    objVec1[1] <- objLogisticEns(data = dat, 
                                 beta = beta_start, 
                                 wVal = w_start, 
                                 lambdaVal = lambda, 
                                 penVec = penaltyVec,
                                 stackInt = TRUE,
                                 stackTune = stackPenalty) # first one is objective from standard stacking
    
    objVec2[1] <- objVec1[1]
    obj0 <- objVec1[1] # starting objective value
    
    #print(paste0("starting objective: ", obj0))
    
    # if(stackInt == TRUE){
    #     # include intercept in stacking model
    #     w_opt = rep(0, K + 1) # add column of ones for intercept to design matrix
    # }else{
    #     w_opt = rep(0, K)
    # }
    # 
    w_opt <- w_start
    
    # etaSVD <- svd(t(as.matrix(dat[,-c(1,2)])) %*% as.matrix(dat[,-c(1,2)]) )
    # eta_beta <- max( etaSVD$d ) # max singular value
   
    cntr <- 0 # counter for accepted descents
    
    while (abs(eps) > tau && itrs <= itrMax){
        
        # with fixed w do inner minimization
        # print("  inneroptLasso   ")
        
        beta_work <- logBetaOpt(dat, beta_opt, 
                                w_opt, lambda, eta_beta, itrs = betaItrs, 
                                stackTune = stackPenalty, penVec = penaltyVec)
        #print(beta_opt)
        # print("  inneroptLasso complete  ")
        # print("  stackFit   ")
        
        obj1 <- objLogisticEns(data = dat, beta = beta_work, 
                               wVal = w_opt, 
                               lambdaVal = lambda, 
                               penVec = penaltyVec,
                               stackInt = TRUE, 
                               stackTune = stackPenalty) # objective value
        
        objVec1[itrs + 1] <- obj1
        
        if(obj1 > obj0){
            # only update beta if it is an improvement
            beta_opt <- beta_work
            
            obj0 <- obj1 # just added in 1/10/20
            
            #print(paste0("Obj val after Beta opt: ", obj1))
            
        }else{
            
            #print(paste0("Decrement in betaOpt, obj: ", obj1))
            
        }
        
        
        
        
        #w_opt <- logisticWOpt(dat, beta_opt, w_start, lambda)
        predsMat <- as.data.frame(as.matrix( dat[,-c(1,2)]) %*% t(as.matrix(beta_opt))) # predicted values
        predsMat <- data.frame(cbind(Y = dat$Y, predsMat))
        #colnames(predsMat)[1] <- "Y"
        
        # use normal logistic regression
        #mod <- glm(Y ~., data = predsMat, family = "binomial")
        #w_opt <- as.vector(mod$coefficients)
        # nnls logistic regression
        
        if(!avgInd){
            
            w_work <- as.vector( (logisticNNLS(predsMat, w_start = w_opt, stackInt = TRUE, nnls = nnlsInd,
                                               psi = stackPsi, simplexInd = FALSE, solverName = solverName, 
                                               itrMax = wItrs,
                                               stackTune = stackPenalty) ) )
            
        }

        
        #print(w_opt)

        obj <- objLogisticEns(data = dat, beta = beta_opt, 
                              wVal = w_work, 
                              lambdaVal = lambda, 
                              penVec = penaltyVec,
                              stackInt = TRUE,
                              stackTune = stackPenalty ) # objective value
        
        objVec2[itrs + 1] <- obj
        
        # logisticWOpt(X, Y, K, beta_start, w, λ, η = 2^(-6), itrs = 5)
        # w_opt = stackFit(X, Y, K, beta_opt, 1, stackPenal, optW);
        
        if(obj > obj0){
            # only update beta if it is an improvement
            w_opt <- w_work
            #print(paste0("Obj val after w opt: ", obj))
            
            ###### moved this here from below if else
            eps <- obj - obj0 # objective difference
            # print("objective difference")
            # print(eps)
            obj0 <- obj # current objective
            
            
        }else{
            eps <- 0 # added this in to prevent negatives
           # print(paste0("Decrement in wOpt, obj: ", obj))
            
        }
        
        # print("  stackFit complete  ")
        
        # print(paste0("Obj val: ", obj))
        
        

        # w_start <- w_opt # update value
        # beta_start <- beta_opt # update value
        itrs <- itrs + 1
        # print(paste0("Iteration: ", itrs))
    }
    
    # beta_opt <- logBetaOpt(dat, beta_start, w_start, lambda, eta_beta , 
    #                        itrs = betaItrs * 3)

    obj <- objLogisticEns(data = dat, 
                          beta = beta_opt, 
                          wVal = w_opt, 
                          lambdaVal = lambda, 
                          penVec = penaltyVec,
                          stackInt = TRUE,
                          stackTune = stackPenalty) # objective value
    #print(paste0("Final Obj: ", obj))
    
    # print(itrs)
    
    if(itrs >= itrMax){
        # print("Maximum Iterations Reached")
    }
    
    # remove place holder values
    objVec1 <- objVec1[objVec1 != 0]
    objVec2 <- objVec2[objVec2 != 0]
    
    return ( list(beta_opt, w_opt, itrs, objVec1, objVec2, obj) );
}

expit <- function(X){
    return(
        
        1 / (1 + exp(-X) )
    )
}


# 
# objLogisticEns <- function(data, beta, wVal, lambdaVal, penVec, stackInt = TRUE, stackTune = 0){
#     # dat-- columns are Study  Y    X_1   X_2 ..... X_p   where first column of X is ones if theres an intercept in each model
#     K <- length(unique(data$Study))
#     obj <- 0
#     
#     # beta <- as.matrix(beta)
#     # dat <- as.matrix(dat)
#     # w <- as.vector(w)
#     
#     if(K == length(wVal)){
#         # print("1")
#         # if no intercept in stacking regression
#         obj <- lambdaVal * ( -sum(as.matrix( data[data$Y <= 0, -c(1,2)]) %*% as.matrix(t(beta)) %*% as.vector(wVal)) - 
#                               sum(expit(-as.matrix( data[,-c(1,2)]) %*% as.matrix(t(beta)) %*% as.vector(wVal))) )  # w has intercept
#     }else if(length(wVal) == K + 1){
#         # print("2")
#         # if there is an intercept
#         obj <- lambdaVal * ( -sum(as.matrix( data[data$Y <= 0, -c(1,2)] )%*% as.matrix(t(beta)) %*% as.vector(wVal[-1]) + as.vector(wVal[1]) ) - 
#                               sum(expit(as.matrix( -data[,-c(1,2)]) %*% as.matrix(t(beta)) %*% as.vector(wVal[-1]) + as.vector(wVal[1])) ))  # w has intercept
#     }else{
#         # print("objective value")
#     }
#     
#     for(j in 1:K){
#         #print("3")
#         # iterate through each study to update objective with study-specific loss (each row of beta matrix is a vector of 
#         # study-specific regression coefficients)
#         studyIndx <- which(data$Study == j) # rows in this study
#         Y_zeros <- which(data$Y == 0) # rows with outcome == 0
#         indx0 <- intersect(Y_zeros, studyIndx) # observations in current study who had 0 outcome
#         
#         if(length(indx0) > 0 ){
#             obj <- obj + (1 - lambdaVal) * ( -sum( as.matrix( data[indx0, -c(1,2)]) %*% beta[j, ] ) - 
#                                               sum(expit( -as.matrix(data[studyIndx,-c(1,2)]) %*% beta[j, ]) ) )
#         }else if(length(indx0) == 0 ){
#             obj <- obj + (1 - lambdaVal) * (-sum(expit(-as.matrix(data[studyIndx,-c(1,2)]) %*% beta[j, ]) ) )
#         }
#     }
#     
#     obj <- obj - stackTune * sum(wVal[-1]^2) - sum( t(as.matrix(beta[,-1])^2) %*% as.vector(penVec) )  # penalization for stacking and SSLs
#     
#     return(obj) # return objective value
# } # original version --- wrong
# 


objLogisticEns <- function(data, beta, wVal, lambdaVal, penVec, stackInt = TRUE, stackTune = 0){
    # dat-- columns are Study  Y    X_1   X_2 ..... X_p   where first column of X is ones if theres an intercept in each model
    K <- length(unique(data$Study))
    obj <- 0
    
    # beta <- as.matrix(beta)
    # dat <- as.matrix(dat)
    # w <- as.vector(w)
    
    if(K == length(wVal)){
        #print("1")
        # if no intercept in stacking regression
        obj <- lambdaVal * objValLogisticIndiv( X = as.matrix( data[, -c(1,2)] )%*% as.matrix(t(beta)), 
                                                Y = data$Y, 
                                                beta = wVal, 
                                                stackTune = 0)
        
    }else if(length(wVal) == K + 1){
        #print("2")
        # if there is an intercept
        obj <- lambdaVal * objValLogisticIndiv( X = cbind(1,   as.matrix( data[, -c(1,2)] )%*% as.matrix(t(beta))    ), 
                                                Y = data$Y, 
                                                beta = wVal, 
                                                stackTune = 0) 
    }else{
        # print("objective value")
    }
    
    for(j in 1:K){
        #print("3")
        # iterate through each study to update objective with study-specific loss (each row of beta matrix is a vector of 
        # study-specific regression coefficients)
        studyIndx <- which(data$Study == j) # rows in this study
        Y_zeros <- which(data$Y == 0) # rows with outcome == 0
        indx0 <- intersect(Y_zeros, studyIndx) # observations in current study who had 0 outcome
        
        obj <- obj + (1 - lambdaVal) * objValLogisticIndiv( X = as.matrix( data[studyIndx, -c(1,2)] ), 
                                                            Y = data$Y[studyIndx], 
                                                            beta = beta[j, ], 
                                                            stackTune = 0)
    }
    
    obj <- obj - stackTune * sum(wVal[-1]^2) - sum( t(as.matrix(beta[,-1])^2) %*% as.vector(penVec) )  # penalization for stacking and SSLs
    
    return(obj) # return objective value
}


objValLogisticIndiv <- function(X, Y, beta, stackTune = 0){
    X <- as.matrix(X)
    Y <- as.vector(Y)
    beta <- as.vector(beta)
    
    # print(dim(X))
    # print(length(Y))
    # print(length(beta))
    
    indx0 <- which(Y == 0) # which outcomes are 0
    #print(paste0("indx: ", length(indx0)))
    cost0 <- log(1 - (exp(X[indx0,] %*% beta) / 
                          (1 + exp(X[indx0,] %*% beta))) )
    
    cost0 <- sum( cost0[is.finite(cost0)] ) # sum those that are finite
    
    cost1 <-  log( exp(X[-indx0,] %*% beta) / (1 + exp(X[-indx0,] %*% beta)))  
    cost1 <- sum( cost1[is.finite(cost1)] ) # sum those that are finite
    
    return( cost1 + cost0 - stackTune * sum(beta[-1]^2) )
} 


##############
# Warm Start
##############
warmStartLogistic <- function(  dat,
                                alpha_vec = sort( c(0, 0.0001, 0.001, 0.01, 
                                                    seq(0.001,1, length = 50), 
                                                            2:10) ),
                                alphaTune = TRUE,
                                alpha = NULL,
                                ridge = FALSE,
                                stackPenal = 0.1,
                                stackInt = TRUE,
                                nnlsInd = TRUE,
                                stackPsi = 0.25,
                                GD = TRUE,
                                simpInd = FALSE,
                                solverName = "glmnet", 
                                lambda = 0.5,
                                stackTune = 0,
                                tuningInd = TRUE){
    
    # GD is indicator of whether to use gradient descent
    # simpInd is a simplex indicator to restrict stacking coefs to sum to 1
    # stack tune is regularization parameter
    # tuningInd == TRUE indicates to tune stacking hyperparameter
    library(glmnet)
    library(ModelMetrics)
    # add intercept term 
    
    dat <- cbind(dat[,c(1,2)], 1, dat[,-c(1,2)]) # add intercept
    colnames(dat)[3] <- "V_0"
    
    trainInd <- TRUE # begins at TRUE-- says to generate matrix of SSL coefficients below
    
    K <- length(unique(dat$Study)) # number of studies
    alphaOpt_SSL <- vector(length = K) # vector of optimized hyperparameters
    
    # print(head(dat))
    
    ######################
    # Stacking with NNLS
    ######################
    
    # ens_betaMat_Int
    # 
    # # global intFlag = int; #starts this way and is updated if autoInt
    # global w_nnls
    # 
    # global w_nnlsInt
    
    n_k <- nrow(dat)
    
    ens_betaMat <- matrix(nrow = K, ncol = ncol(dat) - 2)
    
    errorMat <- matrix(nrow = K, ncol = length(alpha_vec))
    # store error terms
    
    ########################
    # Tune hyperparameter
    ########################
    if (alphaTune == TRUE){
        # print("tune hyperparameter")
        
        trainInd <- FALSE # reset to FALSE so no need to retrain SSLs below
        p <- ncol(dat) - 2 # number of covariates including intercept
        betaMat <- matrix(ncol = p, nrow = length(alpha_vec)) # store coefficients in matrix for each hyperparameter so we can save later
        aucVec <- vector(length = length(alpha_vec) ) # store errors associated with each hyperparameter (tested on other studies)
        
        # set held out
        for(studyNum in 1:K){
            
            # iterate through studies to hold out
            rowIndx <- which(dat$Study == studyNum)
            
            for(lam in 1:length(alpha_vec)){
                # alpha vec is the regularization hyperparameter for Ridge penalty
                tuneParam <- alpha_vec[lam]
                mod <- glmnet(y = as.vector(dat$Y[rowIndx]),
                                x = as.matrix( dat[rowIndx, -c(1,2, 3)] ), # remove intecept term as well 
                                family = "binomial", 
                                alpha = 0,          # no lasso penalty
                                lambda = tuneParam, # ridge penalty
                                standardize = TRUE,
                                intercept = TRUE) 
                
                betaVec <- as.vector( mod$beta )
                betaMat[lam, ] <- c( as.vector( mod$a0 ), betaVec)  # add intercept & save to matrix
                
                preds <- as.matrix( dat[-rowIndx, -c(1,2)] ) %*% as.vector( betaMat[lam, ] ) # predicted linear predictor (systematic component)
                probs <- 1 / ( 1 + exp( -preds ) ) # predicted probabilities  (i.e., expit(preds) )
                aucVec[lam] <- auc( dat$Y[-rowIndx], probs ) # auc
                
            }
            
            alphaOpt_SSL[studyNum] <- alpha_vec[ which.max( aucVec ) ] # choose L2 penalty associated with best AUC
            ens_betaMat[studyNum, ] <- betaMat[which.max( aucVec ), ] # save coefficients associated with best AUC-producing penalty term
            
        }
        
        # print(alphaOpt_SSL)
        ################################## hyper parameter tuning
        
    }else{
        # otherwise use prespecified alpha
        alpha_opt <- rep(alpha, K)
        # print("use preexisting alpha")
    }
    
    ########################
    # Train SSE classifiers
    ########################
    # train individual classifiers
    # print("train SSE models")
    # timeElaps =@elapsed begin
    # trainInd == TRUE means to train, if it is false then training happened above
    if(trainInd){
        
        for(studyNum in 1:K){
            
            # iterate through studies to train classifiers
            
            trainIndx <- which(dat$Study == studyNum)  # indices of study rows -- must be adjusted if uneven rows
            
            if(!GD){
                
                mod <- glm(Y ~. - 1, data = dat[trainIndx, -1], family = "binomial" ) # train logistic regression-- # intercept included in design matrix
                ens_betaMat[studyNum, ] <- as.vector( mod$coefficients ) # save to matrix
                rm(mod)
                
            }else{
                
                ens_betaMat[studyNum, ] <-  logisticNNLS(dat[trainIndx, -1], w_start = NULL, nnls = FALSE, stackInt = FALSE, 
                                                         psi = 0.1, simplexInd = FALSE, solverName = NULL, 
                                                         accelInd = FALSE, tau = 0.000001, alpha = 0, itrMax = 10000)
            }
            
            
            
            # ens_betaMat[studyNum, ] <-  mod$coefficients # save to matrix
            #
        }
    }
    
    
    
    ########################
    # Stacking - Traditional
    ########################
    # print("stacking")
    
    # if an intercept is included in the stacking regression
    if(stackInt == TRUE){
        
        if(tuningInd){
            
            # tune hyperparameter
            p <- ncol(dat) - 2 # number of covariates including intercept
            aucMat <- matrix(ncol = length(alpha_vec), nrow = K ) # store errors associated with each hyperparameter (tested on other studies)
            
            # set held out
            for(studyNum in 1:K){
                
                # iterate through studies to hold out
                rowIndx <- which(dat$Study == studyNum)
                
                # remove current study and corresponding coefficients
                predsMat <- as.matrix( dat[rowIndx, -c(1,2)] ) %*% 
                                    t( as.matrix(ens_betaMat[-studyNum,] ) )
 
                
                predsMat <- cbind(Y = dat$Y[rowIndx], predsMat)

                #colnames(predsMat)[1] <- "Y"
                
                for(lam in 1:length(alpha_vec)){
                    # alpha vec is the regularization hyperparameter for Ridge penalty
                    tuneParam <- alpha_vec[lam]
                    
                    if(nnlsInd){
                        mod <- glmnet(y = as.vector(predsMat[, 1]), 
                                      x = as.matrix(predsMat[, -1]), 
                                      family = "binomial", 
                                      alpha = 0,          # no lasso penalty
                                      lambda = tuneParam, # ridge penalty
                                      standardize = TRUE,
                                      intercept = TRUE,
                                      lower.limits = 0) 
                    }else{
                        mod <- glmnet(y = as.vector(predsMat[, 1]), 
                                      x = as.matrix(predsMat[, -1]), 
                                      family = "binomial", 
                                      alpha = 0,          # no lasso penalty
                                      lambda = tuneParam, # ridge penalty
                                      standardize = TRUE,
                                      intercept = TRUE) 
                    }

                    w <- as.vector( mod$beta ) # stacking weights
                    w <- c( as.vector( mod$a0 ), w ) # add intecept
                    preds <- as.matrix( cbind(1, predsMat[, -1]) ) %*% w # predicted linear predictor (systematic component)
                    probs <- 1 / ( 1 + exp( -preds ) ) # predicted probabilities  (i.e., expit(preds) )
                    aucMat[studyNum, lam] <- auc( as.vector(predsMat[, 1]), probs ) # save AUC for this HO study and this tuning param
                }
    
            }
            
            optIndx <- which.max( colMeans( aucMat ) ) # average across held out studies
            stackTune <- alpha_vec[optIndx]
            alpha_opt <- stackTune
        }
            # fit final stack with optimized parameter
            predsMat <- as.data.frame( 
                            as.matrix( dat[, -c(1,2)] ) %*% 
                            t( as.matrix(ens_betaMat ) )
                        )
            predsMat <- cbind( Y = dat$Y, predsMat)
    
            #colnames(predsMat)[1] <- "Y"
            # 
            # mod <- glm( Y ~., 
            #             data = predsMat, 
            #             family = "binomial", 
            #             alpha = 0,          # no lasso penalty
            #             lambda = optStackParam, # ridge penalty
            #             standardize = FALSE ) # intercept included in design matrix
            # 
            # w_stack <- as.vector( mod$coefficients ) # optimal stacking weights
            # 
            # # include intercept in stacking
            # 
            # # need tochange this if regularization added or w >= 0 (NNLS)--i.e., use CVXR
            # 
            # # predicted values
            # #print("predsmat")
            # predsMat <- as.data.frame( as.matrix(dat[,-c(1,2)]) %*% t(as.matrix( ens_betaMat)))
            # #dim(predsMat)
            # #print("cbind Y")
            # # print(length(Y))
            # predsMat <- cbind(dat$Y, predsMat)
            # #print("colnames Y")
            # colnames(predsMat)[1] <- "Y"
        
        
        ## regular logistic regression without stacking
        # mod <- glm(Y ~., data = predsMat,
        #                 family = "binomial") # train logistic regression-- # intercept included in design matrix
        # w_stack <- mod$coefficients
        # rm(mod)
        #print(predsMat)
        if( is.null(solverName) ){
            w_stack <- as.vector( logisticNNLS(predsMat, w_start = NULL, nnls = nnlsInd, stackInt = TRUE, 
                                               psi = 0.1, simplexInd = simpInd, solverName = NULL, 
                                               accelInd = FALSE, tau = 0.000001, alpha = stackTune, itrMax = 1000) )
        }else if(solverName == "glmnet"){
            w_stack <- as.vector( suppressWarnings( logisticNNLS(predsMat, stackInt = TRUE, 
                                                                 nnls = nnlsInd, 
                                                                 psi = stackPsi,
                                                                 simplexInd = simpInd, solverName = "glmnet",
                                                                 alpha = stackTune) ) ) # solverName = "ECOS"
        }
        
        
        #
        rm(predsMat)
    }else{
        # no stacking intercept
        
        if(tuningInd){
            
            # tune hyperparameter
            p <- ncol(dat) - 2 # number of covariates including intercept
            aucMat <- matrix(ncol = length(alpha_vec), nrow = K ) # store errors associated with each hyperparameter (tested on other studies)
            
            # set held out
            for(studyNum in 1:K){
                
                # iterate through studies to hold out
                rowIndx <- which(dat$Study == studyNum)
                
                # remove current study and corresponding coefficients
                predsMat <- as.data.frame(
                    as.matrix( dat[rowIndx, -c(1,2)] ) %*% 
                        t( as.matrix(ens_betaMat[-studyNum,] ) )
                )
                
                predsMat <- cbind(Y = dat$Y[rowIndx], predsMat)
                
                #colnames(predsMat)[1] <- "Y"
                
                for(lam in 1:length(alpha_vec)){
                    # alpha vec is the regularization hyperparameter for Ridge penalty
                    tuneParam <- alpha_vec[lam]
                    
                    mod <- glmnet( y = as.vector(predsMat[, 1]), 
                                x = as.matrix(predsMat[, -1]), 
                                family = "binomial", 
                                alpha = 0,          # no lasso penalty
                                lambda = tuneParam, # ridge penalty
                                standardize = TRUE,
                                intercept = FALSE) # intercept included in design matrix
                    
                    w <- as.vector( mod$beta ) # stacking weights
                    preds <- as.matrix(predsMat[, -1]) %*% w # predicted linear predictor (systematic component)
                    probs <- 1 / ( 1 + exp( -preds ) ) # predicted probabilities  (i.e., expit(preds) )
                    aucMat[studyNum, lam] <- auc( as.vector(predsMat[, 1]), probs ) # save AUC for this HO study and this tuning param
                    
                }
                
            }
            
            optIndx <- which.max( colMeans( aucMat ) ) # average across held out studies
            stackTune <- alpha_vec[optIndx]
            alpha_opt <- stackTune
        }
        
        predsMat <- as.data.frame( as.matrix(dat[,-c(1,2)]) %*% t(as.matrix( ens_betaMat)))
        predsMat <- cbind(Y = dat$Y, predsMat)
        #colnames(predsMat)[1] <- "Y"
        # print(predsMat)
        
        # normal logistic regression w/o NNLS
        # mod <- glm(Y ~. -1, data = predsMat,
        #                 family = "binomial") # train logistic regression-- # intercept included in design matrix
        # w_stack <- mod$coefficients
        # rm(mod)
        
        
        
        # if(GD){
        #     w_stack <- as.vector(logisticNNLS(predsMat, w_start = NULL, nnls = TRUE, stackInt = FALSE, 
        #                                       psi = 0.1, simplexInd = FALSE, solverName = NULL, 
        #                                       accelInd = FALSE, tau = 0.000001, alpha = 0, itrMax = 10000)) 
        # }else{
        
        w_stack <- as.vector( suppressWarnings(logisticNNLS(predsMat, 
                                                            stackInt = FALSE, 
                                                            nnls = nnlsInd, 
                                                            psi = stackPsi,
                                                            solverName = solverName,
                                                            stackTune = alpha_opt) ) )
        #}
        
        
        rm(predsMat)
        
    }
    
    obj1 <- objLogisticEns(data = dat, beta = ens_betaMat, wVal = w_stack, 
                           lambdaVal = lambda, penVec = alphaOpt_SSL, stackInt = TRUE, 
                           stackTune = alpha_opt) # objective value
    # print(paste0("final obj: ", obj1))
    
    return(list(w_stack, ens_betaMat, alpha_opt, obj1, alphaOpt_SSL))
}


accLogistic <- function(X, Y, beta_r,  alpha_r = 0.5){
    # calculates RMSE
    probs <-  1 / (1 + exp(-as.matrix(X) %*% as.matrix(beta_r))) # predicted probabilities
    preds <- I(probs >= alpha_r) # predicted outcome (0 or 1)
    acc <- sum( I(Y == preds) ) / length(Y) # accuracy
    
    return(acc)
}

accCalc <- function(probs, Y,  alpha_r = 0.5, expitInd = TRUE){
    # calculates RMSE
    if(expitInd){
        preds <- I( expit( probs) >= alpha_r ) # predicted outcome (0 or 1)
    }else{
        preds <- I( probs >= alpha_r ) # predicted outcome (0 or 1)
    }
    
    acc <- sum( I(Y == preds) ) / length(Y) # accuracy
    
    return(acc)
}

RMSEcalc <- function(X, Y, beta, intInd = TRUE, weights = NULL){
    # calculates RMSE
    # beta is K x p matrix
    
    if(is.null(weights)){
        # if no weighrs
        if(!intInd){
            
            sqrt( mean( ((Y) - rowMeans( as.matrix(X) %*% t( as.matrix(beta) ) )     )^2     )     )
            
        }else{
            # if intercept
            sqrt( mean( (Y - rowMeans( as.matrix(cbind(1, X) ) %*% t( as.matrix(beta) ) )   )^2  )   )
        }
    }else{
        # if stacking/CPSW weights
        if(!intInd){
            
            sqrt( mean( ((Y) - as.matrix(X) %*% t( as.matrix(beta) ) %*% as.vector(weights)     )^2     )     )
            
        }else{
            # if intercept
            sqrt( mean( ((Y) - cbind(1, as.matrix(cbind(1, X) ) %*% 
                                                      t( as.matrix(beta) ) ) %*% 
                                                            as.vector(weights) )^2  )  )
        }
        
        
        
    }

    
}


# study regularizer

studyReg <- function(dat, divPen = 0, penVec = NULL){
    
    K <- length(unique(dat[,1]))
    
    if(is.null(penVec)){
        penVec <- rep(0, K) # set each study specific penalty hyperparameter to 0 if not specified by user
    }
    
    ########################
    # Generate C Matrix
    ########################
    
    CMat <- matrix(0, ncol = K, nrow = choose(K, 2) )
    
    cr <- er <- 0 # keep track of current row
    for(study in 1:(K-1) ){
        
        cr <- er + 1 # lower bound of rows
        er <- cr + K - study - 1 # upper bound of rows
        
        CMat[cr:er, study] <- 1
        
        ident <- diag(-1, K - study) # identity matrix of corresponding size
        
        CMat[cr:er, (1 + study):K] <- ident
    }
    
    
    ########################
    # Optimize
    ########################
    dat <- as.matrix(cbind( dat[,c(1,2)], 1, dat[,-c(1,2)]) ) # add ones for intercept
    
    library(CVXR)
    
    beta <- Variable(K, ncol(dat) - 2) # take away 1 for outcome and 1 for study 
    obj <- 0 # set to 0 to begin with and iteratively add in
    
    # study rows
    for(z in 1:K){
        
        j <- unique(dat[,1])[z] # study label
        
        studyIndx <- which(dat[,1] == j) # rows in this study
        obj <- obj + sum(  (dat[studyIndx, 2]  - 
                                dat[studyIndx,-c(1,2)] %*% vec( beta[z, ])  )^2 ) + 
                                    penVec[z] * norm1( vec( beta[z, ]) )
        
    }
    
    
    obj <- obj + divPen * cvxr_norm(CMat %*% beta, 1) / K # study regularizer term
    
    prob <- Problem(Minimize(obj)) # minimize loss
    result <- solve( prob )
    betaOpt <- result$getValue(beta)
    
    return(betaOpt) 
} 


# study regularizer

studyRegStack <- function(dat, divPen = 0, penVec = NULL, 
                          weights = NULL, alpha = 0.5){
    # weights null means average weights
    # alpha is tuning parameter for OEC
    
    K <- length(unique(dat[,1]))
    
    if(is.null(penVec)){
        penVec <- rep(0, K) # set each study specific penalty hyperparameter to 0 if not specified by user
    }
    
    if(is.null(weights)){
        weights <- rep(1/K, K) # set each study specific penalty hyperparameter to 0 if not specified by user
    }else{
        ones <- rep(1, nrow(dat))
    }
    
    
    ########################
    # Generate C Matrix
    ########################
    
    CMat <- matrix(0, ncol = K, nrow = choose(K, 2) )
    
    cr <- er <- 0 # keep track of current row
    for(study in 1:(K-1) ){
        
        cr <- er + 1 # lower bound of rows
        er <- cr + K - study - 1 # upper bound of rows
        
        CMat[cr:er, study] <- 1
        
        ident <- diag(-1, K - study) # identity matrix of corresponding size
        
        CMat[cr:er, (1 + study):K] <- ident
    }
    
    ########################
    # zero out
    ########################
    # zMat <- matrix(1, nrow = nrow(dat), ncol = K)
    # for(study in 1:K ){
    #     
    #     j <- unique(dat[,1])[study] # study label
    #     
    #     studyIndx <- which(dat[,1] == j) # rows in this study
    #     
    #     zMat[studyIndx, study] <- 0
    #     
    # }
    
    ########################
    # Optimize
    ########################
    dat <- as.matrix(cbind( dat[,c(1,2)], 1, dat[,-c(1,2)]) ) # add ones for intercept
    
    library(CVXR)
    
    beta <- Variable(K, ncol(dat) - 2) # take away 1 for outcome and 1 for study 
    obj <- 0 # set to 0 to begin with and iteratively add in
    
    # study rows
    for(z in 1:K){
        
        j <- unique(dat[,1])[z] # study label
        
        studyIndx <- which(dat[,1] == j) # rows in this study
        obj <- obj + sum(  (dat[studyIndx, 2]  - 
                                dat[studyIndx,-c(1,2)] %*% vec( beta[z, ])  )^2 ) + 
            penVec[z] * norm1( vec( beta[z, ]) )
        
    }
    
    
    obj <- obj + divPen * cvxr_norm(CMat %*% beta, 1) / K # study regularizer term
    
    if( length(weights) == K ){
        obj <- (1 - alpha) * obj + alpha * sum( ( dat[,2] -  as.matrix(dat[,-c(1,2)]) %*% t(beta) %*% as.vector( weights ))^2 )
        
    }else{
        
        wInt <- rep(w[1], nrow(dat))
        obj <- (1 - alpha) * obj + alpha * sum( (dat[,2] - 
                                as.matrix(dat[,-c(1,2)]) %*% t(beta) %*% as.vector( weights[-1] ) + wInt )^2 )
        
    }
            
    prob <- Problem(Minimize(obj)) # minimize loss
    result <- solve( prob )
    betaOpt <- result$getValue(beta)
    
    return(betaOpt) 
} 


# same as study Reg bug try to maximize between study differences
# doesnt work because not convex
studyRegMinus <- function(dat, divPen = 0, penVec = NULL){
    
    K <- length(unique(dat[,1]))
    
    if(is.null(penVec)){
        penVec <- rep(0, K) # set each study specific penalty hyperparameter to 0 if not specified by user
    }
    
    ########################
    # Generate C Matrix
    ########################
    
    CMat <- matrix(0, ncol = K, nrow = choose(K, 2) )
    
    cr <- er <- 0 # keep track of current row
    for(study in 1:(K-1) ){
        
        cr <- er + 1 # lower bound of rows
        er <- cr + K - study - 1 # upper bound of rows
        
        CMat[cr:er, study] <- 1
        
        ident <- diag(-1, K - study) # identity matrix of corresponding size
        
        CMat[cr:er, (1 + study):K] <- ident
    }
    
    
    ########################
    # Optimize
    ########################
    dat <- as.matrix(cbind( dat[,c(1,2)], 1, dat[,-c(1,2)]) ) # add ones for intercept
    
    library(CVXR)
    
    beta <- Variable(K, ncol(dat) - 2) # take away 1 for outcome and 1 for study 
    obj <- 0 # set to 0 to begin with and iteratively add in
    
    # study rows
    for(z in 1:K){
        
        j <- unique(dat[,1])[z] # study label
        
        studyIndx <- which(dat[,1] == j) # rows in this study
        obj <- obj + sum(  (dat[studyIndx, 2]  - 
                                dat[studyIndx,-c(1,2)] %*% vec( beta[z, ])  )^2 ) + 
            penVec[z] * norm1( vec( beta[z, ]) )
        
    }
    
    
    #obj <- obj - divPen * log( cvxr_norm(CMat %*% beta, 1) )/ K # study regularizer term
    obj <- obj - divPen * 2 * log( sum(CMat %*% beta ) ) / K # study regularizer term
    
    prob <- Problem(Minimize(obj)) # minimize loss
    result <- solve( prob ) # , solver = "OSQP" 
    betaOpt <- result$getValue(beta)
    
    return(betaOpt) 
} 



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
        # w coefficients 
        # wT <-  as.vector(mod$beta)
        # wT0 <- as.vector(mod$a0)
        # 
        # w0 <- c(wT0, wT)
        
        
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
        
        # if(intercept){
        #     step <- 1/ max(eigen(XX + lambda * Diagonal(x = c(0, rep(1, ncol(XX) - 1) ) ))$values)
        # }else{
        #     step <- 1/ max(eigen(XX + lambda * Diagonal( ncol(XX) )   )$values)
        # }
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
        
        # if(intercept){
        #     step <- 1/ max(eigen(XX + lambda * Diagonal(x = c(0, rep(1, ncol(XX) - 1) ) ))$values)
        # }else{
        #     step <- 1/ max(eigen(XX + lambda * Diagonal( ncol(XX) )   )$values)
        # }
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



#### uses solver --- slow
stack_lim <- function(y,
                      x,
                      low = -Inf,
                      up = Inf,
                      lambda = 0,
                      intercept = TRUE,
                      warmStart = TRUE,
                      simplex = FALSE,
                      solver = "OSQP",
                      u1 = 1 / nrow(x),
                      u2 = 1){
    library(CVXR)
    n <- nrow(x) # sample size
    p <- ncol(x)
    
    if(low == up){
        lambda <- 0
    }
    
    if(lambda == 0){
        # this functionally makes the covariate standardization not happen which is unnecessary when lambda = 0, and when low == up then its necessary
        mean_x <- rep(0, p)
        sd_x <- rep(1, p)
        
    }else{
        mean_x <- colMeans(x)
        sd_x <- sqrt( apply(x, 2, var) * (n-1)/n )
        
        # scale hyperparamter
        sd_y <- as.numeric( sqrt(var(y)*(n-1)/n) ) # standard deviation of y (MLE)
        lambda <- lambda / sd_y # glmnet standardization of penalty
        
        # scale covariates
        

        for(i in 1:p){
            x[,i] <- (x[,i] - mean_x[i])/sd_x[i] 
        }
        
    }
    
    
    w <- Variable( p )
    constraintsList <- list()
    
    if(low == up){
        constraintsList[length(constraintsList) + 1] <- list(w == up) # standardize
    }else{
        
        if(low != -Inf){
            constraintsList[length(constraintsList) + 1] <- list(w >= low * sd_x) # standardize limits like coefficients
        }   
        
        if(up != Inf){
            constraintsList[length(constraintsList) + 1] <- list(w <= up * sd_x) # standardize limits like coefficients
        }
    }
    
    
    if(intercept){
        w0 <- Variable(1)
        obj <- u1 * sum_squares(  y - w0 - x %*% w   ) / 2 + 
            u2 * lambda * sum_squares(w) / 2
        
        prob <- Problem(Minimize(obj), 
                        constraints = constraintsList)
        
        result <- solve( prob, 
                         warm_start = warmStart, 
                         solver = solver )
        
        wOpt <- result$getValue(w)
        wOpt <- wOpt / sd_x # back to original scale like glmnet
        
        w_0 <- result$getValue(w0) - t(wOpt) %*% mean_x # get intercept back to original scale like glmnet
        
        
        wOpt <- c(w_0, wOpt)
    }else{
        obj <- u1 * sum_squares(  y - x %*% w  ) / 2 + 
            u2 * lambda * sum_squares(w) / 2
        
        prob <- Problem(Minimize(obj), 
                        constraints = constraintsList)
        
        result <- solve( prob, 
                         warm_start = warmStart, 
                         solver = solver )
        
        wOpt <- result$getValue(w)
        wOpt <- wOpt / sd_x # back to original scale like glmnet
    }
    
    
    
    
    return(wOpt)
    
    
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
        wLimits <- 0 #c(-Inf, rep(0, p ) )
        
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
    
    #     stackInt <- FALSE
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
        # pcs$rotation <- cbind(pcs$rotation[, ncol(pcs$rotation)], pcs$rotation[, -ncol(pcs$rotation)])
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
    #lambdas <- rep(tuneParam, K)
    
    # if(pcaInd){
    #     if(stackInt){
    #         w <- cbind(1, pcs) %*% w
    #     }else{
    #         w <- pcs %*% w
    #     }
    #     
    # }
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




################
# Was the RidgeAlt before Glm --works just does not have glment capability
################

ridgeAltNoGlm <- function(data,
                     betaStart,
                     wStart,
                     lambdaVec,
                     mu,
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 10,
                     Avg = FALSE){

    # eta is the parameter that determines convex combination of the losses

    library(glmnet)

    eps <- tol + 1 # initialize above tolerance level


    beta <- beta_Star <- betaStart
    w0 <- wStart[1]

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

        # w0 update
        w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow

        # w update
        if(!Avg){
            # only update if not using average weights -- stacking intercept still is updated
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

            }else{

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

#############
# NewWR
#############
# was RidgeAlt until September 14,2020 (was supposed to be faster) BUT doesn't have correct updates (assumes old w[k] for inverting and doesnt update at each iteration)
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        weights <- Diagonal( x = 1 / sigStack[data$Study] ) #Diagonal( x = rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
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
    # XX_k_list <- vector("list", length = K)
    # Xy_k_list <- vector("list", length = K)
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
    # muMat <- Diagonal( x = rep( mu, K )   )
    
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
# switched back to using September 14, 2020
# This was the RidgeAlt until August 20, 2020
# replaced with faster version
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
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
                        # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
        
        # w update
        # if(!Avg){
        #     # if not using average weights
        #     
        #     if(wUpdate == "glmnet"){
        #         # if use glmnet to update parameters
        #         mod <- glmnet(y = y,
        #                       x = as.matrix(X %*% beta), 
        #                       alpha = 0,
        #                       lambda = mu,
        #                       lower.limits = lowLim,
        #                       standardize = standardize,
        #                       intercept = TRUE,
        #                       weights = diag(weights) )
        #         
        #         # coefficients
        #         w <-  as.vector(mod$beta)
        #         w0 <- as.vector(mod$a0)
        #         
        #     }else{
        #         
        #         # w0 update
        #         w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #         
        #         # w update
        #         bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
        #         w <- bInv %*% (
        #             eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
        #             
        #         )
        #         
        #         
        #         
        #     }
        #     
        # }else{
        #     # if average just update the intercept with closed form expression
        #     # w0 update
        #     w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #     
        # }
        # 
        # if(nnlsInd){
        #     # if there is a nonnegativity constraint on the stacking weights
        #     # do not apply to intercept
        #     
        #     if(projs > 0){
        #         
        #         # if project onto unit simplex
        #         
        #         for(pr in 1:projs){
        #             # project onto unit simplex
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #             # Projection onto Sum of 1
        #             w <- w - ((sum(w) - 1) / K)
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #         }
        #         
        #     }else{
        #         # just standard projection onto non negative orthant -- not onto unit simplex
        #         
        #         # Projection onto Non Negative Orthant
        #         w <- sapply(w, function(x) max(0, x) )
        #     }
        #     
        #     
        # }
        # 
        # 
        # 
        
        eps <- obj0 - objImp[cur]
        
        obj0 <- objImp[cur]
        # if(obj_r < min(objVec)){
        #     # if the objective improves save it as current best
        #     beta_Star <- beta
        #     w_Star <- c(w0, w)
        # }else{
        #     # if it is not better, then do not update the best
        #     
        #     if(objCriteria == TRUE){
        #         # if obj criteria is true, use the best iteration as this iterations value
        #         w <- w_Star[-1]
        #         w0 <- w_Star[1]
        #         beta <- beta_Star
        #     }
        # }
        
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
        
            # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
                                           # /   
    }
    
    
    itr <- cur <- 1
    
    # #############################
    # # Calculate study specific value for below
    # #############################
      lambdaList <- vector("list", length = K)
    # XX_k_list <- vector("list", length = K)
    # Xy_k_list <- vector("list", length = K)
    # SSLindx <- c()
    # 
    # for(k in 1:K){
    #     indx <- which(Study == k)
    #     
    #     if(dataSplit < 1){
    #         
    #         n_k <- length(indx) # n_k 
    #         nSplit <- round(n_k * dataSplit) # n_k split
    #         indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
    #         indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking
    #         
    #     }else{
    #         # use all data for both
    #         indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both
    #         
    #     }
    #     
    #     
    #     
    #     lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
    #     Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    # }
    
    
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
            
            # inv <- solve( 
            #     eta * u1 * w[k]^2 * XX + 
            #         (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )  
            # )
            # 
            # betaTemp[,k] <- inv %*% ( 
            #     eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) + 
            #         (1 - eta) * u3[k] * Xy_k_list[[k]]
            # )
            
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
    
    # removed for speed
    # if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
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
        
        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
        # /   
    }
    
    
    itr <- cur <- 1
    
    # #############################
    # # Calculate study specific value for below
    # #############################
    lambdaList <- vector("list", length = K)
    # XX_k_list <- vector("list", length = K)
    # Xy_k_list <- vector("list", length = K)
    # SSLindx <- c()
    # 
    # for(k in 1:K){
    #     indx <- which(Study == k)
    #     
    #     if(dataSplit < 1){
    #         
    #         n_k <- length(indx) # n_k 
    #         nSplit <- round(n_k * dataSplit) # n_k split
    #         indx_k_SSL <- indx[1:nSplit] # indicies corresponding to SSL
    #         indx_k_Stack <- setdiff(indx, indx_k_SSL) # indicies corresponding to stacking
    #         
    #     }else{
    #         # use all data for both
    #         indx_k_SSL <- indx # if its 1 that corresponds to all data be used in both
    #         
    #     }
    #     
    #     
    #     
    #     lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
    #     Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
    # }
    
    
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
            
            
            # inv <- solve( 
            #     eta * u1 * w[k]^2 * XX + 
            #         (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )  
            # )
            # 
            # betaTemp[,k] <- inv %*% ( 
            #     eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) + 
            #         (1 - eta) * u3[k] * Xy_k_list[[k]]
            # )
            
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
    
    wStart[testIndx + 1] <- 0 # zero out coefficient corresponding to test training set (+1 because of intercept)
    betaStart[,testIndx] <- 0 # zero out coefficients corresponding to test training study
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
    
    if(low <= 0 || up < Inf){
        # use glmnet() for optimization
        # if low >0 but   up < Inf then this will default to glmnet which will prioritize the upper limits
        glmnetInd <- TRUE
    }else{
        glmnetInd <- FALSE
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
        
        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
        # /   
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
        
        
        
        # lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k], p - 1)   ) )#Diagonal( x = c(0, rep( lambdaVec[k] * u4[k], p - 1)   )  )  # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
        # 
        # XX_k_list[[k]] <- t(X[indx_k_SSL,]) %*% X[indx_k_SSL,]
        # Xy_k_list[[k]] <- t(X[indx_k_SSL,])  %*% y[indx_k_SSL]
        
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
    
    # for(k in 1:K){
    #     
    #     indx <- which(Study == k)
    #     lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list <- t(X[indx,]) %*% X[indx,]
    #     Xy_k_list <- t(X[indx,])  %*% y[indx]
    #     
    #     inv_list[[k]] <- as.matrix( solve( 
    #         eta * u1 * w[k]^2 * XX + 
    #             (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )  
    #     ) )
    #     
    #     
    #     mat[,k] <- as.vector(inv_list[[k]] %*% ( 
    #         ( eta * u1 * w[k] *  Xy  ) + 
    #             (1 - eta) * u3[k] * Xy_k_list 
    #     ) )
    #     
    #     mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )
    #     
    #     rm(XX_k_list, Xy_k_list)
    #     
    # }
    
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
            
            # betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + 
            #     mat[,k] - w0 * mat2[,k]
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
        
        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
        # /   
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
        
        # 
        # lambdaList[[k]] <- diag( c(0, rep( lambdaVec[k] * u4[k], p - 1)   ) )#Diagonal( x = c(0, rep( lambdaVec[k] * u4[k], p - 1)   )  )  # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
       
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
    
    # for(k in 1:K){
    #     
    #     indx <- which(Study == k)
    #     lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list <- t(X[indx,]) %*% X[indx,]
    #     Xy_k_list <- t(X[indx,])  %*% y[indx]
    #     
    #     inv_list[[k]] <- as.matrix( solve( 
    #         eta * u1 * w[k]^2 * XX + 
    #             (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )  
    #     ) )
    #     
    #     
    #     mat[,k] <- as.vector(inv_list[[k]] %*% ( 
    #         ( eta * u1 * w[k] *  Xy  ) + 
    #             (1 - eta) * u3[k] * Xy_k_list 
    #     ) )
    #     
    #     mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )
    #     
    #     rm(XX_k_list, Xy_k_list)
    #     
    # }
    
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
            
            # betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + 
            #     mat[,k] - w0 * mat2[,k]
            
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
# USE PCA on the LHS (for stacking)
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
    
    # removed for speed
    # if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    # wStart[-1] <- pcaMat %*% wStart[-1] # PCA 
    
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
        
        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
        # /   
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
    # inv_list <- vector("list", length = K)
    # mat <- mat2 <- matrix(ncol = K, nrow = p )
    
    # for(k in 1:K){
    #     
    #     indx <- which(Study == k)
    #     lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list <- t(X[indx,]) %*% X[indx,]
    #     Xy_k_list <- t(X[indx,])  %*% y[indx]
    #     
    #     inv_list[[k]] <- as.matrix( solve( 
    #         eta * u1 * w[k]^2 * XX + 
    #             (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )  
    #     ) )
    #     
    #     
    #     mat[,k] <- as.vector(inv_list[[k]] %*% ( 
    #         ( eta * u1 * w[k] *  Xy  ) + 
    #             (1 - eta) * u3[k] * Xy_k_list 
    #     ) )
    #     
    #     mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )
    #     
    #     rm(XX_k_list, Xy_k_list)
    #     
    # }
    
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
            
            # betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + 
            #     mat[,k] - w0 * mat2[,k]
            
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
                    #wT <-  wT[-1] #as.vector(mod$beta)

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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
    studyVec <- studies <- sort( unique(data$Study) ) # added in sort 9/10/20 because of study indexing problems)
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
        
        # *******if dataSplit != then you need to adjust these weights (with stackIndex) to equal 1 *******
        # weights[Stackindx, Stackindx] <- weights[Stackindx, Stackindx] / sum(diag(weights[Stackindx, Stackindx]))
        
        # /   
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
    inv_list <- vector("list", length = K)
    mat <- mat2 <- matrix(ncol = K, nrow = p )
    
    # for(k in 1:K){
    #     
    #     indx <- which(Study == k)
    #     lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list <- t(X[indx,]) %*% X[indx,]
    #     Xy_k_list <- t(X[indx,])  %*% y[indx]
    #     
    #     inv_list[[k]] <- as.matrix( solve( 
    #         eta * u1 * w[k]^2 * XX + 
    #             (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )  
    #     ) )
    #     
    #     
    #     mat[,k] <- as.vector(inv_list[[k]] %*% ( 
    #         ( eta * u1 * w[k] *  Xy  ) + 
    #             (1 - eta) * u3[k] * Xy_k_list 
    #     ) )
    #     
    #     mat2[,k] <- eta * u1 * w[k] * as.vector(inv_list[[k]] %*% ( X_rowSums ) )
    #     
    #     rm(XX_k_list, Xy_k_list)
    #     
    # }
    
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
            
            # betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + 
            #     mat[,k] - w0 * mat2[,k]
            
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
                
                # if(glmnetInd){
                #     # use glmnet for optimization
                #     mod <- glmnet(y = y,
                #                   x = as.matrix(X %*% beta), 
                #                   alpha = 0,
                #                   lambda = mu,
                #                   lower.limits = low,
                #                   upper.limits = up,
                #                   standardize = standardize,
                #                   intercept = TRUE,
                #                   weights = weights, 
                #                   thresh = 1e-10)
                #     
                #     # w coefficients 
                #     wT <-  as.vector(mod$beta)
                #     wT0 <- as.vector(mod$a0)
                #     
                # }else{
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
    
    # removed for speed
    # if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
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
                # wT0 <- mod[1]
                
                #     glmnet(y = y,
                #               x = as.matrix(X %*% beta), 
                #               alpha = 0,
                #               lambda = mu,
                #               lower.limits = lowLim,
                #               standardize = standardize,
                #               intercept = TRUE,
                #               weights = diag(weights), 
                #               thresh = 1e-10)
                # 
                # # w coefficients 
                # wT <-  as.vector(mod$beta)
                # wT0 <- as.vector(mod$a0)
                
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                
                # wT0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow # eta remove 8/14/20
                
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
                    #w0 <- wT0
                    
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
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
# fix stacking weights
ridgeAltFix2 <- function(data, 
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
    
    library(CVXR)
    # solver for w update
    
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
            
            inv <- solve( 
                eta * u1 * w[k]^2 * XX + 
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )  
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
            
            
        } # end for loop for beta updates
        
        # w update
        if(!Avg){
            
            mod <-       stackLim(y = y,
                                  x = as.matrix(X %*% beta),
                                  low = low,
                                  # up = up,  doesnt currently work for upper lims
                                  lambda = mu,
                                  w0 = c(w0, w), # warm start
                                  tol = stackTol,  # better results from 1e-12
                                  u1 = u1,
                                  u2 = u2,
                                  weights = weights)
            
            # mod <-              stack_lim(y = y,
            #                               x = as.matrix(X %*% beta),
            #                               low = low,
            #                               up = up,
            #                               lambda = mu,
            #                               intercept = TRUE,
            #                               warmStart = TRUE,
            #                               simplex = simplex,
            #                               solver = solverNm,
            #                               u1 = u1,
            #                               u2 = u2)
            
            wT <- mod[-1]
            wT0 <- mod[1]
            #     glmnet(y = y,
            #               x = as.matrix(X %*% beta), 
            #               alpha = 0,
            #               lambda = mu,
            #               lower.limits = lowLim,
            #               standardize = standardize,
            #               intercept = TRUE,
            #               weights = diag(weights), 
            #               thresh = 1e-10)
            # 
            # # w coefficients 
            # wT <-  as.vector(mod$beta)
            # wT0 <- as.vector(mod$a0)
            
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
            # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
        # end w update
        
        ######################################
        
        eps <- obj0 - objImp[cur]
        
        obj0 <- objImp[cur]

        
    }         
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))
    
}


#############
# NewWR
#############
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
    # #############################
    # # Calculate study specific value for below
    # #############################
    # lambdaList <- vector("list", length = K)
    # # XX_k_list <- vector("list", length = K)
    # # Xy_k_list <- vector("list", length = K)
    # SSLindx <- c()
    # Stackindx <- c()
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
    
    # if(sampSzWeight < 6){
    #     weights <- NULL
    # }else{
    #     weights <- diag(weights) # vectorize this
    # }
    
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
    
    # #########################
    # # Inverse Matrix
    # #########################
    # inv_list <- vector("list", length = K)
    # mat <- matrix(ncol = K, nrow = p )
    # 
    # for(k in 1:K){
    #     
    #     indx <- which(Study == k)
    #     lambdaList[[k]] <- Diagonal( x = c(0, rep( lambdaVec[k], p - 1)   )  ) # put 0 in first position so that we do not penalize intercept for study specific L2 penalty
    #     
    #     XX_k_list <- t(X[indx,]) %*% X[indx,]
    #     Xy_k_list <- t(X[indx,])  %*% y[indx]
    #     
    #     inv_list[[k]] <- as.matrix( solve( 
    #         eta * u1 * w[k]^2 * XX + 
    #             (1 - eta) * ( XX_k_list * u3[k] + lambdaList[[k]] * u4[k] )  
    #     ) )
    #     
    #     
    #     mat[,k] <- as.vector(inv_list[[k]] %*% ( 
    #         ( eta * u1 * w[k] * ( Xy - w0 * X_rowSums ) ) + 
    #             (1 - eta) * u3[k] * Xy_k_list 
    #     ) )
    #     
    #     rm(XX_k_list, Xy_k_list)
    #     
    # }
    
    if(dataSplit == 1)    Stackindx <- NULL
    
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
            
            
            # betaTemp[,k] <- -eta * u1 * w[k] * inv_list[[k]] %*% ( XX %*% beta[,-k] %*% w[-k] ) + mat[,k]
            
            inv <- solve(
                eta * u1 * w[k]^2 * XX +
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]]  ) # * u4[k] above
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
                                      # up = up,  doesnt currently work for upper lims
                                      lambda = mu,
                                      w0 = c(w0, w), # warm start
                                      tol = stackTol,  # better results from 1e-12
                                      u1 = u1,
                                      u2 = u2,
                                      weights = weights)
                
                # mod <-              stack_lim(y = y,
                #                               x = as.matrix(X %*% beta),
                #                               low = low,
                #                               up = up,
                #                               lambda = mu,
                #                               intercept = TRUE,
                #                               warmStart = TRUE,
                #                               simplex = simplex,
                #                               solver = solverNm,
                #                               u1 = u1,
                #                               u2 = u2)
                
                    wT <- mod[-1]
                    wT0 <- mod[1]
                }
                #     glmnet(y = y,
                #               x = as.matrix(X %*% beta), 
                #               alpha = 0,
                #               lambda = mu,
                #               lower.limits = lowLim,
                #               standardize = standardize,
                #               intercept = TRUE,
                #               weights = diag(weights), 
                #               thresh = 1e-10)
                # 
                # # w coefficients 
                # wT <-  as.vector(mod$beta)
                # wT0 <- as.vector(mod$a0)
                
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
# stopped using this on September 14, 2020 because doesnt change based on w[k]
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
# fix stacking weights
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
    # XX_k_list <- vector("list", length = K)
    # Xy_k_list <- vector("list", length = K)
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

    # if(sampSzWeight < 6){
    #     weights <- NULL
    # }else{
    #     weights <- diag(weights) # vectorize this
    # }
    
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

            # inv <- solve( 
            #     eta * u1 * w[k]^2 * XX + 
            #         (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )  
            # )
            # 
            # betaTemp[,k] <- inv %*% ( 
            #     eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) + 
            #         (1 - eta) * u3[k] * Xy_k_list[[k]]
            # )
            
            
            
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
                                      # up = up,  doesnt currently work for upper lims
                                      lambda = mu,
                                      w0 = c(w0, w), # warm start
                                      tol = stackTol,  # better results from 1e-12
                                      u1 = u1,
                                      u2 = u2,
                                      weights = weights)
                
                # mod <-              stack_lim(y = y,
                #                               x = as.matrix(X %*% beta),
                #                               low = low,
                #                               up = up,
                #                               lambda = mu,
                #                               intercept = TRUE,
                #                               warmStart = TRUE,
                #                               simplex = simplex,
                #                               solver = solverNm,
                #                               u1 = u1,
                #                               u2 = u2)
                
                wT <- mod[-1]
                wT0 <- mod[1]
                #     glmnet(y = y,
                #               x = as.matrix(X %*% beta), 
                #               alpha = 0,
                #               lambda = mu,
                #               lower.limits = lowLim,
                #               standardize = standardize,
                #               intercept = TRUE,
                #               weights = diag(weights), 
                #               thresh = 1e-10)
                # 
                # # w coefficients 
                # wT <-  as.vector(mod$beta)
                # wT0 <- as.vector(mod$a0)
                
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# started using July 9, 2020 after benchmarking
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
        
        # w update
        # if(!Avg){
        #     # if not using average weights
        #     
        #     if(wUpdate == "glmnet"){
        #         # if use glmnet to update parameters
        #         mod <- glmnet(y = y,
        #                       x = as.matrix(X %*% beta), 
        #                       alpha = 0,
        #                       lambda = mu,
        #                       lower.limits = lowLim,
        #                       standardize = standardize,
        #                       intercept = TRUE,
        #                       weights = diag(weights) )
        #         
        #         # coefficients
        #         w <-  as.vector(mod$beta)
        #         w0 <- as.vector(mod$a0)
        #         
        #     }else{
        #         
        #         # w0 update
        #         w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #         
        #         # w update
        #         bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
        #         w <- bInv %*% (
        #             eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
        #             
        #         )
        #         
        #         
        #         
        #     }
        #     
        # }else{
        #     # if average just update the intercept with closed form expression
        #     # w0 update
        #     w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #     
        # }
        # 
        # if(nnlsInd){
        #     # if there is a nonnegativity constraint on the stacking weights
        #     # do not apply to intercept
        #     
        #     if(projs > 0){
        #         
        #         # if project onto unit simplex
        #         
        #         for(pr in 1:projs){
        #             # project onto unit simplex
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #             # Projection onto Sum of 1
        #             w <- w - ((sum(w) - 1) / K)
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #         }
        #         
        #     }else{
        #         # just standard projection onto non negative orthant -- not onto unit simplex
        #         
        #         # Projection onto Non Negative Orthant
        #         w <- sapply(w, function(x) max(0, x) )
        #     }
        #     
        #     
        # }
        # 
        # 
        # 
        
        eps <- obj0 - objImp[cur]
        
        obj0 <- objImp[cur]
        # if(obj_r < min(objVec)){
        #     # if the objective improves save it as current best
        #     beta_Star <- beta
        #     w_Star <- c(w0, w)
        # }else{
        #     # if it is not better, then do not update the best
        #     
        #     if(objCriteria == TRUE){
        #         # if obj criteria is true, use the best iteration as this iterations value
        #         w <- w_Star[-1]
        #         w0 <- w_Star[1]
        #         beta <- beta_Star
        #     }
        # }
        
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
        u4 <- 1 / sigK # standardize the hyperparamter because of tuning
    }else if(sampSzWeight == 6){
        # just variance adjustment on LHS
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
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
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] ) # * u4[k]   --- above
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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
        
        # w update
        # if(!Avg){
        #     # if not using average weights
        #     
        #     if(wUpdate == "glmnet"){
        #         # if use glmnet to update parameters
        #         mod <- glmnet(y = y,
        #                       x = as.matrix(X %*% beta), 
        #                       alpha = 0,
        #                       lambda = mu,
        #                       lower.limits = lowLim,
        #                       standardize = standardize,
        #                       intercept = TRUE,
        #                       weights = diag(weights) )
        #         
        #         # coefficients
        #         w <-  as.vector(mod$beta)
        #         w0 <- as.vector(mod$a0)
        #         
        #     }else{
        #         
        #         # w0 update
        #         w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #         
        #         # w update
        #         bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
        #         w <- bInv %*% (
        #             eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
        #             
        #         )
        #         
        #         
        #         
        #     }
        #     
        # }else{
        #     # if average just update the intercept with closed form expression
        #     # w0 update
        #     w0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
        #     
        # }
        # 
        # if(nnlsInd){
        #     # if there is a nonnegativity constraint on the stacking weights
        #     # do not apply to intercept
        #     
        #     if(projs > 0){
        #         
        #         # if project onto unit simplex
        #         
        #         for(pr in 1:projs){
        #             # project onto unit simplex
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #             # Projection onto Sum of 1
        #             w <- w - ((sum(w) - 1) / K)
        #             
        #             # Projection onto Non Negative Orthant
        #             w <- sapply(w, function(x) max(0, x) )
        #         }
        #         
        #     }else{
        #         # just standard projection onto non negative orthant -- not onto unit simplex
        #         
        #         # Projection onto Non Negative Orthant
        #         w <- sapply(w, function(x) max(0, x) )
        #     }
        #     
        #     
        # }
        # 
        # 
        # 
        
        eps <- obj0 - objImp[cur]
        
        obj0 <- objImp[cur]
        # if(obj_r < min(objVec)){
        #     # if the objective improves save it as current best
        #     beta_Star <- beta
        #     w_Star <- c(w0, w)
        # }else{
        #     # if it is not better, then do not update the best
        #     
        #     if(objCriteria == TRUE){
        #         # if obj criteria is true, use the best iteration as this iterations value
        #         w <- w_Star[-1]
        #         w0 <- w_Star[1]
        #         beta <- beta_Star
        #     }
        # }
        
    }         
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr, objImp = objImp))
    
}


##########################
# Zero Out OEC
##########################
# randomized order
# only uses glmnet -- does not allow for average weights
# does not allow for stackIndx or obsW
# stopped using 10/4/20

ridgeAlt0_OLD <- function(data, 
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
                     weights = NULL){
    
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
    Imat <- Diagonal(N)
    
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
    I_k  <- vector("list", length = K)
    XX_k_list <- vector("list", length = K)
    Xy_k_list <- vector("list", length = K)
    SSLindx <- c()
    Stackindx <- c()
    
    for(k in 1:K){
        
        # zero out Identity matrices
        indx <- which(Study == k)
        I_k[[k]] <- Imat
        I_k[[k]][indx, indx] <- 0 # zero out diagonal entries associated with this study
        
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
    
    rm(Imat)
    
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
    muMat <- diag( rep( mu, K )   )
    y <- y[Stackindx]
    
    # NOW SET TO NULL for objective evaluation below -- above Stackindx was set to 1:nrow(X) -- keep both!
    stkIndex <- Stackindx
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
            
            sumMat <- 0
            # for a part below
            for(sK in seq(1,K)[-k] ){
                sumMat <- sumMat + w[sK] * t(X[stkIndex,]) %*% 
                    I_k[[k]] %*% I_k[[sK]] %*% X[stkIndex,] %*% beta[,sK]
            }
            
            inv <- solve( 
                eta * u1 * w[k]^2 * t(X[stkIndex,]) %*% I_k[[k]] %*% X[stkIndex,] + 
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k] )  
            )
            
            betaTemp[,k] <- as.vector( inv %*% ( 
                                        eta * u1 * w[k] * ( t(X[stkIndex,]) %*% I_k[[k]] %*% y[stkIndex]  - 
                                        w0 * rowSums( t(X[stkIndex,] ) %*% I_k[[k]] ) - 
                                        sumMat ) + 
                                        (1 - eta) * u3[k] * Xy_k_list[[k]]
                                        )
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
                
                # zero out
                predsMat <- as.matrix(X %*% beta)
                
                for(studyNum in 1:K){ 
                    sIndex <- which(data$Study == studyNum)
                    predsMat[sIndex, studyNum ] <- 0 # zero out 
                }
                
                
                mod <- glmnet(y = y,
                              x = predsMat, 
                              alpha = 0,
                              lambda = mu,
                              lower.limits = lowLim,
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights), 
                              thresh = 1e-20)
                
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
                # wT0 <- eta * ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
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

#############
# glmnet
#############
# allows obsWeights and StackIndex but old convergence criteria # usedJuly 9, 2020
# update after going through all K studies
ridgeAltNew <- function(data, 
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
                          weights = NULL){
    
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
    
    if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
    }
    
    
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
    
    objVec <- objOEC(data = data, 
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
        for(k in 1:K){
            
            inv <- solve( 
                eta * u1 * w[k]^2 * XX + 
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k])  
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
                            weights = weights
            )
            
            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp

            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }
            
            objVec <- c(objVec, obj_k) # add objective even if no improvement
            
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
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights) )
                
                # coefficients
                w <-  as.vector(mod$beta)
                w0 <- as.vector(mod$a0)
                
            }else{
                
                # w0 update
                w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                
                # w update
                bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
                w <- bInv %*% (
                    eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
                    
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
                        StackIndx = Stackindx,
                        sampSzWeight = sampSzWeight,
                        weights = weights
        )
        
        
        
        if(obj_r < min(objVec)){
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
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))
    
}

#############
# glmnet
#############
# # added July 9, 2020 
# used this before benchmarking and changing to randomize order and update w at each time
# and before checking objective at each round
# THIS ALLOWS FOR obsWeights and StackIndx (using some observations for stacking)
# BUT you would have to change the objective function to the objOECweighted
 # updates beta after each w
ridgeAltNewW <- function(data, 
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
                        weights = NULL){
    
    # eta is the parameter that determines convex combination of the losses
    
    # # rename studies from 1:K 
    # studies <- unique(data$Study)
    # StudyV <- vector(length = K)
    # for(i in 1:length(data$Study)){
    #     StudyV[i] <-  which(studies == data$Study[i]) # switch to numerical label from 1:K
    # }
    # 
    # data$Study <- as.numeric(StudyV)
    # rm(StudyV)
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
    
    if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
    }
    
    
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
    
    objVec <- objOEC(data = data, 
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
        for(k in 1:K){
            
            inv <- solve( 
                eta * u1 * w[k]^2 * XX + 
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k])  
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
                            weights = weights
            )
            
            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                beta_Star <- beta <- betaTemp
                
            }else{
                betaTemp <- beta # set betaTemp to what it was previously before update
            }
            
            objVec <- c(objVec, obj_k) # add objective even if no improvement
            
            
            # w update
            mod <- glmnet(y = y,
                          x = as.matrix(X %*% beta), 
                          alpha = 0,
                          lambda = mu,
                          lower.limits = lowLim,
                          standardize = standardize,
                          intercept = TRUE,
                          weights = diag(weights) )
            
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
                            weights = weights
            )
            
            if(obj_k < min(objVec) ){
                # if the objective improves save it as current best
                w_Star <- c(wT0, wT)
                w <- wT
                w0 <- wT0
                
            }
            
            objVec <- c(objVec, obj_k) # add objective even if no improvement
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
                              standardize = standardize,
                              intercept = TRUE,
                              weights = diag(weights) )
                
                # coefficients
                w <-  as.vector(mod$beta)
                w0 <- as.vector(mod$a0)
                
            }else{
                
                # w0 update
                w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                
                # w update
                bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
                w <- bInv %*% (
                    eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
                    
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
                        StackIndx = Stackindx,
                        sampSzWeight = sampSzWeight,
                        weights = weights
        )
        
        
        
        if(obj_r < min(objVec)){
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
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))
    
}

#############
# glmnet
#############
# update after going through all K studies
# changed July 9, 2020 
# used this before benchmarking and changing to randomize order and update w at each time
# and before checking objective at each round
# THIS ALLOWS FOR obsWeights and StackIndx (using some observations for stacking)
# BUT you would have to change the objective function to the objOECweighted
ridgeAltobsWeight <- function(data, 
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
                     weights = NULL){
    
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
    
    if(is.null(weights))      weights <- matrix( diag( nrow( data ) ) ) # make the sample size weights the identity if its null
    
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
    }

    
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
    
    objVec <- objOEC(data = data, 
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
            # beta update
            for(k in 1:K){
                
                inv <- solve( 
                    eta * u1 * w[k]^2 * XX + 
                        (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k])  
                )
                
                beta[,k] <- inv %*% ( 
                    eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) + 
                        (1 - eta) * u3[k] * Xy_k_list[[k]]
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
                                  standardize = standardize,
                                  intercept = TRUE,
                                  weights = diag(weights) )
                    
                    # coefficients
                    w <-  as.vector(mod$beta)
                    w0 <- as.vector(mod$a0)
                    
                }else{
                    
                    # w0 update
                    w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                    
                    # w update
                    bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
                    w <- bInv %*% (
                        eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) %*% weights ) ) 
                        
                    )
                    
                    
                    
                }
                
            }else{
                # if average just update the intercept with closed form expression
                # w0 update
                w0 <-  ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                
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
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight,
                            weights = weights
                            )
            
            
            
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
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))
    
}


#############
# glmnet
#############
 # no sample specific weights for stacking portion of loss
ridgeAltUnweighted <- function(data, 
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
                     sampSzWeight = 4){
    
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
  
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
    }
    
    
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
    Xnrow <- N <- nrow(X) # number of rows in all studies
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
                     StackIndx = Stackindx,
                     sampSzWeight = sampSzWeight)
    
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
        # beta update
        for(k in 1:K){
            
            inv <- solve( 
                eta * u1 * w[k]^2 * XX + 
                    (1 - eta) * ( XX_k_list[[k]] * u3[k] + lambdaList[[k]] * u4[k])  
            )
            
            beta[,k] <- inv %*% ( 
                eta * u1 * w[k] * ( Xy - w0 * X_rowSums - XX %*% beta[,-k] %*% w[-k] ) + 
                    (1 - eta) * u3[k] * Xy_k_list[[k]]
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
                              standardize = standardize,
                              intercept = TRUE)
                
                # coefficients
                w <-  as.vector(mod$beta)
                w0 <- as.vector(mod$a0)
                
            }else{
                
                # w0 update
                w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                
                # w update
                bInv <- solve(    eta * ( u1 * t(beta) %*% XX %*% beta + u2 * muMat )    )
                w <- bInv %*% (
                    eta * u1 * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) ) ) 
                    
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
                        StackIndx = Stackindx,
                        sampSzWeight = sampSzWeight)
        
        
        
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
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))
    
}


# used this before the multi-objective
ridgeAlt_singleObj <- function(data, 
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
                     sampSzWeight = TRUE){
    
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
    
    library(glmnet)
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
    
    nVec <- as.vector( table(data$Study) ) # number of rows in each study
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
    Xnrow <- N <- nrow(X) # number of rows in all studies
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
                     StackIndx = Stackindx,
                     sampSzWeight = sampSzWeight)
    
    #############################
    # Iteratively update
    #############################
    
    # check to see how objective is weighted
    # sampSzWeight
    
   nVec # the regularization penalty is scaled so it is n_k * lambda_k (elementwise product)
    
    if(sampSzWeight){
        # if scale LHS of objective by 1 / N to account for variable sample sizes
        while(eps > tol){
            
            #print(paste("iter", itr))
            itr <- itr + 1 
            # beta update
            for(k in 1:K){
                
                inv <- solve( 
                    eta * w[k]^2 / N * XX + 
                        (1 - eta) * ( XX_k_list[[k]] / N + lambdaList[[k]] * nVec[k] / N)  
                )
                
                beta[,k] <- inv %*% ( 
                    eta * w[k] * ( Xy / N - w0 / N * X_rowSums - XX %*% beta[,-k] %*% w[-k] / N ) + 
                        (1 - eta)  *  Xy_k_list[[k]] / N
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
                                  standardize = standardize,
                                  intercept = TRUE)
                    
                    # coefficients
                    w <-  as.vector(mod$beta)
                    w0 <- as.vector(mod$a0)
                    
                }else{
                    
                    # w0 update
                    w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                    
                    # w update
                    bInv <- solve(    eta * ( t(beta) %*% XX %*% beta / N + muMat )    )
                    w <- bInv %*% (
                        eta / N * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) ) ) 
                        
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
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight)
            
            
            
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
            
    }else{
        # original objective if you do NOT scale LHS of objective by 1 / N to account for variable sample sizes
        
        while(eps > tol){
            
            #print(paste("iter", itr))
            itr <- itr + 1 
            # beta update
            for(k in 1:K){
                
                inv <- solve( 
                    eta * w[k]^2 / N * XX + 
                        (1 - eta) * ( XX_k_list[[k]] / nVec[k] + lambdaList[[k]] )  
                )
                
                beta[,k] <- inv %*% ( 
                    eta * w[k] * ( Xy / N - w0 / N * X_rowSums - XX %*% beta[,-k] %*% w[-k] / N ) + 
                        (1 - eta)  *  Xy_k_list[[k]] / nVec[k]
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
                                  standardize = standardize,
                                  intercept = TRUE)
                    
                    # coefficients
                    w <-  as.vector(mod$beta)
                    w0 <- as.vector(mod$a0)
                    
                }else{
                    
                    # w0 update
                    w0 <- ( y_sum - sum( X %*% beta %*% w ) ) / Xnrow
                    
                    # w update
                    bInv <- solve(    eta * ( t(beta) %*% XX %*% beta / N + muMat )    )
                    w <- bInv %*% (
                        eta / N * ( t(beta) %*% Xy - w0 * rowSums( t(beta) %*% t(X) ) ) 
                        
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
                            StackIndx = Stackindx,
                            sampSzWeight = sampSzWeight)
            
            
            
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
    }
    
    
    
    
    return(list(beta = beta_Star, w = w_Star, obj = objVec, itrs = itr))
    
}



#############
# glmnet
#############
# was used until 5 -22- 20 until conversation with Rahil about scaline the loss function by 1/N or 1?n_k
ridgeAltUnstandardized <- function(data, 
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
                     standardize = FALSE){
    
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
                              standardize = standardize,
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

ridgeAltConst <- function(data, 
                     betaStart, 
                     wStart, 
                     lambdaVec, 
                     mu, 
                     nnlsInd = TRUE,
                     tol = 0.001,
                     objCriteria = FALSE,
                     eta = 0.5,
                     dataSplit = 1,
                     projs = 3,
                     wUpdate = "glmnet"){
    
    # same as other function but constrained so that the ws sum to 1
    
    # eta is the parameter that determines convex combination of the losses
    
    library(glmnet)
    
    # set for glmnet -- if using nnls for stacking
    if(nnlsInd){
        lowLim <- 0
    }else{
        lowLim <- -Inf
    }
    
    
    eps <- tol + 1 # initialize above tolerance level
    
    
    beta <- beta_Star <- betaStart
    w0 <- wStart[1]
    
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
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
        
    }
    
    
    # if intercept in stacking regression add column of ones to beta
    if(stackInt){
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2] - w[1] - cbind(1, data[,-c(1,2)]) %*% as.matrix( beta ) %*% w[-1] )^2 )
        # obj <- u1 * eta / 2 * sum( ( data[,2] - cbind(1, data[,-c(1,2)]) %*%  cbind(1, beta) %*% w )^2 )   DOESNT MAKE SENSE way calculated objective
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( cbind(1, beta) ) %*% as.vector(w))
    }else{
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2]  -  cbind(1, data[,-c(1,2)]) %*% beta %*% w)^2     )
        
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( beta ) %*% as.vector(w)) 
    }
    
    for(j in 1:K){
        
        indx <- which(Study == j) # rows from this study
        
        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% as.matrix( beta[-1, j] ) )^2     )
        
        
        #     objLinear(X = as.matrix( cbind(1, data[indx,-c(1,2)] ) ), 
        #                                            Y = as.vector( data[indx, 2] ), 
        #                                            beta = as.matrix(beta[,j],
        #                                                             weights = 1
        #                                            )
        # )
        
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
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)
        
        
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / (sigK * nVec) # K x 1 vector
        u4 <- 1 / sigK
        
    }
    
    
    # if intercept in stacking regression add column of ones to beta
    if(stackInt){
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2] - w[1] - cbind(1, data[,-c(1,2)]) %*% beta %*% pcaMat %*% w[-1] )^2 )
        # obj <- u1 * eta / 2 * sum( ( data[,2] - cbind(1, data[,-c(1,2)]) %*%  cbind(1, beta) %*% w )^2 )   DOESNT MAKE SENSE way calculated objective
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( cbind(1, beta) ) %*% as.vector(w))
    }else{
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[,2]  -  cbind(1, data[,-c(1,2)]) %*% beta %*% pcaMat %*% w)^2     )
        
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( beta ) %*% as.vector(w)) 
    }
    
    for(j in 1:K){
        
        indx <- which(Study == j) # rows from this study
        
        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% beta[-1, j] )^2     )
        
        
        #     objLinear(X = as.matrix( cbind(1, data[indx,-c(1,2)] ) ), 
        #                                            Y = as.vector( data[indx, 2] ), 
        #                                            beta = as.matrix(beta[,j],
        #                                                             weights = 1
        #                                            )
        # )
        
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
        u1 <- 1 / length(StackIndx)  # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)
        
        
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(StackIndx)  # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
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
        # obj <- u1 * eta / 2 * sum( ( data[,2] - cbind(1, data[,-c(1,2)]) %*%  cbind(1, beta) %*% w )^2 )   DOESNT MAKE SENSE way calculated objective
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( cbind(1, beta) ) %*% as.vector(w))
    }else{
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2]  -  cbind(1, data[StackIndx,-c(1,2)]) %*% as.matrix( beta ) %*% w)^2     )
        
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( beta ) %*% as.vector(w)) 
    }
    
    for(j in 1:K){
        
        indx <- which(Study == j) # rows from this study
        
        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% as.matrix( beta[-1, j] ) )^2     )
        
        
        #     objLinear(X = as.matrix( cbind(1, data[indx,-c(1,2)] ) ), 
        #                                            Y = as.vector( data[indx, 2] ), 
        #                                            beta = as.matrix(beta[,j],
        #                                                             weights = 1
        #                                            )
        # )
        
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
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)
        
        
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / length(StackIndx)  # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
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
        # obj <- u1 * eta / 2 * sum( ( data[,2] - cbind(1, data[,-c(1,2)]) %*%  cbind(1, beta) %*% w )^2 )   DOESNT MAKE SENSE way calculated objective
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( cbind(1, beta) ) %*% as.vector(w))
    }else{
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / 2 * sum( u1 * ( data[StackIndx,2]  -  cbind(1, data[StackIndx,-c(1,2)]) %*% beta %*% pcaMat %*% w)^2     )
        
        
        # objLinear(X = as.matrix( cbind(1, data[,-c(1,2)]) ), 
        #                             Y = as.matrix( data[,2] ), 
        #                             beta = as.matrix( beta ) %*% as.vector(w)) 
    }
    
    for(j in 1:K){
        
        indx <- which(Study == j) # rows from this study
        
        #divide by weights to make objective consider sample sizes
        obj <- obj + (1 - eta) / 2 * u3[j] * sum( ( data[indx, 2] -  beta[1, j] - data[indx,-c(1,2)] %*% beta[-1, j] )^2     )
        
        
        #     objLinear(X = as.matrix( cbind(1, data[indx,-c(1,2)] ) ), 
        #                                            Y = as.vector( data[indx, 2] ), 
        #                                            beta = as.matrix(beta[,j],
        #                                                             weights = 1
        #                                            )
        # )
        
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
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
        u2 <- 1 # LHS Objective weight
        
        # RHS weights
        u3 <- 1 / nVec # K x 1 vector
        u4 <- rep(1, K)
        
      
    }else if(sampSzWeight == 7){
        # variance adjustment on both sides
        # LHS weights
        u1 <- 1 / N # 1 / (rep(sigStack, nVec) * N ) # LHS Objective weight   - N x 1 vector
        # weights <- diag( rep(sigStack, nVec) ) # use weights so its easier for matrix multiplication instead of using u1 scalar multiplication
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

objOEC_old <- function(data, 
                   beta, 
                   w, 
                   mu = NULL, 
                   lambdaVec = NULL, 
                   stackInt = TRUE, 
                   eta = 0.5,
                   StackIndx = NULL,
                   sampSzWeight = FALSE){
    # returns objective of OEC Linear
    # sampSzWeight determines whether obejective weights sample sizes on RHS
    
    Study <- data$Study
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each study
    data <- as.matrix(data)
    K <- length(unique(Study))
        
    
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
    
    # if intercept in stacking regression add column of ones to beta
    if(stackInt){
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / (2 * N) * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ), 
                                   Y = as.matrix( data[SIndx,2] ), 
                                   beta = as.matrix( cbind(1, beta) ) %*% as.vector(w) )
    }else{
        
        # stacking function -- divide by N to make objective consistent with glmnet
        obj <- eta / (2 * N) * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ), 
                                   Y = as.matrix( data[SIndx,2] ), 
                                   beta = as.matrix( beta ) %*% as.vector(w) ) 
    }
    
    if(sampSzWeight){
        for(j in 1:K){
            
            indx <- which(Study == j) # rows from this study
            sslRows <- intersect( SSLIndx, indx ) # rows from this study that are used for SSL training
            
            #divide by N to make objective consider sample sizes
            obj <- obj + (1 - eta) / (2 * N) * objLinear(X = as.matrix( cbind(1, data[sslRows,-c(1,2)] ) ), 
                                                         Y = as.vector( data[sslRows, 2] ), 
                                                         beta = as.matrix(beta[,j])
            )
            
        }
        
        lambdaVec <- lambdaVec * nVec / N # scale L2 penalty on study specific betas by sample sizes
        # add penalties
        obj <- obj + (1 - eta) * sum( beta^2 %*% lambdaVec ) / 2 + eta * mu * sum( w^2 ) / 2 # the 2s are to make objective consistent with glmnet 
        
        
    }else{
        # original way
        for(j in 1:K){
            
            indx <- which(Study == j) # rows from this study
            sslRows <- intersect( SSLIndx, indx ) # rows from this study that are used for SSL training
            
            #divide by n_k to make objective consistent with glmnet
            obj <- obj + (1 - eta) / (2 * nVec[j]) * objLinear(X = as.matrix( cbind(1, data[sslRows,-c(1,2)] ) ), 
                                                         Y = as.vector( data[sslRows, 2] ), 
                                                         beta = as.matrix(beta[,j])
            )
            
        }
        
        # add penalties
        obj <- obj + (1 - eta) * sum( beta^2 %*% lambdaVec ) / 2 + eta * mu * sum( w^2 ) / 2 # the 2s are to make objective consistent with glmnet 
        
    }
    
    return(obj)
    
}



# used objOEC_unscaled until conversation with Rahul on May 22, 2020 -- made consistent with GLMNET
objOEC_unscaled <- function(data, 
                   beta, 
                   w, 
                   mu = NULL, 
                   lambdaVec = NULL, 
                   stackInt = TRUE, 
                   eta = 0.5,
                   StackIndx = NULL){
    # returns objective of OEC Linear
    
    Study <- data$Study
    data <- as.matrix(data)
    K <- length(unique(Study))
    N <- nrow(data) # sample size of Merged
    nVec <- as.vector( table(data$Study) ) # vector of sample sizes of each
    
    
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
    
    # if intercept in stacking regression add column of ones to beta
    if(stackInt){
        
        # stacking function
        obj <- eta * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ), 
                         Y = as.matrix( data[SIndx,2] ), 
                         beta = as.matrix( cbind(1, beta) ) %*% as.vector(w) )
    }else{
        
        # stacking function
        obj <- eta * objLinear(X = as.matrix( cbind(1, data[SIndx,-c(1,2)]) ), 
                         Y = as.matrix( data[SIndx,2] ), 
                         beta = as.matrix( beta ) %*% as.vector(w) ) 
    }
    
    
    for(j in 1:K){
        
        indx <- which(Study == j) # rows from this study
        sslRows <- intersect( SSLIndx, indx ) # rows from this study that are used for SSL training
        
        obj <- obj + (1 - eta) * objLinear(X = as.matrix( cbind(1, data[sslRows,-c(1,2)] ) ), 
                               Y = as.vector( data[sslRows, 2] ), 
                               beta = as.matrix(beta[,j])
        )
        
    }
    
    # add penalties
    obj <- obj + (1 - eta) * sum( beta^2 %*% lambdaVec ) + eta * mu * sum( w^2 )
    
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
        # if(is.null(coefWeights)){
        #     diagMat <- as.matrix( diag( ncol(x) ) )
        #     diagMat[1,1] <- 0 # do not penalize intercept
        # }
        
    }else{
        x <- as.matrix(x)
        # if(is.null(coefWeights)){
        #     diagMat <- as.matrix( diag( ncol(x) ) )
        # }
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
        #wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
        
        for(tune in 1:nrow(tune.grid) ){
            
            # print(paste0("study ", study, "_tune_", tune))
            
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
            tuneParam <- sslLambdas #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
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
        ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies
        
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
                
                # print(paste0("study ", study, "_tune_", tune))
                
                stackParam <- warmRes$mu
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
            tuneParam <- sslLambdas[-study] #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
            
            for(tune in 1:length(tune.grid) ){

                # print(paste0("study ", study, "_tune_", tune))
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
                                      ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies
        
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
            ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies
            
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
            ) #createFolds(factor(data$Study), k = nfolds) # make folds with balanced studies
            
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
                
                # print(paste0("study ", study, "_tune_", tune))
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
            tuneParam <- sslLambdas[-study] #rep( tune.grid[tune, 2],  length(unique(data$Study[indxList]) )  ) # use the same parameter for all studies
            
            for(tune in 1:length(tune.grid) ){
                
                # print(paste0("study ", study, "_tune_", tune))
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


# 
# 
# logisticWOpt <- function(dat, beta, w_stack, lambda, solverName = "ECOS"){
#     # beta is K x p matrix
#     # dat is dataframe with following columns: Study   Y   X_1   X_2  ....  X_p
#     # lambda is tradeoff between stacking and individual loss functions
#     K <- length(unique(dat$Study))
#     #print(K)
#     
#     # beta <- as.matrix(beta)
#     # dat <- as.matrix(dat)
#     # w <- as.vector(w)
#     
#     # number of studies
#     
#     len <- length(w_stack)
#     #print(len)
#     
#     if( is.null(solverName) ){
#         # if to do with cvxr
#         
#         # variable to optimize
#         w <- Variable( len ) # make matrix of variables K x p matrix (p includes intercept)
#         oneVec0 <- rep(1, nrow(dat[dat$Y <= 0, ])) # just rows with outcome = 0
#         oneVec <- rep(1, nrow(dat)) # all rows
#         
#         if(len == K){
#             # print("1")
#             # if no intercept in stacking regression
#             obj <- lambda * ( -sum( as.matrix(dat[dat$Y <= 0, -c(1,2)]) %*% t(as.matrix(beta)) %*% (w) ) - 
#                                   sum( logistic(- as.matrix(dat[,-c(1,2)]) %*% t(as.matrix(beta)) %*% (w)) ) )  # w has intercept
#         }else if(len == K + 1){
#             # print("2")
#             # print(dim(t(beta)))
#             # print(dim(dat[,-c(1,2)]))
#             # if there is an intercept then length(w) == K + 1
#             obj <- lambda * ( -sum( as.matrix( dat[dat$Y <= 0, -c(1,2)] ) %*% t(as.matrix(beta)) %*% (w[2:len]) + oneVec0 * w[1]) - 
#                                   sum(logistic(as.matrix( -dat[,-c(1,2)] ) %*% t(as.matrix(beta)) %*% (w[2:len]) + oneVec * w[1]) ) )  # w has intercept
#         }else{
#              print("dim error")
#         }
#         
#         # oneVec0 *
#         # oneVec *
#         
#         prob <- Problem(Maximize(obj)) # maximize likelihood
#         result <- solve( prob, solver = solverName )
#         #print(result$status)
#         w_res <- result$getValue(w)
#     }else{
#         # else do with PGD
#        # print("NULL NO SOLVE")
#         
#     }
#     
#     return(w_res)
# }

