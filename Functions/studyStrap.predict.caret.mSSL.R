#' Study Strap Prediction Function: Makes predictions on object of class "ss"
#'
#' @param ss.obj A model object (of class "ss") fit with studyStrap package (e.g., ss, cmss, sse, merge).
#' @param X A dataframe of the study to make predictions on. Must include covariates with the same names as those used to train models.
#' @return Matrix of predictions. Each column are predictions with different weighting schemes.
#' @examples
#' ##########################
#' ##### Simulate Data ######
#' ##########################
#'
#' set.seed(1)
#' # create half of training dataset from 1 distribution
#' X1 <- matrix(rnorm(2000), ncol = 2) # design matrix - 2 covariates
#' B1 <- c(5, 10, 15) # true beta coefficients
#' y1 <- cbind(1, X1) %*% B1
#'
#' # create 2nd half of training dataset from another distribution
#' X2 <- matrix(rnorm(2000, 1,2), ncol = 2) # design matrix - 2 covariates
#' B2 <- c(10, 5, 0) # true beta coefficients
#' y2 <- cbind(1, X2) %*% B2
#'
#' X <- rbind(X1, X2)
#' y <- c(y1, y2)
#'
#' study <- sample.int(10, 2000, replace = TRUE) # 10 studies
#' data <- data.frame( Study = study, Y = y, V1 = X[,1], V2 = X[,2] )
#'
#' # create target study design matrix for covariate profile similarity weighting and
#' # accept/reject algorithm (covaraite-matched study strap)
#' target <- matrix(rnorm(1000, 3, 5), ncol = 2) # design matrix
#' colnames(target) <- c("V1", "V2")
#'
#' ##########################
#' ##### Model Fitting #####
#' ##########################
#'
#' # Fit model with 1 Single-Study Learner (SSL): PCA Regression
#' ssMod1 <- ss(data = data, formula = Y ~.,
#'             target.study = target,
#'             bag.size = length(unique(data$Study)), straps = 5, stack = "standard",
#'             sim.covs = NA, ssl.method = list("pcr"),
#'             ssl.tuneGrid = list(data.frame("ncomp" = 2)),
#'             sim.mets = TRUE,
#'             model = TRUE, customFNs = list() )
#'
#' #########################
#' #####  Predictions ######
#' #########################
#'
#' preds <- studyStrap.predict(ssMod1, target)
#' @importFrom stats predict
#' @export

studyStrap.predict <- function(ss.obj, X){
    #ss.obj is the object
    # X is the design matrix and must have the same covariates

    num.SSLs <- length(ss.obj$modelInfo$SSL)


    if( is.matrix(ss.obj$simMat) ){
        # check to see if similarity metrics were calculated
        preds.mat <- matrix(nrow = nrow(X), ncol = ncol(ss.obj$simMat) + 2)
        # 23 similarity measures and average and stacking weights
        colnames(preds.mat) <- c("Avg", paste0(ss.obj$modelInfo$stack.type, "_Stacking"),
                                 colnames(ss.obj$simMat) )

    }else{
        preds.mat <- matrix(nrow = nrow(X), ncol = 2) # just avg and stacking weights
        colnames(preds.mat) <- c("Avg", paste0(ss.obj$modelInfo$stack.type, "_Stacking") )
    }

    raw.preds <- matrix(nrow = nrow(X),
                        ncol = ss.obj$modelInfo$numStraps * num.SSLs ) # raw predictions

    counter <- 1
    for(SSL in 1:num.SSLs){
        # iterate through SSL methods
        for(mod in 1:ss.obj$modelInfo$numStraps){
            # iterate through models within each SSL method

            raw.preds[,counter] <- predict(ss.obj$models[[SSL]][[mod]], X)  # make predictions on each model
            counter <- counter + 1
            }
    }

    if( is.matrix(ss.obj$simMat) ){
        if( num.SSLs > 1){
            # if multiple SSLs need to make sim.matrix into weights by replicating and scaling appropriately
            ss.obj$simMat <- do.call(rbind, replicate(num.SSLs, ss.obj$simMat, simplify=FALSE))
            ss.obj$simMat <- prop.table( ss.obj$simMat, 2 ) # ensure weights sum to 1
            # rbind copies of the sim.matrix and scale by number of SSLs so add to 1
        }
    }

        preds.mat[,1] <- rowMeans(raw.preds) # simple average
        preds.mat[,2] <- cbind(1, raw.preds) %*% ss.obj$stack.coefs # stacking includes intercept

        if( is.matrix(ss.obj$simMat) ){
            preds.mat[, 3:ncol(preds.mat) ] <- raw.preds %*% ss.obj$simMat # use similarity weights
        }

    return(preds.mat)
}
