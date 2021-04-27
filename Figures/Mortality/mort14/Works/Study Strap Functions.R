#------------------------------------------------------
# Wiggle Function 
#------------------------------------------------------
# this variation only calculates the distance between the 
wiggle.fn <- function(base.CV, win = 2, wiggle = 10){
  #train.mat -- the matrix that forms the base points (can be the same as test.mat)
  # test.mat -- what is analyzed (can be the same as train.mat)
  # win      -- the window around the points that are included in the output in addition to the base points/wiggle points
  # wiggle   -- the window around which inflection points are detected
  # takes in a matrix of CVs: n x 1000
  # averages the CVs to use as a base CV for comparison
  # generates a window of 5 points to the left and right 
  # and if a new CVs inflection points fall in that range it uses either those inflection points
  # or the inflection points of the base CV otherwise use 0
  # 
  ####### removed to save memory/increase speed
  # train.mat <- as.matrix(train.mat) #necessary for diff function
  # test.mat <- as.matrix(test.mat) #necessary for diff function
  # base.CV <- colMeans(train.mat) #average CV
  ####### 
  
  infl <- c(FALSE, diff(diff(base.CV)>0)!=0) #inflection points of average CV ---ORIGINAL WAY
  #infl <- diff(sign(diff(base.CV)))
  #base.points <- which(infl == -2) #inflection points of average CV
  base.points <- which(infl == TRUE) #inflection points of average CV
  #-------------------------
  # Eliminate similar points within window of 30 # added this in
  #-------------------------
  inf.points <- c()
  
  for(y in 1:length(base.points)){
    if (any(abs(base.points[y] - base.points) <= 30)){
      if (any(abs(base.points[y] - inf.points) <= 30)){ #make sure there arent already these values in the new list
      }
      else{
        inf.points <- c(inf.points, mean(base.points[which(abs(base.points[y] - base.points) <= 30)], base.points[y])) #if any are within a window of 30 take the mean
      }
    }
    else{
      inf.points <- c(inf.points, base.points[y])
    }
  }
  base.points <- inf.points #change it back for consistency below
  
  return(c(base.points, base.CV[base.points]))
}


sim.metrics <- function(dat1, dat2){
    # currently 21 metrics
    
    # Takes in 2 design matrices (outcome vector is not attached) to compare and returns vector of 
    # similarity metrics that when combined for each iteration that can be used to calculate weights
    # 
    
    # best seem to be 
    # sum of covariance matrix of CCA
    # mean correlation between corresponding variables
    # r1 and r3
    similarity.vec <- c()
    
    
    # number of covariates
    covs <- ncol(dat1)
    n <- nrow(dat1)
    p <- ncol(dat1)
    q <- ncol(dat2)
    
    ############################################
    # Correlation Coefficient Similarity Metric
    ############################################
    # https://stats.idre.ucla.edu/r/dae/canonical-correlation-analysis/
    ## UCLA
    
    # correlation matrix (make sure rows match up)
    mc1 <- matcor(dat1[1:min(nrow(dat1), nrow(dat2)),], dat2[1:min(nrow(dat1), nrow(dat2)),]) # correlation matrix
    
    #mc1[[3]][1:p,(q+1):(p+q)] == t(mc1[[3]][(p+1):(p+q),1:q]) # TRUE
    cfs <- diag( mc1[[3]][1:p, (q+1):(p+q)] ) # correlation coefficients between corresponding values
    
    # similarity metric based on the mean correlation coefficients between predictors (assuming same predictors in test (target) and training sets)
    #similarity.metric.corr <- mean(cfs) # has negative values
    similarity.metric.corr.abs <- mean(abs(cfs)) # the higher the better
    similarity.metric.corr.sq <- sum(mc1[[3]][1:p, (q+1):(p+q)]) # the higher the better
    similarity.metric.corr.sq.abs <- sum(abs(mc1[[3]][1:p, (q+1):(p+q)]))
    # add to vector of similarity metrics
    similarity.vec <- c(similarity.vec, similarity.metric.corr.abs, 
                        similarity.metric.corr.sq, similarity.metric.corr.sq.abs)
    # IDEA -- could weight the correlation coefficients of variables based on importance of variables (VIM or B coefficients)
    # when coming up with a mean similarity metric (as opposed to just an average)
    
    
    ############################################
    # Canonical Correlations Similarity Metric
    ############################################
    # can do this way or way below
    
    #cc1 <- cc(dat1[1:min(nrow(dat1), nrow(dat2)),], dat2[1:min(nrow(dat1), nrow(dat2)),])
    # canonical correlations
    # similarity.metric.can <- mean(cc1$cor)
    # similarity.metric.can.abs <- mean(abs(cc1$cor))
    # similarity.metric.can.sq <- mean((cc1$cor)^2)
    # 
    # # add to vector of similarity metrics
    # similarity.vec <- c(similarity.vec, similarity.metric.can, similarity.metric.can.abs, 
    #                     similarity.metric.can.sq)
    
    ###################################################
    # Canonical Correlations Loadings Similarity Metric
    ###################################################
    # ******************maybe use canonical loadings in situations where different covariates are present******************
    
    
    # compute canonical loadings
    #cc2 <- comput(dat1, dat2, cc1)
    
    # display canonical loadings
    #cc2[3:6]
    
    ############################################
    # Covariance of Canonical Variables
    ############################################
    # http://users.stat.umn.edu/~helwig/notes/cancor-Notes.pdf
    # Nathaniel E. Helwig
    # Univ of Minnesota
    
    # standardize data
    Xs <- scale(dat1[1:min(nrow(dat1), nrow(dat2)),])
    Ys <- scale(dat2[1:min(nrow(dat1), nrow(dat2)),])
    
    # cca (the normal way)
    Sx <- cov(Xs)
    Sy <- cov(Ys)
    Sxy <- cov(Xs,Ys)
    Sxeig <- eigen(Sx, symmetric=TRUE)
    Sxisqrt <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
    Syeig <- eigen(Sy, symmetric=TRUE)
    Syisqrt <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
    Xmat <- Sxisqrt %*% Sxy %*% solve(Sy) %*% t(Sxy) %*% Sxisqrt
    Ymat <- Syisqrt %*% t(Sxy) %*% solve(Sx) %*% Sxy %*% Syisqrt
    Xeig <- eigen(Xmat, symmetric=TRUE)
    Yeig <- eigen(Ymat, symmetric=TRUE)
    
    
    #-------------------------------------------
    # Canonical Correlations Similarity Metric
    #-------------------------------------------
    # can also use cc1 <- cc(dat1[1:min(nrow(dat1), nrow(dat2)),], dat2[1:min(nrow(dat1), nrow(dat2)),])
    # using this way so do not have to recalculate twice
    rho <- sqrt(Xeig$values)
    
    similarity.metric.can <- mean(rho)
    similarity.metric.can.abs <- mean(abs(rho))
    similarity.metric.can.sq <- mean((rho)^2)
    
    # add to vector of similarity metrics
    similarity.vec <- c(similarity.vec, similarity.metric.can, 
                        similarity.metric.can.sq)
    #-------------------------------------------
    
    #-------------------------------------------
    # Covariance of Canonical Variables
    #-------------------------------------------
    
    #
    Ahat <- Sxisqrt %*% Xeig$vectors
    Bhat <- Syisqrt %*% Yeig$vectors
    Ainv <- solve(Ahat)
    Binv <- solve(Bhat)
    
    # canonical variables
    U <- Xs %*% Ahat
    V <- Ys %*% Bhat
    
    # covariance of canonical variables (U and V)
    rhomat <- cbind(diag(rho), matrix(0, p, q-p))
    similarity.metric.can.UV.rho <- sum( (cov(U, V) - rhomat)^2 )
    similarity.metric.can.UV.cov <- sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv)^2 )
    sim4 <- sum( (cov(U, V) - rhomat) )
    sim5 <- sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv) ) #So far the best 
    sim6 <- sum(diag((cov(U, V) - rhomat)^2))
    sim7 <- sum(diag((Sxy - crossprod(Ainv, rhomat) %*% Binv)))
    sim8 <- sum(diag((Sxy - crossprod(Ainv, rhomat) %*% Binv))^2)
    
    mean.cor <- abs(cor(colMeans(dat1), colMeans(dat2)))
    
    similarity.vec <- c(similarity.vec, similarity.metric.can.UV.rho, similarity.metric.can.UV.cov, sim4, sim5, sim6, sim7, sim8, mean.cor)
    
    # similarity.vec <- c(similarity.vec, similarity.metric.can.UV.rho, similarity.metric.can.UV.cov, sim4, sim5, sim6, )
    
    
    #-------------------------------------------
    # Matrix Correlations
    #-------------------------------------------
    
    # added from MatrixCorrelations package
    ## number of components is arbitrary--just chose number of columns of each matrix to try
    comp1 <- p
    comp2 <- q
    #similarity.vec <- c(similarity.vec, (allCorrelations(Xs, Ys, ncomp1 = comp1, ncomp2 = comp2))^2 )
    similarity.vec <- c(similarity.vec, (allCorrelations(Xs, Ys, ncomp1 = comp1, ncomp2 = comp2)) )
    
    return(similarity.vec)
}

# Trim fat
fatTrim = function(cm) {
  # modified from http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$scores = c()  
  cm$loadings = c()
  cm$weights = c()
  cm$Yloadings = c()
  cm$Xtotvar = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}



#------------------------------------------------------
# Wiggle Room 
#------------------------------------------------------
# this variation has no base 
point.finder <- function(train.mat, test.mat, win = 2, wiggle = 10){
  #train.mat -- the matrix that forms the base points (can be the same as test.mat)
  # test.mat -- what is analyzed (can be the same as train.mat)
  # win      -- the window around the points that are included in the output in addition to the base points/wiggle points
  # wiggle   -- the window around which inflection points are detected
  # takes in a matrix of CVs: n x 1000
  # averages the CVs to use as a base CV for comparison
  # generates a window of 5 points to the left and right 
  # and if a new CVs inflection points fall in that range it uses either those inflection points
  # or the inflection points of the base CV otherwise use 0
  # 
  train.mat <- as.matrix(train.mat) #necessary for diff function
  test.mat <- as.matrix(test.mat) #necessary for diff function
  base.CV <- colMeans(train.mat) #average CV
  infl <- c(FALSE, diff(diff(base.CV)>0)!=0) #inflection points of average CV ---ORIGINAL WAY
  #infl <- diff(sign(diff(base.CV)))
  #base.points <- which(infl == -2) #inflection points of average CV
  base.points <- which(infl == TRUE) #inflection points of average CV
  #-------------------------
  # Eliminate similar points within window of 30 # added this in
  #-------------------------
  inf.points <- c()
  
  for(y in 1:length(base.points)){
    if (any(abs(base.points[y] - base.points) <= 30)){
      if (any(abs(base.points[y] - inf.points) <= 30)){ #make sure there arent already these values in the new list
      }
      else{
        inf.points <- c(inf.points, mean(base.points[which(abs(base.points[y] - base.points) <= 30)], base.points[y])) #if any are within a window of 30 take the mean
      }
    }
    else{
      inf.points <- c(inf.points, base.points[y])
    }
  }
  base.points <- inf.points #change it back for consistency below
  #-------------------------
  
  window0 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- if a point falls in a window,  use the current at that point, and 0 if it falls outside the window
  window.base <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- if a point falls in a window, and use the current at the base point if it falls outside the window
  base0 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- use the base point if any inflection point falls in a base point window, 0s otherwise
  base.base <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points)) #- use base points regardless
  
  #gives window of 10 voltage points around modes
  window10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- if a point falls in a window,  use the current at that point, and 0 if it falls outside the window
  window.base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- if a point falls in a window, and use the current at the base point if it falls outside the window
  base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- use the base point if any inflection point falls in a base point window, 0s otherwise
  base.base10 <- matrix(NA, nrow = nrow(test.mat), ncol = length(base.points) * (win * 2 + 1)) #- use base points regardless
  
  
  window.points <- c()
  for (x in 1:length(base.points)){
    window.points <- c(window.points, (base.points[x] - wiggle):(base.points[x] + wiggle)) #create window of points where
    #inflection points can fall
  }
  d <- 0
  skips <- 0 # counts the number of times there are no inflection points in the curve and use previous CV's inflection points
  infl.vec <-c()
  for (i in 1:nrow(test.mat)){ #NEED TO FIGURE OUT WHAT TO DO IF INFL is SHORTER THAN BASE.POINTS
    #print(i)
    ## volt.mat[i,] <- voltammograms.m[seq(1,17001, by = 1000)[i]:seq(1000,18000, by = 1000)[i],1]
    #infl <- c(FALSE, diff(diff(test.mat[i,])>0)!=0) #inflection points for each CV ---ORIGINAL WAY
    #infl <- which(infl == TRUE) #add inflection point point numbers (voltage potential number) to vector ---ORIGINAL WAY
    #if (sum(diff(sign(diff(test.mat[i,]))) == -2) == 0){
    if (length(which(c(FALSE, diff(diff(test.mat[i,])>0)!=0)==TRUE)) == 0){
      skips <- skips + 1 #if there are no inflection points then use the inflection points from the previous CV
    }
    else{
      infl <- c(FALSE, diff(diff(test.mat[i,])>0)!=0) #inflection points for each CV ---ORIGINAL WAY
      infl <- which(infl == TRUE) #add inflection point point numbers (voltage potential number) to vector ---ORIGINAL WAY
      #-------------------------
      # Eliminate similar points within window of 30 # added this in
      #-------------------------
      inf.points <- c()
      
      for(y in 1:length(infl)){
        if (any(abs(infl[y] - infl) <= 30)){
          if (any(abs(infl[y] - inf.points) <= 30)){ #make sure there arent already these values in the new list
          }
          else{
            inf.points <- c(inf.points, mean(infl[which(abs(infl[y] - infl) <= 30)], infl[y])) #if any are within a window of 30 take the mean
          }
        }
        else{
          inf.points <- c(inf.points, infl[y])
        }
      }
      infl <- inf.points #change it back for consistency below
      
      
      #infl <- c(FALSE, diff(sign(diff(test.mat[i,]))))
      #infl <- which(infl == -2)
    }
    
    #points.keep <- c()
    #print(i)
    
    #add conditional to determine if base.points and infl are of different lengths and if they arent then use either 0 or basepoints if its shorter (maybe use furthest point if its basepoints)
    w10.vec <- c()
    wb10.vec <- c()
    b10.vec <- c()
    bb10.vec <- c()
    for (x in 1:length(base.points)){
      #print(x)
      if(abs(infl[which.min(abs(infl - base.points[x]))] - base.points[x]) <= wiggle){ # if the biggest difference between a point another point is <=10, keep it
        #print(1)
        #points.keep <- c(points.keep, infl[which.min(abs(infl - y))])
        window0[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        w10.vec <- c(w10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)]) #adds values of a window around the inflection point
        # w10.vec <- c(w10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x])) - window):(infl[which.min(abs(infl - base.points[x]))] + window)])
        
        # print(2)
        window.base[i,x] <- test.mat[i,infl[which.min(abs(infl - base.points[x]))]] #if the point fits then use the inflection point specific to each CV
        wb10.vec <- c(wb10.vec, test.mat[i,(infl[which.min(abs(infl - base.points[x]))] - win):(infl[which.min(abs(infl - base.points[x]))] + win)])
        
        # print(3)
        # base0[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point fits, then use the closest base point
        # base.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]]
        base0[i,x] <- test.mat[i, base.points[x]] #if the point fits, then use the closest base point
        b10.vec <- c(b10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        # print(4)
        base.base[i,x] <- test.mat[i, base.points[x]] #either way use base inflection points
        bb10.vec <- c(bb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
      }
      
      else{
        #print("else")
        d <- d + 1
        window0[i,x] <- 0 #0s if inflection points don't fall into window
        w10.vec <- c(w10.vec, rep(0, win * 2 + 1))
        
        #window.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point doesn't fit then use the base inflection points
        window.base[i,x] <- test.mat[i, base.points[x]] #if the point doesn't fit then use the base inflection points
        wb10.vec <- c(wb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
        
        base0[i,x] <- 0 #0s if inflection points don't fall into window, use 0
        b10.vec <- c(b10.vec, rep(0,win * 2 + 1))
        
        # base.base[i,x] <- test.mat[i, base.points[which.min(base.points - infl[x])]] #if the point fits then use the base inflection points
        base.base[i,x] <- test.mat[i, base.points[x]] #either way use the base inflection points
        bb10.vec <- c(bb10.vec, test.mat[i, (base.points[x] - win):(base.points[x] + win)])
        
      }
    }
    
    # window points added to matrix after done iterating through base points
    window10[i, ] <- w10.vec
    window.base10[i,] <- wb10.vec
    base10[i,] <- b10.vec
    base.base10[i,] <- bb10.vec
  }
  
  print("Point counter")
  print(d)
  print(paste0("skip counter: ", skips))
  return(list(window0, window.base, base0, base.base, window10, window.base10, base10, base.base10, c(base.points, base.CV[base.points]), base.CV))
}
print("wiggle loaded")
#############################################################

####################################
# Bag Generator Function
####################################

bagGen <- function(full.CVs, target.CVs, test.elec, bag.size = nrow(full.CVs) - 1){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  
  target.CVs <- as.matrix(target.CVs)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs) # number of studies not including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  A <- A[,-test.elec] # remove   
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = nrow(target.CVs)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = nrow(target.CVs)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("Prop_E",study.num)
  row.names(bag.recipe) <- paste0("Target",1:nrow(target.CVs))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  # create w vector to estimate (i.e., weights to learn)
  # for (x in 1:bag.size){
  #   assign(paste0("y", x), Int(1))
  # }
  # # need to change this for more generalizable bag.size
  # # for(x in 2:bag.size){
  # #   w <- vstack(w, get(paste0("y", x)))
  # # }
  # w <- vstack(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10,y11, y12, y13, y14) ## Create x expression
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  for (i in 1:nrow(bag.mat)){
    elec.target <- target.CVs[i,] # set coordinates of target electrode

    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 

    # Constrained Optimization (Integer Least Squares)
    
    problem <- Problem(objective, constraints = list(w >= 0, 
                                                     sum(w[1:bag.size]) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
      bag <- c(bag, rep(x, w.solution[x]))
    }
    #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
    # so need to use list of elecs (study num) and index them by bag
    bag.mat[i,] <- study.num[bag] # add bag to matrix.
    bag.recipe[i, ] <- w.solution  # add bag recipe to 
  }
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}



####################################
# Bag Generator Function Based on distances
####################################

bagGen.dist <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  
  distances <- as.vector(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:length(distances))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
 
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) + lambda * sum(w^2) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
      }
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.w_Norm2 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) + lambda * norm2(w) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}

bagGen.w_NormInf <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) + lambda * norm_inf(w) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.dist <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) + lambda * sum(w^2) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


######################################################
# Generate Target Electrode Inflection Points
######################################################
target.fn <- function(full.CVs, test.elec, straps, dist = 1, sprd = - 2.5, volts = 2){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # dist is a scalar indicating the method of selecting the distances from the test electrode 
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  # volts = 1 is calculating voltages of target CVs based upon a linear model
  # volts = 2 is uses voltages of test elec for all target CVs 
  # dist = 1 is evenly spaced intervals for all inflection points
  # dist = 2 is densely spaced intervals using the "sprd" variable to specify the degree of density
  # dist = 3 is using just the distance metric 
  #------------------------------------------
  
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  ##################
  # fits
  ##################
  A.1 <- data.frame(t(A[c(1,5),]))
  colnames(A.1) <- c("y","x")
  lm.1 <- lm(y ~ x, data = A.1)
  
  A.2 <- data.frame(t(A[c(2,6),]))
  colnames(A.2) <- c("y","x")
  lm.2 <- lm(y ~ x, data = A.2)
  
  A.3 <- data.frame(t(A[c(3,7),]))
  colnames(A.3) <- c("y","x")
  lm.3 <- lm(y ~ x, data = A.3)
  
  A.4 <- data.frame(t(A[c(4,8),]))
  colnames(A.4) <- c("y","x")
  lm.4 <- lm(y ~ x, data = A.4)
  
  ###########################################   
  # CV boundaries (heightest and lowest)
  ########################################### 
  # 4 inflection point boundaries
  CV.boundaries <- matrix(NA, ncol = 2, nrow = 4) 
  colnames(CV.boundaries) <- c("Min", "Max")
  row.names(CV.boundaries) <- paste0("Inflect_", 1:4)
  for (i in 1:3){
    CV.boundaries[i, 1] <- min(A[4 + i, ]) # minimum height (current)
    CV.boundaries[i, 2] <- max(A[4 + i, ]) # maximum height (current)
  }
  # reverse for last metric because it is negative
  CV.boundaries[4, 1] <- max(A[8, ]) # minimum height (current)
  CV.boundaries[4, 2] <- min(A[8, ]) # minimum height (current)
  #---------------------------------------
  
  # find min/max elecs associated with first hump
  
  source("Study Strap Functions.R")
  sim.mat <- matrix(NA, ncol = 3, nrow = K) # matrix of similarities to test elec and Heights
  colnames(sim.mat) <- c("Elec", "Similarity", "Height")
  row.names(sim.mat) <- paste0("Elec_", 1:K)
  sim.mat[, 1] <- 1:K # elec numbers
  s.test <- wiggle.fn(full.CVs[test.elec,-c(1,2)]) # inflection points of test elec
  for (i in seq(1,K)){
    s <- wiggle.fn(full.CVs[i, -c(1,2)]) # inflection points of each elec
    sim.mat[i, 2] <- sum( (s - s.test)^2  ) # L2 norm of each elec to test elec
    sim.mat[i, 3] <- sum( (A[5:8,i])^2  ) # Squared "heights" of inflection points
  }
  
  # highest and lowest CVs (i.e., those on the extremes)
  height.max <- which.max(sim.mat[,3])
  height.min <- which.min(sim.mat[,3])
  
  ###############################################################
  # Generate Sequences (currents of target electrodes)
  ###############################################################
  # straps <- 15 # number of straps for sequence
  seq.mat <- matrix(NA, nrow = 4, ncol = straps) # matrix of currents to make predictions on
  colnames(seq.mat) <- paste0("I_Strap_", 1:straps)
  row.names(seq.mat) <- paste0("Inflect_", 1:4)
  
  # dist is different methods of generating target electrode distances
  if (dist ==  1){
      # first check to see if test elec is on extremes
      if (height.max == test.elec | height.min == test.elec){
        # use CV.boundaries as boundaries for sequence of points 
        for (i in 1:4){
          # make a sequence from the lowest possible inflection point to the highest for each inflection point
          seq.mat[i,] <- seq(CV.boundaries[i, 1], CV.boundaries[i, 2], length = straps)
        }
        # A[4 + i, test.elec]
      }else{
        # if test elec is not on boundary (is in between training elecs), then make sequence with
        # number of training elecs above and below
        
        # number of training elecs above and below based on sum of squared heights of all 
        # inflection points
        num.above <- sum(sim.mat[test.elec, 3] <  sim.mat[-test.elec, 3] )
        num.below <- sum(sim.mat[test.elec, 3] >  sim.mat[-test.elec, 3] )
        
        # calculate corresponding amount of straps proportional to number of training
        # straps above and below actual strap
        straps.above <- round( straps * num.above / (num.above + num.below) )
        straps.below <- straps - straps.above
        
        offset <- 0 # amount of current to offset around test.elec
        
        for (i in 1:4){
          #FIX THIS
          # make a sequence from the lowest possible inflection point to the test electrode 
          # of length equal to the number of straps calculated for "Straps below"
          seq.mat[i,1:straps.below] <- seq(from = CV.boundaries[i, 1], 
                                           to = A[4 + i, test.elec], length = straps.below)
          
          # make a sequence from the test electrode (plus an offset so we dont end up 
          # with duplicates) to the highest possible inflection point 
          # of length equal to the number of straps calculated for "Straps above"
          
          seq.mat[i,(straps.below + 1):straps] <- seq( from = (A[4 + i, test.elec] + offset), 
                                                       to = CV.boundaries[i, 2], 
                                                       length = straps.above)
        }
        
        
        # determine whether more space above or below
        # determine proportions of straps above or below based on that
        # generate sequence. Make sure for the evenly spaced to give a little distance from test elec on either side so get close elctrodes on both sides
      }
  }else if (dist == 2){
        # first check to see if test elec is on extremes
        if (height.min == test.elec){
            # if test elec is the lowest
            # use CV.boundaries as boundaries for sequence of points 
            d <- 10^(seq(sprd, 0, length = straps))
            d[length(d)] <- 0 # make the largest value (the most distal strap) equal to 0 (i.e., the target electrode's inflection points)
            for (i in 1:4){
                dist.diff <-  CV.boundaries[i, 2] - CV.boundaries[i, 1] # distance between test.elec and boundary
                # make a sequence from the lowest possible inflection point to the highest for each inflection point
                seq.mat[i,] <- CV.boundaries[i, 1] + dist.diff * d # from the minimum (which is also the test elec to the max)
            }
        }else if (height.max == test.elec ){
          # if test elec is the highest
          # use CV.boundaries as boundaries for sequence of points
            d <- 10^(seq(sprd, 0, length = straps))
            d[length(d)] <- 0 # make the largest value (the most distal strap) equal to 0 (i.e., the target electrode's inflection points)
            for (i in 1:4){
                dist.diff <-  CV.boundaries[i, 2] - CV.boundaries[i, 1] # distance between test.elec and boundary
                # make a sequence from the lowest possible inflection point to the highest for each inflection point
                seq.mat[i,] <- CV.boundaries[i, 2] - dist.diff * d  # from the max (which is also the test elec to the min)
                }
          # A[4 + i, test.elec]
        }else{ 
          # if test elec is not on boundary (is in between training elecs), then make sequence with
          # number of training elecs above and below
          
          # number of training elecs above and below based on sum of squared heights of all 
          # inflection points
          num.above <- sum(sim.mat[test.elec, 3] <  sim.mat[-test.elec, 3] )
          num.below <- sum(sim.mat[test.elec, 3] >=  sim.mat[-test.elec, 3] )
          
          # calculate corresponding amount of straps proportional to number of training
          # straps above and below actual strap
          straps.above <- round( straps * num.above / (num.above + num.below) )
          straps.below <- straps - straps.above

          offset <- 0 # amount of current to offset around test.elec
          
          for (i in 1:4){
              # make a sequence from the lowest possible inflection point to the test electrode 
              # of length equal to the number of straps calculated for "Straps below"
              # dist.diff.below <-  A[4 + i, test.elec] - CV.boundaries[i, 1]
              # dist.diff.above <-  A[4 + i, test.elec] - CV.boundaries[i, 2]
            
              # calculate d (length of sequence associated with electrodes above or below)
              # calcuylate distance difference between boundaries and test.elec
              if (straps.above > straps.below){
                d <- 10^(seq(sprd, 0, length = straps.above))
                dist.diff <- CV.boundaries[i, 2] - A[4 + i, test.elec] # test elec to max boundary
                #dist.diff.above <- CV.boundaries[i, 2] - A[4 + i, test.elec]
              }else{
                d <- 10^(seq(sprd, 0, length = straps.below))
                dist.diff <-  A[4 + i, test.elec] - CV.boundaries[i, 1] # test elec to min boundary
                #dist.diff <-  A[4 + i, test.elec] - CV.boundaries[i, 1]
              }
              
              d[which.max(d)] <- 0 # make the largest value (the most distal strap) equal to 0 (i.e., the target electrode's inflection points)
              
            
            ### ******* note dist.diff here is calculated in reference to the side (above/below) with more elecs******
            ### ******* We could make the dist.diff.below and dist.diff.above and use that 
            
            # of length equal to the number of straps calculated for "Straps below"
              below.seq <- 1:length(straps.below)
              seq.mat[i,1:straps.below] <-  A[4 + i, test.elec] - abs(dist.diff) * d[1:straps.below]
              # seq.mat[i,1:straps.below] <-  A[4 + i, test.elec] - dist.diff.below * d[length(straps.below)]

            
            # of length equal to the number of straps calculated for "Straps above"
              seq.mat[i,(straps.below + 1):straps] <- A[4 + i, test.elec] + abs(dist.diff) * d[1:(straps - straps.below)]
              #seq.mat[i,(straps.below + 1):straps] <- A[4 + i, test.elec] + dist.diff.above * d[straps.below - straps]
          }
          
          
          # determine whether more space above or below
          # determine proportions of straps above or below based on that
          # generate sequence. Make sure for the evenly spaced to give a little distance from test elec on either side so get close elctrodes on both sides
        }
    
  }else if(dist == 3){
          # just use distance metric
        dist.vec <- vector(length = ncol(A)) #vector of distance metrics
        for (z in 1:ncol(A)){
          dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
        }
        dist.diff <- max(dist.vec) # maximum distance metric from the test elec
        d <- 10^(seq(sprd, 0, length = straps)) # strap dist
        d[which.max(d)] <- 0
        seq.vec <-  abs(dist.diff) * d

        }

  # matrix of target inflection points of all target electrodes
  target.coord <- matrix(NA, ncol = 8, nrow = straps)
  colnames(target.coord) <- paste0("Inflect_", 1:8)
  row.names(target.coord) <- paste0("Strap_", 1:straps)
  target.coord[,5:8] <- t(seq.mat) 
  
  ###############################################################
  # Generate Predictions (Predicted voltage potentials)
  ###############################################################
  
  # make predictions (voltage potentials) on the corresponding model fit (for corresponding 
  # inflection point) on sequence of currents generated  -- round since voltage potentials are discrete
  
  
  for (i in 1:4){
    if (volts == 1){
      # target CVs voltage potentials based on a linear model 
      target.coord[,i] <- round( predict( get(paste0("lm.",i)), data.frame(x = seq.mat[i ,]) ) )
    }else if (volts == 2){
      # target CVs voltage potentials based on target elec's
      target.coord[,i] <- A[i, test.elec] # voltage potentials of the test elec
    }
  }
  
  if(dist == 3){
    # if we use distance measure only (not inflection points) then "coordinates" are just the distances from test elec
    # it would be just a vector
    target.coord <- seq.vec
  }
  
  return(target.coord)
  
}


######################################################
# Generate Target Distances
######################################################
targetDist <- function(full.CVs, test.elec, straps, sprd = - 2.5){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
    # just use distance metric
    dist.vec <- vector(length = ncol(A)) #vector of distance metrics
    for (z in 1:ncol(A)){
      dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
    }
    dist.diff <- max(dist.vec) # maximum distance metric from the test elec
    d <- 10^(seq(sprd, 0, length = straps)) # strap dist
    d[which.max(d)] <- 0
    seq.vec <-  abs(dist.diff) * d

  return(seq.vec)
  
}



######################################################
# Generate Target Distances
######################################################
targetDist.MinMax <- function(full.CVs, test.elec, straps, sprd = - 2.5){
    # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
    # test.elec is the number of the electrode that is the test electrode (integer \in [K])
    # straps is the number of study straps to generate
    # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
    # that the target electrodes should be
    
    source("Study Strap Functions.R") # needed for wiggle.fn
    #------------------------------------------
    # make Coordinate matrix A
    K <- nrow(full.CVs) 
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_",1:K)
    row.names(A) <- paste0("Coordinate_",1:8)
    
    for (i in 1:nrow(full.CVs)){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    # just use distance metric
    dist.vec <- vector(length = ncol(A)) #vector of distance metrics
    for (z in 1:ncol(A)){
        dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
    }
    dist.diff <- max(dist.vec) # maximum distance metric from the test elec
    d <- 10^(seq(sprd, 0, length = straps)) # strap dist
    d[which.max(d)] <- 0
    seq.vec <-  abs(dist.diff) * d
    
    return(seq.vec)
    
}


######################################################
# Generate Target Distances
######################################################
targetDist2 <- function(full.CVs, test.elec, straps, sprd = - 2.5, bag.size = nrow(full.CVs)){
    # starting values are the minimum distances (when d = 0) instead of 0
    # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
    # test.elec is the number of the electrode that is the test electrode (integer \in [K])
    # straps is the number of study straps to generate
    # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
    # that the target electrodes should be
    
    source("Study Strap Functions.R") # needed for wiggle.fn
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    #------------------------------------------
    # make Coordinate matrix A
    K <- nrow(full.CVs) 
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_",1:K)
    row.names(A) <- paste0("Coordinate_",1:8)
    
    for (i in 1:nrow(full.CVs)){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    # just use distance metric
    dist.vec <- vector(length = ncol(A)) #vector of distance metrics
    for (z in 1:ncol(A)){
        dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    A <- A[,-test.elec] # remove for optimization below

    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    # Constrained Optimization (Integer Least Squares)
    problem <- Problem(objective, constraints = list(w >= 0, 
                                                     sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    min.dist <- sum((elec.target - X %*% w.solution )^2)
    
    dist.diff <- abs(max(dist.vec) - min.dist) # difference between maximum distance metric from closest feasible study strap to  test elec
    d <- 10^(seq(sprd, 0, length = straps)) # strap dist
    d[which.max(d)] <- 0 # make the furthest distance a 0 since the study strap would be just the furthest electrode
    seq.vec <-  min.dist + abs(dist.diff) * d
    
    return(seq.vec)
    
}

######################################################
# Generate Target Distances
######################################################
targetDist.exp <- function(full.CVs, test.elec, straps, sprd = - 2.5){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  dist.diff <- max(dist.vec) # maximum distance metric from the test elec
  d <- exp(seq(sprd, 0, length = straps)) # strap dist
  d[which.max(d)] <- 0
  seq.vec <-  abs(dist.diff) * d
  
  return(seq.vec)
  
}





max.bag <- function(bag.mat, r){
  # provide a matrix where each row is a recipe for the bag and choose a subset of r bags that are as different as possible
  
  ###############
  # Preprocess
  ##############
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  final.bags <- unique(bag.mat[[2]]) # remove duplicate rows to save computation time
  bag.mat <- unique(bag.mat[[1]]) # remove duplicate rows to save computation time
  B <- bag.mat %*% t(bag.mat) # matrix of inner products between all studies
  diag(B) <- rep(0, nrow(B)) # remove all diagonal terms since these are just L2 norm of bag recipes. could keep but provides unclear regularization on bags that are less varied
  z <- B %*% rep(1, nrow(B)) # vector of sum of innter products between all studies
  
  bags <- 1:nrow(bag.mat) # indices of bags to include
  
  # initialize matrices
  B.work <- B # initialize working B matrix
  bag.work <- bag.mat # initialize working bag matrix
  
  while( nrow(B.work) > r){
    # until we have selected r studies, keep iterating
    #-----------
    s <- Int( nrow(B.work) ) # vector to optimize of length of number of bags to consider. Constrain to be integer
    
    objective <- Minimize(t(z) %*% s)  # 1/K part is to 
    problem <- Problem(objective, constraints = list(s >= 0, s <= 1, sum(s) == (nrow(B.work) - 1) ) ) # sum so that only get rid of one study at a time
    result <- solve(problem)
    s.solution <- round(result$getValue(s)) # round since the 0s are not exactly 0
    
    # reinitialize matrices
    elim <- which(s.solution == 0) # check which study to eliminate
    bags <- bags[-elim] # eliminate the bag associated with the row that was deleted
    bag.work <- bag.work[-elim, ] # remove from bag matrix the row/bag that was eliminated
    B.work <- bag.work %*% t(bag.work) # make inner product matrix
    diag(B.work) <- rep(0, nrow(B.work)) # remove all diagonal terms since these are just L2 norm of bag recipes. could keep but provides unclear regularization on bags that are less varied
    z <- B.work %*% rep(1, nrow(B.work)) # vector of sum of inner products between all studies
    
  }
  
  return(list(final.bags[bags,], bag.mat[bags,]))
}


######################################################
# Generate Target Epsilons to bound 
######################################################
targetEps <- function(full.CVs, test.elec, eta = 100, bag.size = nrow(full.CVs) - 1){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  # eta is step size
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  elec.target <- A[,test.elec] # coordinates of test electrode
  A <- A[,-test.elec] # delete test.elec column
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  
  max.dist <- max(dist.vec) # maximum distance metric from the test elec
  
  ### find minimum distance
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  # Constrained Optimization (Integer Least Squares)
  problem <- Problem(objective, constraints = list(w >= 0, 
                                                     sum(w[1:bag.size]) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  
  min.dist <- sum((elec.target - X %*% w.solution )^2)
  seq.vec <-  seq(min.dist, max.dist, by = eta)
  
  return(seq.vec)
  
}



####################################
# Bag Generator that minimizes ||w||^2
####################################

bagGen.eps <- function(full.CVs, epsilon.vec, test.elec, bag.size = nrow(full.CVs) - 1, 
                                straps = nrow(full.CVs)){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  
  epsilon.vec <- as.vector(epsilon.vec)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  flag <- TRUE
  i <- 1 # accepted bag counter
  z <- 1 # epsilon index counter
  # target electrode 
  
  while(flag == TRUE){
    # integer program
    eps <- epsilon.vec[y]
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    objective <- Minimize(norm2(w))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == K - 1,
                                                     sum((elec.target - X %*% w )^2) <= eps))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0

    # ensure do not proceed if no solution found (only length 1 if solution is NA)
    if (length(w.solution) > 1){ 
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      
      if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x), 
                                                              as.numeric(want))), w.solution))){
          # if the bag has not yet been added to the bag matrix (i.e., it is unique) 
          # then add to the matrix
          
          # make bag
          bag <- c()
          for (x in 1:length(w.solution)){
              bag <- c(bag, rep(x, w.solution[x]))
          }
          
          # add to matrices
          bag.mat[i,] <- study.num[bag] # add bag to matrix.
          bag.recipe[i, ] <- w.solution  # add bag recipe to 
          i <- i + 1 # update bag counter
      }
    }
    
    z <- z + 1 # update epsilon index counter
    
    if (y > length(epsilon.vec) | i > straps){
      # break while loop if no epsilons left or we have achieved desired strap number
        flag <- FALSE
    }
    
    
  }
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


###############
# Bag Selector
##############
# selects bags that maximize between-bag variability
bagSelect <- function(bag.mat, r){
  # provide a matrix where each row is a recipe for the bag and choose a 
  # subset of r bags that are as different as possible
  
  bag.recipes <- unique(bag.mat[[2]]) # remove duplicate rows to save computation time
  bag.mat <- unique(bag.mat[[1]]) # remove duplicate rows of bag recpies--each row is a recipe of a seperate bag
  
  ###############
  # Preprocess
  ##############
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  B <- matrix(NA, ncol = nrow(bag.mat), nrow = nrow(bag.mat)) # matrix of l2 norms of differences
  
  # matrix of differences of l2 norms of bag recipes
  for ( i in 1:nrow(bag.recipes)){
      for (j in 1:nrow(bag.recipes)){
          B[i,j] <- sum( (bag.recipes[i, ] - bag.recipes[j,])^2 )
      }
  }
  bags <- 1:nrow(bag.mat) # indices of bags to include
  
  # initialize matrices
  # B.work <- B # initialize working B matrix
  
  while( nrow(B) > r){
      # quad_form(w, B)  --- could use this
    # until we have selected r studies, keep iterating
    #elim <- which.max(rowSums(B)) # round since the 0s are not exactly 0
    elim <- which.min(rowSums(B))
    # reinitialize matrices
    bags <- bags[-elim] # eliminate the bag associated with the row that was deleted
    B <- B[-elim, -elim]
  
  }
 

  #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
  # so need to use list of elecs (study num) and index them by bag
  
  
    # B <- as.matrix(B)
    # s <- Int( nrow(B) ) # vector of bags to select. Constrain each element to be an integer \in {0,1}
    # B.s <- B %*% s
    # 
    # B.mat <- quad_form(s, B.s)
    # 
    # objective <- Minimize(-1 * B.s %*% s ) # 1/K part is to 
    # #objective <- Minimize(norm2( B %*% s) ) # 1/K part is to 
    # problem <- Problem(objective, constraints = list(s >= 0, s <= 1, sum(s) == r) ) # sum so that end up with r studies
    # result <- solve(problem)
    # s.solution <- round(result$getValue(s)) # round since the 0s are not exactly 0
    
    
  
    # return rows of selected bags
  return(list(bag.mat[bags,], bag.recipes[bags,]))
}


mat_diff <- function(mat){
  mat <- mat[rowSums(!is.na(mat)) > 0,] # remove NAs
  B <- matrix(NA, nrow = nrow(mat), ncol = nrow(mat))
  for ( i in 1:nrow(mat)){
    for (j in 1:nrow(mat)){
      B[i,j] <- sum( (mat[i, ] - mat[j,])^2 )
    }
  }
  return( sum(B) / nrow(B)   )  # average similarity
}



##############
# Stacking
##############
# this assumes everything is in global environment
stackCoefs <- function(x){
  # make stacking matrix
  for (ss_mat in 1:row.total){ # iterate through Study Strap design matrices
    
    #---------------------------------------------------
    for (SSL in 1:row.total){ # iterate through Single Study Learners (models)
      # stack.mat is main matrix
      # current.ss is current study strap matrix which is appended to main stack.mat
      if ( SSL == 1 ){
        # first learner of the current study strap -- bind labels (Y - outcome) with predictions made with first SSL
        current.ss <- cbind( get(paste0("study_", y, "_ss_", ss_mat))$DA, 
                             predict( get(paste0("study_", y, "_learner_", SSL)), 
                                      get(paste0("study_", y, "_ss_", ss_mat))[,-1]) )
      }
      else{
        # subsequent learners of the current study strap -- bind predictions made with current SSL
        current.ss <- cbind( current.ss, 
                             predict( get(paste0("study_", y, "_learner_", SSL)), 
                                      get(paste0("study_", y, "_ss_", ss_mat))[,-1]) )            
      }
    }
    
    rm(list = paste0("study_", y, "_ss_", ss_mat)) # Delete current Study Strap to save memory
    
    
    ####
    # Create Stack Matrix based on predictions made on a single Study Strap Design Matrix based on models from all SSLs (s)
    ####
    if (ss_mat == 1){
      # if its the first study strap, then make the stack.mat matrix based on the first matrix of predictions
      stack.mat <- current.ss
      rm(current.ss) # save memory
    }
    else{
      # if its not the first study strap, then rbind the stack.mat matrix with current iteration
      stack.mat <- rbind(stack.mat, current.ss)
      rm(current.ss) # save memory
    }
    ###################
    #---------------------------------------------------
    
    # delete current study strap to save memory:
    rm(list = paste0("study_",y,"_ss_", ss_mat))
    
  }
  # end stacking matrix prep.
  
  # delete models to save memory:
  rm(list = paste0("study_", y, "_learner_", 1:straps))
  # delete study straps to save memory:
  #rm(list = paste0("study_",y,"_ss_", 1:straps))
  
  ###################
  # Stack Regression Weights 
  ###################
  
  # nnls
  print("Stack nnls")
  library(nnls)
  nnls.coefs <- coef(nnls(A = as.matrix(stack.mat[,-1]), b = as.matrix(stack.mat[,1])))
  nnls.coefs.int <- coef(nnls(A = as.matrix(cbind(1,stack.mat[,-1])), b = as.matrix(stack.mat[,1])))
  
  print(paste0("dim RF Preds: ", dim(RF.preds)))
  print(paste0("dim stack.mat: ", dim(stack.mat)))
  print(paste0("length of stack nnls: ", length(nnls.coefs)))
  
  library(glmnet)
  print("Stack ridge")
  
  # ridge
  lambdas <- 10^seq(3, -2, by = -.2) # -.1 gives 50 lamba values
  # #cv_fit <- cv.glmnet(x, y, alpha = 0, lambda = lambdas)
  cv_fit <- cv.glmnet(x = (as.matrix(stack.mat[,-1])), y = as.matrix(stack.mat[,1]),
                      family = "gaussian", alpha = 0, lambda = lambdas, 
                      x = FALSE, model = FALSE, y = FALSE)
  
  ridge.coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
  #ridge.coefs <- nnls.coefs.int # save time
  rm(stack.mat, cv_fit) # save memory
  print(paste0("length of stack ridge: ", length(ridge.coefs)))
  return(list(nnls.coefs, nnls.coefs.int, ridge.coefs))
}



##############
# Stacking
##############
# this assumes everything is in global environment
stack.fn <- function(x){
  # make stacking matrix
  for (ss_mat in 1:row.total){ # iterate through Study Strap design matrices
    
    #---------------------------------------------------
    for (SSL in 1:row.total){ # iterate through Single Study Learners (models)
      # stack.mat is main matrix
      # current.ss is current study strap matrix which is appended to main stack.mat
      if ( SSL == 1 ){
        # first learner of the current study strap -- bind labels (Y - outcome) with predictions made with first SSL
        current.ss <- cbind( full[get(paste0("study_", y, "_ss_", ss_mat)),1], 
                             predict( get(paste0("study_", y, "_learner_", SSL)), 
                                      full[get(paste0("study_", y, "_ss_", ss_mat)),-1]) )
      }
      else{
        # subsequent learners of the current study strap -- bind predictions made with current SSL
        current.ss <- cbind( current.ss, 
                             predict( get(paste0("study_", y, "_learner_", SSL)), 
                                      full[get(paste0("study_", y, "_ss_", ss_mat)),-1]) )            
      }
    }
    
    rm(list = paste0("study_", y, "_ss_", ss_mat)) # Delete current Study Strap to save memory
    
    
    ####
    # Create Stack Matrix based on predictions made on a single Study Strap Design Matrix based on models from all SSLs (s)
    ####
    if (ss_mat == 1){
      # if its the first study strap, then make the stack.mat matrix based on the first matrix of predictions
      stack.mat <- current.ss
      rm(current.ss) # save memory
    }
    else{
      # if its not the first study strap, then rbind the stack.mat matrix with current iteration
      stack.mat <- rbind(stack.mat, current.ss)
      rm(current.ss) # save memory
    }
    ###################
    #---------------------------------------------------
    
    # delete current study strap to save memory:
    rm(list = paste0("study_",y,"_ss_", ss_mat))
    
  }
  # end stacking matrix prep.
  
  # delete models to save memory:
  rm(list = paste0("study_", y, "_learner_", 1:straps))
  # delete study straps to save memory:
  #rm(list = paste0("study_",y,"_ss_", 1:straps))
  
  ###################
  # Stack Regression Weights 
  ###################
  
  # nnls
  print("Stack nnls")
  library(nnls)
  nnls.coefs <- coef(nnls(A = as.matrix(stack.mat[,-1]), b = as.matrix(stack.mat[,1])))
  nnls.coefs.int <- coef(nnls(A = as.matrix(cbind(1,stack.mat[,-1])), b = as.matrix(stack.mat[,1])))
  
  print(paste0("dim RF Preds: ", dim(RF.preds)))
  print(paste0("dim stack.mat: ", dim(stack.mat)))
  print(paste0("length of stack nnls: ", length(nnls.coefs)))
  
  library(glmnet)
  print("Stack ridge")
  
  # ridge
  lambdas <- 10^seq(3, -2, by = -.2) # -.1 gives 50 lamba values
  # #cv_fit <- cv.glmnet(x, y, alpha = 0, lambda = lambdas)
  cv_fit <- cv.glmnet(x = (as.matrix(stack.mat[,-1])), y = as.matrix(stack.mat[,1]),
                      family = "gaussian", alpha = 0, lambda = lambdas, 
                      x = FALSE, model = FALSE, y = FALSE)
  
  ridge.coefs <- as.vector(coef(cv_fit, s = "lambda.min"))
  #ridge.coefs <- nnls.coefs.int # save time
  rm(stack.mat, cv_fit) # save memory
  print(paste0("length of stack ridge: ", length(ridge.coefs)))
  return(list(nnls.coefs, nnls.coefs.int, ridge.coefs))
}


# calculate minimum distance based on d = 0
# add a constant \eta to be your epsilon
# 

# bagGen.CF(full.CVs = full.CVs, eta = 100, test.elec = y, bag.size = bag.size, 
#           straps = straps)

bagGen.CF <- function(full.CVs, eta = 100, test.elec, bag.size = nrow(full.CVs) - 1, 
                       straps = nrow(full.CVs)){
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  dists <- vector(length= (K - 1))
  for (i in seq(1,K)){
    dists[i] <- sum( (A[,i] - A[,test.elec])^2)
  }
  dist.se <- sd(dists[-test.elec]) #sqrt(var(dists)) # standard error of distances from test elec
  print(dist.se)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  flag <- TRUE
  i <- 1 # accepted bag counter
  z <- 1 # epsilon index counter
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  min.dist <- sum((elec.target - X %*% w.solution )^2)
  eps <- 1 + min.dist # radius of the campfire
  
  #### add this to be the first bag
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }
  
  # add to matrices
  bag.mat[1,] <- study.num[bag] # add bag to matrix.
  bag.recipe[1, ] <- w.solution  # add bag recipe to 
  
  
  # target electrode 
  
  while(flag == TRUE){
    # integer program
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    B <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
    objective <- Minimize( sum( B %*% w ))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size,
                                                     sum((elec.target - X %*% w )^2) <= eps))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # ensure do not proceed if no solution found (only length 1 if solution is NA)
    if (length(w.solution) > 1){ 
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      
      if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x), 
                                                                       as.numeric(want))), w.solution))){
        # if the bag has not yet been added to the bag matrix (i.e., it is unique) 
        # then add to the matrix
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        
        # add to matrices
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
        i <- i + 1 # update bag counter
        print(i)
      }
    }else{
      # if NA then expand the campfire radius
      eps <- eps + min.dist # expand radius by factor of how big the minimum distance is
    }
    
    # z <- z + 1 # update epsilon index counter
    
    if (i > straps){
      # break while loop if no epsilons left or we have achieved desired strap number
      flag <- FALSE
    }
    
    
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.dist2 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0, limit = 20){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  
  distances <- as.vector(c(0,distances[1:limit]))
  print(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:length(distances))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  
  for (i in 1:length(distances)){
    # integer program
    d <- distances[i]
    print(d)
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    objective <- Minimize( pos( norm2(elec.target - X %*% w) - d) + lambda * norm2(w) )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))     
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
      bag <- c(bag, rep(x, w.solution[x]))
    }
    #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
    # so need to use list of elecs (study num) and index them by bag
    bag.mat[i,] <- study.num[bag] # add bag to matrix.
    bag.recipe[i, ] <- w.solution  # add bag recipe to 
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}



bagGen.L2Con <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, 
                      straps = nrow(full.CVs)){
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  dists <- vector(length= (K - 1))
  for (i in seq(1,K)){
    dists[i] <- sum( (A[,i] - A[,test.elec])^2)
  }
  dist.se <- sd(dists[-test.elec]) #sqrt(var(dists)) # standard error of distances from test elec
  print(dist.se)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  flag <- TRUE
  i <- 1 # accepted bag counter
  z <- 1 # epsilon index counter
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  # min.dist <- sum((elec.target - X %*% w.solution )^2)
  # eps <- 1 + min.dist # radius of the campfire
  w.min <- sum(w.solution^2) - 1 # subtract one so we do not repeat it
  w.vec <- seq(w.min, K-1)
  
  #### add this to be the first bag
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }
  
  # add to matrices
  bag.mat[1,] <- study.num[bag] # add bag to matrix.
  bag.recipe[1, ] <- w.solution  # add bag recipe to 
  
  
  # target electrode 
  
  while(flag == TRUE){
    # integer program
    w.limit <- w.vec[z]
    objective <- Minimize(sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size, sum(w^2) <= w.limit)) 
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # ensure do not proceed if no solution found (only length 1 if solution is NA)
    if (length(w.solution) > 1){ 
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      
      if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x), 
                                                                       as.numeric(want))), w.solution))){
        # if the bag has not yet been added to the bag matrix (i.e., it is unique) 
        # then add to the matrix
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        
        # add to matrices
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
        i <- i + 1 # update bag counter
        print(i)
        }else{
          # if NA then move to next w limit
          z <- z + 1 
        }
      
    }else{
      # if repeat then move to next w limit
      z <- z + 1 # expand radius by factor of how big the minimum distance is
    }
    
    
    if (i > straps | z > length(w.vec) ){
      # break while loop if no w limits left or we have achieved desired strap number
      flag <- FALSE
    }
    
    
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.L2Constrained <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, 
                         straps = nrow(full.CVs)){
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  dists <- vector(length= (K - 1))
  for (i in seq(1,K)){
    dists[i] <- sum( (A[,i] - A[,test.elec])^2)
  }
  dist.se <- sd(dists[-test.elec]) #sqrt(var(dists)) # standard error of distances from test elec
  print(dist.se)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  flag <- TRUE
  i <- 1 # accepted bag counter
  z <- 1 # epsilon index counter
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  # min.dist <- sum((elec.target - X %*% w.solution )^2)
  # eps <- 1 + min.dist # radius of the campfire
  w.min <- sum((w.solution)^2)
  w.limit <- w.min # initiate at this value
  
  #### add this to be the first bag
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }
  
  # add to matrices
  bag.mat[1,] <- study.num[bag] # add bag to matrix.
  bag.recipe[1, ] <- w.solution  # add bag recipe to 
  
  
  # target electrode 
  
  while(flag == TRUE){
    # integer program

    objective <- Minimize(sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size, sum(w^2) < w.limit)) 
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # ensure do not proceed if no solution found (only length 1 if solution is NA)
    if (length(w.solution) > 1){ 
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      
      if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x), 
                                                                       as.numeric(want))), w.solution))){
        # if the bag has not yet been added to the bag matrix (i.e., it is unique) 
        # then add to the matrix
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        
        # add to matrices
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
        i <- i + 1 # update bag counter
        w.limit <- sum((w.solution)^2) #update w threshold
      }else{
        # if NA then break loop
        z <- z + 1 
        w.limit <- w.limit - 1
      }
      
    }else{
      # if repeat then move to next w limit
      z <- z + 1 # expand radius by factor of how big the minimum distance is
    }
    
    
    if (i > straps | w.limit < bag.size){
      # break while loop if no w limits left or we have achieved desired strap number
      flag <- FALSE
      print("break ||w||^2 threshold below bag.size")
      print(w.limit)
      
    }
    
    
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}




####################################
# Bag Generator Function Based on distances
####################################

bagGen.dist10 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0, voltWeights = 10){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # multiple voltage potentials by voltWeights to weight their importance
  
  distances <- as.vector(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:length(distances))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  
  for (i in 1:length(distances)){
    # integer program
    d <- distances[i]
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    objective <- Minimize( pos( norm2(elec.target - X %*% w) - d) + lambda * norm2(w) )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))     
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    if( length(w.solution) > 1){
      # if solution is not NA
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}



######################################################
# Generate Target Distances
######################################################
targetDist10 <- function(full.CVs, test.elec, straps, sprd = - 2.5, voltWeights = 10){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  dist.diff <- max(dist.vec) # maximum distance metric from the test elec
  d <- 10^(seq(sprd, 0, length = straps)) # strap dist
  d[which.max(d)] <- 0
  seq.vec <-  abs(dist.diff) * d
  
  return(seq.vec)
  
}


######################################################
# Generate Target Distances
######################################################
targetDistCV <- function(full.CVs, test.elec, straps, sprd = - 2.5, voltWeights = 10){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  # A <- matrix(NA, nrow = 8, ncol = K)
  # colnames(A) <- paste0("Elec_",1:K)
  # row.names(A) <- paste0("Coordinate_",1:8)
  
  
  #A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  A <- t(full.CVs)
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  dist.diff <- max(dist.vec) # maximum distance metric from the test elec
  d <- 10^(seq(sprd, 0, length = straps)) # strap dist
  d[which.max(d)] <- 0
  seq.vec <-  abs(dist.diff) * d
  
  return(seq.vec)
  
}


####################################
# Bag Generator Function Based on distances
####################################

bagGen.distCV <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0, voltWeights = 10){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # multiple voltage potentials by voltWeights to weight their importance
  
  distances <- as.vector(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  # A <- matrix(NA, nrow = 8, ncol = K)
  # colnames(A) <- paste0("Elec_", 1:K)
  # row.names(A) <- paste0("Coordinate_", 1:8)
  # for (i in 1:K){
  #   A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  # }
  # 
  # A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  A <- t(full.CVs)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:length(distances))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  
  for (i in 1:length(distances)){
    # integer program
    d <- distances[i]
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    objective <- Minimize( pos( norm2(elec.target - X %*% w) - d) + lambda * norm2(w) )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))     
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    if( length(w.solution) > 1){
      # if solution is not NA
      # make bag
      bag <- c()
      for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
      }
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      bag.mat[i,] <- study.num[bag] # add bag to matrix.
      bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}






# diff

######################################################
# Generate Target Distances
######################################################
targetDistDeriv <- function(full.CVs, test.elec, straps, sprd = - 2.5, voltWeights = 1){
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  # A <- matrix(NA, nrow = 8, ncol = K)
  # colnames(A) <- paste0("Elec_",1:K)
  # row.names(A) <- paste0("Coordinate_",1:8)
  
  #A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  full.CVs <- t(diff(t(full.CVs))) # derivative
  A <- t(full.CVs)
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  dist.diff <- max(dist.vec) # maximum distance metric from the test elec
  d <- 10^(seq(sprd, 0, length = straps)) # strap dist
  d[which.max(d)] <- 0
  seq.vec <-  abs(dist.diff) * d
  
  return(seq.vec)
  
}


####################################
# Bag Generator Function Based on distances
####################################

bagGen.distDeriv <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0, voltWeights = 1){
  # input is coordinates if inflection point coordinates are fed
  # input is distance if distances are fed
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # multiple voltage potentials by voltWeights to weight their importance
  
  distances <- as.vector(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  # A <- matrix(NA, nrow = 8, ncol = K)
  # colnames(A) <- paste0("Elec_", 1:K)
  # row.names(A) <- paste0("Coordinate_", 1:8)
  # for (i in 1:K){
  #   A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  # }
  # 
  # A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
  full.CVs <- t(diff(t(full.CVs))) # derivative
  A <- t(full.CVs)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:length(distances))
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # target electrode 
  
  for (i in 1:length(distances)){
    # integer program
    d <- distances[i]
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    objective <- Minimize( pos( norm2(elec.target - X %*% w) - d)  ) # + lambda * norm2(w)
    #objective <- Minimize( sum((elec.target - X %*% w )^2))     
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    if( length(w.solution) > 1){
      # if solution is not NA
      # make bag
      bag <- c()
      for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
      }
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      bag.mat[i,] <- study.num[bag] # add bag to matrix.
      bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


######################################################
# Generate Target Distances
######################################################
targetDistMin <- function(full.CVs, test.elec, straps, sprd = - 4, bag.size = nrow(full.CVs) - 1){
  # like target dist but starting values is from closest distance to test electrode (instead of 0)
  # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
  # test.elec is the number of the electrode that is the test electrode (integer \in [K])
  # straps is the number of study straps to generate
  # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
  # that the target electrodes should be
  
  source("Study Strap Functions.R") # needed for wiggle.fn
  #------------------------------------------
  # make Coordinate matrix A
  K <- nrow(full.CVs) 
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_",1:K)
  row.names(A) <- paste0("Coordinate_",1:8)
  
  for (i in 1:nrow(full.CVs)){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  elec.target <- A[,test.elec] # the "target" in the objective
  A <- A[,-test.elec]
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  min.dist <- sum((elec.target - X %*% w.solution )^2) # closest distnaces to test electrode
  
  # just use distance metric
  dist.vec <- vector(length = ncol(A)) #vector of distance metrics
  for (z in 1:ncol(A)){
    dist.vec[z] <- sum( (A[,z] - A[,test.elec])^2 )
  }
  dist.diff <- max(dist.vec) # maximum distance metric from the test elec
  d <- 10^(seq(sprd, min.dist, length = straps)) # strap dist
  d[which.max(d)] <- 0
  seq.vec <-  abs(dist.diff) * d
  
  return(seq.vec)
  
}




bagGen.CFDist <- function(full.CVs, eta = 20, test.elec, distances, bag.size = nrow(full.CVs) - 1, 
                      straps = nrow(full.CVs), step.size = 10){
  # like bagGen.CF but takes in distances for constraints
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  # eta is the radius
  # 
  
  distances <- as.vector(distances)
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  dists <- vector(length= (K - 1))
  for (i in seq(1,K)){
    dists[i] <- sum( (A[,i] - A[,test.elec])^2)
  }
  dist.se <- sd(dists[-test.elec]) #sqrt(var(dists)) # standard error of distances from test elec
  print(dist.se)
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = K - 1, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  flag <- TRUE
  i <- 1 # accepted bag counter
  z <- 1 # epsilon index counter
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  min.dist <- (sum((elec.target - X %*% w.solution )^2)) # sqrt because constraint uses l2 norm not the squared norm
  eps <- min.dist + eta # radius of the campfire
  
  #### add this to be the first bag
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }
  
  # add to matrices
  bag.mat[1,] <- study.num[bag] # add bag to matrix.
  bag.recipe[1, ] <- w.solution  # add bag recipe to 
  
  
  # target electrode 
  
  while(flag == TRUE){
    d <- distances[i]
    # integer program
    # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
    B <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
    objective <- Minimize( sum( B %*% w ))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size,
                                                     pos( sum((elec.target - X %*% w)^2) - d) <= eps))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # ensure do not proceed if no solution found (only length 1 if solution is NA)
    if (length(w.solution) > 1){ 
      #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
      # so need to use list of elecs (study num) and index them by bag
      
      if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x), 
                                                                       as.numeric(want))), w.solution))){
        # if the bag has not yet been added to the bag matrix (i.e., it is unique) 
        # then add to the matrix
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
          bag <- c(bag, rep(x, w.solution[x]))
        }
        
        # add to matrices
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
        i <- i + 1 # update bag counter
        eps <- min.dist + eta # restart radius
        print(i)
      }else{
        eps <- eps + step.size # expand radius by factor of step.size
      }
    }else{
      # if NA then expand the campfire radius
      eps <- eps + step.size # expand radius by factor of step.size
    }
    
    # z <- z + 1 # update epsilon index counter
    
    if (i > straps){
      # break while loop if no epsilons left or we have achieved desired strap number
      flag <- FALSE
    }
    
    
  }
  bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
  bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


#########################
# Only Closest Elec
#########################

bagGen0 <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, straps = 1){
  # like bagGen.CF but takes in distances for constraints
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  # eta is the radius
  # 
 
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
 
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  
  #### add this to be the first bag
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }
  
  # add to matrices
  bag.mat <- study.num[bag] # add bag to matrix.
  bag.recipe <- w.solution  # add bag recipe to 
  
  return( list( (bag.mat), (bag.recipe)  ) )  # not unique because it is a single vector
}


#########################
# Only Closest Elec Penalize dot products (between bag variability)
#########################

bagGen_B <- function(full.CVs, test.elec, distances, bag.size = nrow(full.CVs) - 1, straps = 5, lambda = 0.001){
  # like bagGen.CF but takes in distances for constraints
  # input is coordinates if inflection point coordinates are fed
  # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
  # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
  # and average CV matrix without test elec and outputs a matrix of study bags
  # straps are number to go until stopping
  # eta is the radius
  # 
  
  # import wiggle.fn
  source("Study Strap Functions.R")
  
  
  K <- nrow(full.CVs)# number of studies including test elec
  study.num <- seq(1, K)[-test.elec] # the elecs included
  # make objective matrix A
  A <- matrix(NA, nrow = 8, ncol = K)
  colnames(A) <- paste0("Elec_", 1:K)
  row.names(A) <- paste0("Coordinate_", 1:8)
  for (i in 1:K){
    A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
  }
  
  elec.target <- A[,test.elec] # the "target" in the objective
  
  A <- A[,-test.elec] # remove for optimization below
  
  # bag matrix
  bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
  bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
  colnames(bag.recipe) <- paste0("P_",study.num)
  row.names(bag.recipe) <- paste0("Target",1:straps)
  
  #### Prepare Constrained Optimization ####
  # load package
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  
  w <- Int(K - 1)
  X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
  
  # minimum dist
  objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
  problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
  result <- solve(problem)
  w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
  B <- matrix(NA, ncol = bag.size, nrow = 1) # matrix of l2 norms of differences
  B[1,] <- w.solution
  
  # make bag
  bag <- c()
  for (x in 1:length(w.solution)){
    bag <- c(bag, rep(x, w.solution[x]))
  }

  bag.mat[1,] <- study.num[bag] # add bag to matrix.
  bag.recipe[1,] <- w.solution  
  
  itr <- 2
  for(itr in 2:straps){
    # minimum dist
    lambda.norm <- lambda / (itr - 1)  # update this so we don't penalize more with increasing numbers of additions 
    objective <- Minimize(sum((elec.target - X %*% w )^2) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
    B <- rbind(B, w.solution) # update B
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
      bag <- c(bag, rep(x, w.solution[x]))
    }
    # add to matrices
    bag.mat[itr,] <- study.num[bag] # add bag to matrix.
    bag.recipe[itr,] <- w.solution  # add bag recipe to 
  }
  
  #### add this to be the first bag

  return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


#########################
# Only Closest Elec Penalize dot products (between bag variability)
#########################

bagGen_B2 <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, straps = 5, lambda = 0.001){
    # like bagGen_B but includes min distance from test elec so not overly penalized for electrodes that are far
    # input is coordinates if inflection point coordinates are fed
    # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    # straps are number to go until stopping
    # eta is the radius
    # 
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    print("A matrix")
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:straps)
    
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    print("opt step")
    # minimum dist
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    B <- matrix(NA, ncol = K - 1, nrow = 1) # matrix of l2 norms of differences
    B[1,] <- w.solution
    min.dist <- sum((elec.target - X %*% w.solution )^2)
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    
    bag.mat[1,] <- study.num[bag] # add bag to matrix.
    bag.recipe[1,] <- w.solution  
    
    print("iterate opt steps")
    for(itr in 2:straps){
        # minimum dist
        lambda.norm <- lambda / (itr - 1)  # update this so we don't penalize more with increasing numbers of additions 
        objective <- Minimize(pos(sum((elec.target - X %*% w )^2) - min.dist) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
        result <- solve(problem)
        w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
        B <- rbind(B, w.solution) # update B
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        # add to matrices
        print("add opt solutions")
        bag.mat[itr,] <- study.num[bag] # add bag to matrix.
        bag.recipe[itr,] <- w.solution  # add bag recipe to 
    }
    
    #### add this to be the first bag
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}




bagGen_B3 <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, 
                     straps = 20, lambda = 1, lim = 50, phi = 1.25, voltWeights = 1){
    # like bagGen.CF but takes in distances for constraints
    # input is coordinates if inflection point coordinates are fed
    # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    # straps are number to go until stopping
    # phi is the update rate
    # lim is the limit of tries
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    d.vec <- c()
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:straps)
    
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # minimum dist
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    d.vec <- c(d.vec,   sum((elec.target - X %*% w.solution )^2)  )
    
    B <- matrix(NA, ncol = bag.size, nrow = 1) # matrix of l2 norms of differences
    B[1,] <- w.solution
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    
    bag.mat[1,] <- study.num[bag] # add bag to matrix.
    bag.recipe[1,] <- w.solution  
    itr <- 2
    flag <- TRUE
    z <- 0
    
    while(flag == TRUE){
        # minimum dist
        lambda.norm <- lambda / (itr - 1)  # update this so we don't penalize more with increasing numbers of additions 
        objective <- Minimize(sum((elec.target - X %*% w )^2) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
        result <- solve(problem)
        w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
        
        if (length(w.solution) > 1){
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            
            if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x),
                                                                             as.numeric(want))), w.solution))){
                # if the bag has not yet been added to the bag matrix (i.e., it is unique)
                # then add to the matrix
                
                B <- rbind(B, w.solution) # update B
                d.vec <- c(d.vec,  sum((elec.target - X %*% w.solution )^2)  )
                # make bag
                bag <- c()
                for (x in 1:length(w.solution)){
                    bag <- c(bag, rep(x, w.solution[x]))
                }
                
                # add to matrices
                bag.mat[itr,] <- study.num[bag] # add bag to matrix.
                bag.recipe[itr, ] <- w.solution  # add bag recipe to
                itr <- itr + 1 # update bag counter
            }else{
                # if NA then move to next w limit
                lambda <- lambda * phi
                z <- z + 1
            }
            
        }else{
            # if repeat bag then move to next w limit
            lambda <- lambda * phi # expand radius by factor of how big the minimum distance is
            z <- z + 1
        }
        
        
        if (itr > straps | z > lim){
            # break while loop if limit reached or we have achieved desired strap number
            flag <- FALSE
        }
    }
    print(paste0("lambda for test elec_", test.elec,"_", lambda))
    print(paste0("d.vec for test elec_", test.elec))
    print(d.vec)
    
    #### add this to be the first bag
    
    return( list( unique(bag.mat), unique(bag.recipe), d.vec  ) ) 
}

# 
# 
# bagGen_B4 <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, 
#                       straps = 20, lambda = 1, lim = 50, phi = 1.25, voltWeights = 1){
#     # like bagGen_B3 but includes a min dist in the objective so penalization is all with respect to closest elec
#     # input is coordinates if inflection point coordinates are fed
#     # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
#     # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
#     # and average CV matrix without test elec and outputs a matrix of study bags
#     # straps are number to go until stopping
#     # phi is the update rate
#     # lim is the limit of tries
#     
#     # import wiggle.fn
#     source("Study Strap Functions.R")
#     
#     d.vec <- c()
#     
#     K <- nrow(full.CVs)# number of studies including test elec
#     study.num <- seq(1, K)[-test.elec] # the elecs included
#     # make objective matrix A
#     A <- matrix(NA, nrow = 8, ncol = K)
#     colnames(A) <- paste0("Elec_", 1:K)
#     row.names(A) <- paste0("Coordinate_", 1:8)
#     for (i in 1:K){
#         A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
#     }
#     A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
#     
#     elec.target <- A[,test.elec] # the "target" in the objective
#     
#     A <- A[,-test.elec] # remove for optimization below
#     
#     # bag matrix
#     bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
#     bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
#     colnames(bag.recipe) <- paste0("P_",study.num)
#     row.names(bag.recipe) <- paste0("Target",1:straps)
#     
#     #### Prepare Constrained Optimization ####
#     # load package
#     suppressWarnings(library(CVXR, warn.conflicts=FALSE))
#     
#     w <- Int(K - 1)
#     X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
#     
#     # minimum dist
#     objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
#     problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
#     result <- solve(problem)
#     w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
#     min.dist <- sum((elec.target - X %*% w.solution )^2)
#     d.vec <- c(d.vec, min.dist )
#     
#     B <- matrix(NA, ncol = bag.size, nrow = 1) # matrix of l2 norms of differences
#     B[1,] <- w.solution
#     
#     # make bag
#     bag <- c()
#     for (x in 1:length(w.solution)){
#         bag <- c(bag, rep(x, w.solution[x]))
#     }
#     
#     bag.mat[1,] <- study.num[bag] # add bag to matrix.
#     bag.recipe[1,] <- w.solution  
#     itr <- 2
#     flag <- TRUE
#     z <- 0
#     
#     while(flag == TRUE){
#         # minimum dist
#         lambda.norm <- lambda / (itr - 1)  # update this so we don't penalize more with increasing numbers of additions 
#         objective <- Minimize(pos(sum((elec.target - X %*% w )^2) - min.dist) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
#         problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
#         result <- solve(problem)
#         w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
#         
#         if (length(w.solution) > 1){
#             #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
#             # so need to use list of elecs (study num) and index them by bag
#             
#             if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x),
#                                                                              as.numeric(want))), w.solution))){
#                 # if the bag has not yet been added to the bag matrix (i.e., it is unique)
#                 # then add to the matrix
#                 
#                 B <- rbind(B, w.solution) # update B
#                 d.vec <- c(d.vec,  sum((elec.target - X %*% w.solution )^2)  )
#                 # make bag
#                 bag <- c()
#                 for (x in 1:length(w.solution)){
#                     bag <- c(bag, rep(x, w.solution[x]))
#                 }
#                 
#                 # add to matrices
#                 bag.mat[itr,] <- study.num[bag] # add bag to matrix.
#                 bag.recipe[itr, ] <- w.solution  # add bag recipe to
#                 itr <- itr + 1 # update bag counter
#             }else{
#                 # if NA then move to next w limit
#                 lambda <- lambda * phi
#                 z <- z + 1
#             }
#             
#         }else{
#             # if repeat bag then move to next w limit
#             lambda <- lambda * phi # expand radius by factor of how big the minimum distance is
#             z <- z + 1
#         }
#         
#         
#         if (itr > straps | z > lim){
#             # break while loop if limit reached or we have achieved desired strap number
#             flag <- FALSE
#         }
#     }
#     print(paste0("lambda for test elec_", test.elec,"_", lambda))
#     print(paste0("d.vec for test elec_", test.elec))
#     print(d.vec)
#     
#     #### add this to be the first bag
#     
#     return( list( unique(bag.mat), unique(bag.recipe), d.vec  ) ) 
# }


####################################
# Bag Generator Function Based on distances Penalize dot products (between bag variability)
####################################

# bagGen_B.dist <- function(full.CVs, test.elec, distances, bag.size = nrow(full.CVs) - 1, straps = 5, lambda = 0.001){
#     # like bagGen.CF but takes in distances for constraints
#     # input is coordinates if inflection point coordinates are fed
#     # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
#     # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
#     # and average CV matrix without test elec and outputs a matrix of study bags
#     # straps are number to go until stopping
#     # eta is the radius
#     # 
#     
#     # import wiggle.fn
#     source("Study Strap Functions.R")
#     
#     distances <- as.vector(distances)
#     
#     K <- nrow(full.CVs)# number of studies including test elec
#     study.num <- seq(1, K)[-test.elec] # the elecs included
#     # make objective matrix A
#     A <- matrix(NA, nrow = 8, ncol = K)
#     colnames(A) <- paste0("Elec_", 1:K)
#     row.names(A) <- paste0("Coordinate_", 1:8)
#     for (i in 1:K){
#         A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
#     }
#     
#     elec.target <- A[,test.elec] # the "target" in the objective
#     
#     A <- A[,-test.elec] # remove for optimization below
#     
#     # bag matrix
#     bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
#     bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
#     colnames(bag.recipe) <- paste0("P_",study.num)
#     row.names(bag.recipe) <- paste0("Target",1:straps)
#     
#     #### Prepare Constrained Optimization ####
#     # load package
#     suppressWarnings(library(CVXR, warn.conflicts=FALSE))
#     
#     w <- Int(K - 1)
#     X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
#     
#     # minimum dist
#     objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
#     problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
#     result <- solve(problem)
#     w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
#     B <- matrix(NA, ncol = bag.size, nrow = 1) # matrix of l2 norms of differences
#     B[1,] <- w.solution
#     
#     # make bag
#     bag <- c()
#     for (x in 1:length(w.solution)){
#         bag <- c(bag, rep(x, w.solution[x]))
#     }
#     
#     bag.mat[1,] <- study.num[bag] # add bag to matrix.
#     bag.recipe[1,] <- w.solution  
#     
#     for(itr in 2:straps){
#         # minimum dist
#         d <- distances[itr]
#         lambda.norm <- lambda / (itr - 1) # update this so we don't penalize more with increasing numbers of additions 
#         objective <- Minimize(pos(sum((elec.target - X %*% w )^2) - d) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
#         problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
#         result <- solve(problem)
#         w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
#         B <- rbind(B, w.solution) # update B
#         
#         # make bag
#         bag <- c()
#         for (x in 1:length(w.solution)){
#             bag <- c(bag, rep(x, w.solution[x]))
#         }
#         # add to matrices
#         bag.mat[itr,] <- study.num[bag] # add bag to matrix.
#         bag.recipe[itr,] <- w.solution  # add bag recipe to 
#     }
#     
#     #### add this to be the first bag
#     
#     return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
# }


bagGen_B4 <- function(full.CVs, test.elec, distances, bag.size = nrow(full.CVs) - 1, 
                      straps = 20, lambda = 1, lim = 50, phi = 1.25, voltWeights = 1){
    # like bagGen_B3 but includes a min dist in the objective so penalization is all with respect to closest elec
    # input is coordinates if inflection point coordinates are fed
    # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    # straps are number to go until stopping
    # phi is the update rate
    # lim is the limit of tries
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    distances <- as.vector(distances)
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    A[1:4,] <- A[1:4,] * voltWeights # upweight voltage potentials
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = straps) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:straps)
    
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # minimum dist
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    min.dist <- sum((elec.target - X %*% w.solution )^2)
    d.vec <- c( min.dist )
    
    B <- matrix(NA, ncol = bag.size, nrow = 1) # matrix of l2 norms of differences
    B[1,] <- w.solution
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    
    bag.mat[1,] <- study.num[bag] # add bag to matrix.
    bag.recipe[1,] <- w.solution  
    # itr <- 2
    # flag <- TRUE
    # z <- 0
    
    for(itr in 2:straps){
        # minimum dist
        d <- distances[itr]
        lambda.norm <- lambda / (itr - 1)  # update this so we don't penalize more with increasing numbers of additions 
        objective <- Minimize(pos(sum((elec.target - X %*% w )^2) - d) + lambda.norm * sum(B %*% w) )  # 1/K part is to 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
        result <- solve(problem)
        w.solution <- as.vector(round(result$getValue(w))) # round since the 0s are not exactly 0
        
        if (length(w.solution) > 1){
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            
            if (!any(apply(bag.recipe, 1, function(x, want) isTRUE(all.equal(as.numeric(x),
                                                                             as.numeric(want))), w.solution))){
                # if the bag has not yet been added to the bag matrix (i.e., it is unique)
                # then add to the matrix
                
                B <- rbind(B, w.solution) # update B
                d.vec <- c(d.vec,  sum((elec.target - X %*% w.solution )^2)  )
                # make bag
                bag <- c()
                for (x in 1:length(w.solution)){
                    bag <- c(bag, rep(x, w.solution[x]))
                }
                
                # add to matrices
                bag.mat[itr,] <- study.num[bag] # add bag to matrix.
                bag.recipe[itr, ] <- w.solution  # add bag recipe to
                itr <- itr + 1 # update bag counter
            }
        }
            
        #     }else{
        #         # if NA then move to next w limit
        #         lambda <- lambda * phi
        #         z <- z + 1
        #     }
        #     
        # }else{
        #     # if repeat bag then move to next w limit
        #     lambda <- lambda * phi # expand radius by factor of how big the minimum distance is
        #     z <- z + 1
        # }
        
        # 
        # if (itr > straps | z > lim){
        #     # break while loop if limit reached or we have achieved desired strap number
        #     flag <- FALSE
        # }
    }
    print(paste0("lambda for test elec_", test.elec,"_", lambda))
    print(paste0("d.vec for test elec_", test.elec))
    print(d.vec)
    
    bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,] #eliminate NAs
    bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,] #eliminate NAs
    
    return( list( unique(bag.mat), unique(bag.recipe), d.vec  ) ) 
}

####################################
# Bag Generator Function Based on distances
####################################

bagGen.distTest <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
       # objective <- Minimize(  pos(norm2(elec.target - X %*% w ) - d  ) )  
        objective <- Minimize( pos( norm2((elec.target - X %*% w)^2) - d) + lambda * sum(w^2)  )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}




bagGen.dist3 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d)  )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}

# 
# for (i in 1:nrow(bag.mat)){
# 
#     print( sum( (X %*% bag.recipe[i,] - elec.target)^2     )   )
# }


####################################
# Bag Generator Function Based on distances
####################################

bagGen.distConvx <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Variable(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round_preserve_sum(result$getValue(w)) # round since the 0s are not exactly 0
        print(paste0("sum w.sol", sum(w.solution)))
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.distContin <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Variable(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- (result$getValue(w)) # Do not round
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


round_preserve_sum <- function(x, digits = 0) {
         up <- 10 ^ digits
         x <- x * up
        y <- floor(x)
        indices <- tail(order(x-y), round(sum(x)) - sum(y))
        y[indices] <- y[indices] + 1
        y / up
}


accCalc <- function(full, error.type = 1, bag, elecs, test.elec, study.strap.rows, 
                            sampProp = 1, test.prop = 1,  best.comp = best.comp){
    # accuracy calculator
    # takes in objects from global environment and calculates RMSE (training), OOB error metrics
    # returns matrix of errors in cols and models in rows
    # sampProp is proportion of full dataset to fit models 
    # test.prop is proportion of full dataset to test models on
    # rows are the subset of rows used to train the study strap on
    # best.comp is number of components used
    # bag is the study bag as usual
    # elecs are the elecs that can be OOB or trained on (not held out in anyway)
    # error.type is the type of error to calculate
    library(pls)
    OOB <- setdiff(elecs, unique(bag))
    colmns <- (ncol(full)-1000):ncol(full)
    X.colmns <- (ncol(full)-1000):ncol(full)
    
    results.vec <- vector(length = 4)
    names(results.vec) <- c("Train.error", "oob.elecs.merged", "oob.rows", "in.bag.elecs.oob")
    
    # fit model
    print(" acc tune model")
    fit.tune <- pcr(DA ~., data = full[study.strap.rows , colmns], ncomp = best.comp, 
                    model = FALSE, y = FALSE)
    
    ##################
    # Training RMSE
    ##################
    if( error.type == 1){
        print(" acc resids")

        results.vec[1] <- mean( (fit.tune$residuals[, 1, best.comp])^2 ) # training RMSE
    
    }else if (error.type == 2){
        ##################
        # oob.elecs.merged
        ##################
        # All OOB Rows
        # take a sub sample of all OOB rows (length equal to the length of the study strap to save memory)
        
        rows.samp <- sample( seq(1, nrow(full))[-study.strap.rows], length(study.strap.rows), 
                             replace = FALSE )  
        oob.pred <- predict(fit.tune, full[rows.samp, X.colmns], ncomp = best.comp) 
        results.vec[2] <- sqrt(mean(( full$DA[rows.samp] - oob.pred)^2)) #only one metric
        rm(oob.pred, rows.samp)
    
    }else if ( error.type == 3){
    
        #######################
        # All OOB elecs Merged
        #######################
        # Metric 3: All OOB elecs Merged
           # take a sub sample of all OOB elecs merged (length equal to the length of the study strap to save memory)
        if (length(OOB) == 0){ 
            # if there are no OOB Elecs (i.e., the bag is all the elecs) 
            # then just do OOB rows on in bag elecs
            in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
            # take the rows of in bag electrodes, but are out-of-sample
            in.out <- setdiff(in.bag.rows, study.strap.rows) 
            rows.samp <- sample( in.out, min(length(study.strap.rows), length(in.out)), replace = FALSE ) # sub sample to save memory
            oob.pred <- predict(fit.tune, full[rows.samp, X.colmns], 
                                ncomp = best.comp) # predict on model based on current OOB electrode covariates
            results.vec[3] <- sqrt(mean(  (full$DA[rows.samp] - oob.pred)^2  )) #  
            rm(oob.pred, in.bag.rows, in.out, rows.samp, fit.tune) # save memory 
        }else{
            # if there are OOB Elecs
            oob.rows <- seq(1, nrow(full))[is.element(full$Electrode, OOB)]
            rows.samp <- sample( oob.rows, min(length(study.strap.rows), length(oob.rows)), replace = FALSE )  
            oob.pred <- predict(fit.tune, full[rows.samp, X.colmns],
                                ncomp = best.comp) # predict on model based on current OOB electrode covariates
            results.vec[3] <- sqrt( mean(   (full$DA[rows.samp] - oob.pred)^2  ) )
            rm(oob.pred, rows.samp, oob.rows)
        }
    
    }else if (error.type == 4){
        ##################################
        # OOB-obs on in-bag-elecs
        ##################################
        # Metric 4: OOB - Merged all OOB observations on in bag electrodes**********
        in.bag.rows <- which(is.element(full$Electrode, bag)) # all the rows of the in bag electrodes (not necessarily rows used in Study Strap)
        # take the rows of in bag electrodes, but are out-of-sample
        in.out <- setdiff(in.bag.rows, study.strap.rows) 
        rows.samp <- sample( in.out, min(length(study.strap.rows), length(in.out)), replace = FALSE ) # sub sample to save memory
        oob.pred <- predict(fit.tune, full[rows.samp, X.colmns], 
                            ncomp = best.comp) # predict on model based on current OOB electrode covariates
        results.vec[4] <- sqrt(mean(  (full$DA[rows.samp] - oob.pred)^2  )) #  
        rm(oob.pred, in.bag.rows, in.out, rows.samp, fit.tune) # save memory 
        }
    return(results.vec[error.type])
}


bagSampler <- function(bag.mat, tune.fraction = 1, bag.size, full,
                       proportion = 1){
    # takes in a matrix of bags and returns a list of lists of vectors of rows of study straps
    # first list in returned list is the rows for full study straps
    # second list in returned list is the rows fo tuning (i.e., a subset of the above with proportion sampled denoted by tune.fraction)
    # tune.fraction is the proportion of rows to take (1 is the normal amount used for the actual ensemble models)
    
    # length of bag.mat (check to see if its just a vector or a matrix)
    if(class(bag.mat) == "integer" | class(bag.mat) == "numeric"){
        len <- 1
    }else{ 
        len <- nrow(bag.mat) 
        }

    results <- list(length = len) # bag rows for real straps
    tune.straps <- list(length = len) # bag rows for tuning, OOB, accuracy for bag selection
    
    # iterate through rows of matrix and generate bags
    for(bag.num in 1:len){
        bag <- as.vector(bag.mat[bag.num,]) # set current bag
        strap.table <- as.data.frame(table(bag))
        
        indx.tune <- c()
        indx <- c()
        # sub sample rows
        for(i in 1:nrow(strap.table)){
            
            sub.elec <- as.numeric(as.character(strap.table[i,1])) # current electrode being sampled from
            elec.indx <- which(full$Electrode == sub.elec) # rows corresponding to electrode sampled in current bag
            proportion <- as.numeric(as.character(strap.table[i,2])) # proportion of rows to sample = proportion/bag.size
            num.obs <- floor(length(elec.indx) / (bag.size)) * proportion # number of rows to sample times the number of times the electrode shows up in the current study strap bag
            #num.obs <- round( (length(elec.indx) / bag.size) * as.numeric(as.character(strap.table[i,2])) )
            elec.indx <- elec.indx[sample.int(length(elec.indx), num.obs, replace = FALSE)] #sub-sample rows
            elec.tune.indx <- sample(elec.indx, round(length(elec.indx) * tune.fraction), replace = FALSE) # sub sample rows to train classifier on
            indx <- c(indx, elec.indx ) # sample as many rows as indicated in elec.obs and add to indices vector
            indx.tune <- c(indx.tune, elec.tune.indx)
            # bag.recipe[r, which(elecs == sub.elec)] <- proportion # put the proportion in the corresponding index. which(elecs == sub.elec) is so it does not go out of range
        }
        
        results[[bag.num]] <- indx # rows for bag
        tune.straps[[bag.num]] <- indx.tune # for tuning (i.e., its just a subset of the rows in "results" according to tune.frac)
    }
    return( list(results, tune.straps) )
}


divAcc.opt <- function(alpha = 0.5, error.vec, div.mat, final.straps){
    # takes in tuning parameter alpha, error vector and divesity matrix and 
    # selects a subset of length "final.straps" that minimizes a convex combination 
    # of error and diversity
    # returns vector of indicators of whether a bag is selected
    
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    K <- length(error.vec)
    w <- Int(K)
    alpha <- alpha / K # scale this because quadratic in div mat and linear in error vec
    
    print("objective")
    # do 1/ diversity matrix since we are minimizing -- could also do cosine similarity metric
    objective <- Minimize( alpha * quad_form(w, div.mat ) + (1 - alpha) * t(error.vec) %*% w)
    final.straps <- min(final.straps, K) # in case there are fewer bags in bag.mat then final.straps
    problem <- Problem(objective, constraints = list(w >= 0, w <= 1, sum(w) == final.straps))
    print("solve objective")
    result <- solve(problem)
    w.solution <- round(result$getValue(w))
    
    return(w.solution)
}


diversity.calc <- function(mat){
    # provides list of 1) squared l2 norm of difference between vectors
    # and 2) cosine similarity measure
    mat <- mat[rowSums(!is.na(mat)) > 0,] # remove NAs
    B <- matrix(NA, nrow = nrow(mat), ncol = nrow(mat))
    for ( i in 1:nrow(mat)){
        for (j in 1:nrow(mat)){
            B[i,j] <- sum( (mat[i, ] - mat[j,])^2 )
        }
    }
    B.cos <- matrix(NA, nrow = nrow(mat), ncol = nrow(mat))
    # cosine similarity metric
    for ( i in 1:nrow(mat)){
        for (j in 1:nrow(mat)){
            B.cos[i,j] <-  ( mat[i, ] %*% mat[j,] ) / (norm(mat[i, ], "2") * norm(mat[j, ], "2"))
        }
    }
    B <- (1/B) # 1/B because we are minimizing later
    diag(B) <- 1 # delete infinite values
    
    return( list(B , B.cos ) )  
}


paramTune <- function(elecs, full, straps, bags.start, full.CVs, sprd = -2.5, best.comp, bag.size){
    # does a hold one study out CV ensemble
    # elecs are all the electrodes that are not the test electrode in the full algorithm
    # full is the full dataset
    # straps is the number of study straps that we test
    
    error.mat <- matrix(NA, ncol = straps, nrow = length(elecs))
    colnames(error.mat) <- paste0("Strap_", 1:straps)
    
    print("error.vec")
    error.mat <- matrix(NA, ncol = straps, nrow = length(elecs)) # each element is error for each bag (OOB or training error)
    
    for (test.indx in 1:length(elecs) ){
        test <- elecs[test.indx] # test here is the held out test elec (held out in the tuning sense)
        studies <- elecs[-test.indx]
        train.elecs <- elecs[-which(elecs == test)]
        # need to remake targetDist and bagGen.dist so that can take in elecs that are fair game
        print("targetDist.crossVal")
        distances <- targetDist.crossVal(full.CVs = full.CVs, test.elec = test, 
                                         straps = bags.start, sprd = sprd) # generate target study strap inflection points
        print("bagGen.dist.crossVal")
        bag.mat <- bagGen.dist.crossVal(full.CVs = full.CVs, distances = distances, 
                                        test.elec = test, bag.size = bag.size, lambda = 0, 
                                            train.elecs = train.elecs, straps = straps) # generate bag matrix
        bag.recipe <- bag.mat[[2]]
        bag.mat <- bag.mat[[1]]
        print(paste0("bag.mat pre-opt_", y))
        print(bag.mat)
        
        print(paste0("bagSampler_", y, "_CVelec_",test.indx))
        bagRows.mat <- bagSampler(bag.mat, tune.fraction = 1, bag.size = bag.size,
                                  full = full[full$Electrode != test,], proportion = 1)[[2]] # 2nd item is tune.straps rows
        
        print(paste0("accCalc.strapsTuner_", y, "_CVelec_",test.indx))
        error.mat[test.indx,] <- accCalc.strapsTuner(full = full, 
                                    error.type = error.type, bag.mat = bag.mat, elecs = train.elecs, 
                                          test.elec = test, bagRows.mat = bagRows.mat, 
                                            sampProp = 1, test.prop = 1,  best.comp = best.comp)

    }
    print("error.mat")
    print(error.mat)
    print("colmeans error.mat")
    print(colMeans(error.mat))
    best.straps <- which.min(colMeans(error.mat))
    return(best.straps)
}



######################################################
# Generate Target Distances
######################################################
targetDist.crossVal <- function(full.CVs, test.elec, straps, sprd = - 2.5){
    # elecs is the list of elecs that are used for training (does not include test.elec)
    # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
    # test.elec is the number of the electrode that is the test electrode (integer \in [K])
    # straps is the number of study straps to generate
    # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
    # that the target electrodes should be
    
    source("Study Strap Functions.R") # needed for wiggle.fn
    #------------------------------------------
    # make Coordinate matrix A
    print("convert to dataframe")
    full.CVs <- as.data.frame(full.CVs)
    print("test.elec index")
    test.elec.indx <- which(full.CVs$Electrode == test.elec)
    full.CVs <- as.matrix(full.CVs)
    
    K <- nrow(full.CVs) 
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_",1:K)
    row.names(A) <- paste0("Coordinate_",1:8)
    print("A matrix")
    for (i in 1:nrow(full.CVs)){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    # just use distance metric
    dist.vec <- vector(length = ncol(A)) #vector of distance metrics
    for (z in 1:ncol(A)){
        dist.vec[z] <- sum( (A[,z] - A[,test.elec.indx])^2 )
    }
    dist.diff <- max(dist.vec) # maximum distance metric from the test elec
    d <- 10^(seq(sprd, 0, length = straps)) # strap dist
    d[which.max(d)] <- 0
    seq.vec <-  abs(dist.diff) * d
    
    return(seq.vec)
    
}


paramTune2 <- function(elecs, full, straps, bags.start, full.CVs, sprd = -2.5, best.comp, bag.size){
    # same as paramTune but uses targetDist2 instead of targetDist so uses the minimum distance 
    # as the starting point for calculating distances
    # does a hold one study out CV ensemble
    # elecs are all the electrodes that are not the test electrode in the full algorithm
    # full is the full dataset
    # straps is the number of study straps that we test
    
    error.mat <- matrix(NA, ncol = straps, nrow = length(elecs))
    colnames(error.mat) <- paste0("Strap_", 1:straps)
    
    print("error.vec")
    error.mat <- matrix(NA, ncol = straps, nrow = length(elecs)) # each element is error for each bag (OOB or training error)
    
    for (test.indx in 1:length(elecs) ){
        test <- elecs[test.indx] # test here is the held out test elec (held out in the tuning sense)
        studies <- elecs[-test.indx]
        train.elecs <- elecs[-which(elecs == test)]
        # need to remake targetDist and bagGen.dist so that can take in elecs that are fair game
        print("targetDist.crossVal")
        distances <- targetDist.crossVal2(full.CVs = full.CVs, test.elec = test, 
                                         straps = bags.start, sprd = sprd, bag.size = bag.size - 1) # generate target study strap inflection points
        # bag.size = bag.size - 1 since one more is held out
        print("bagGen.dist.crossVal")
        bag.mat <- bagGen.dist.crossVal(full.CVs = full.CVs, distances = distances, 
                                        test.elec = test, bag.size = bag.size, lambda = 0, 
                                        train.elecs = train.elecs, straps = straps) # generate bag matrix
        bag.recipe <- bag.mat[[2]]
        bag.mat <- bag.mat[[1]]
        print(paste0("bag.mat pre-opt_", y))
        print(bag.mat)
        
        print(paste0("bagSampler_", y, "_CVelec_",test.indx))
        bagRows.mat <- bagSampler(bag.mat, tune.fraction = 1, bag.size = bag.size,
                                  full = full[full$Electrode != test,], proportion = 1)[[2]] # 2nd item is tune.straps rows
        
        print(paste0("accCalc.strapsTuner_", y, "_CVelec_",test.indx))
        error.mat[test.indx,] <- accCalc.strapsTuner(full = full, 
                                                     error.type = error.type, bag.mat = bag.mat, elecs = train.elecs, 
                                                     test.elec = test, bagRows.mat = bagRows.mat, 
                                                     sampProp = 1, test.prop = 1,  best.comp = best.comp)
        
    }
    print("error.mat")
    print(error.mat)
    print("colmeans error.mat")
    print(colMeans(error.mat))
    best.straps <- which.min(colMeans(error.mat))
    return(best.straps)
}


######################################################
# Generate Target Distances
######################################################
targetDist.crossVal2 <- function(full.CVs, test.elec, straps, sprd = - 2.5, bag.size = 13){
    # like targetDist.crossVal but we start distances at min.dist
    # elecs is the list of elecs that are used for training (does not include test.elec)
    # fullCVs is a k x 1000 matrix with all the average CVs INCLUDING test elec
    # test.elec is the number of the electrode that is the test electrode (integer \in [K])
    # straps is the number of study straps to generate
    # sprd is a spread parameter indexing spread of elecs for dist method 2 with support \in (-\infty, 0)
    # that the target electrodes should be
    
    source("Study Strap Functions.R") # needed for wiggle.fn
    #------------------------------------------
    # make Coordinate matrix A
    print("convert to dataframe")
    full.CVs <- as.data.frame(full.CVs)
    print("test.elec index")
    test.elec.indx <- which(full.CVs$Electrode == test.elec)
    full.CVs <- as.matrix(full.CVs)
    
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    K <- nrow(full.CVs) 
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_",1:K)
    row.names(A) <- paste0("Coordinate_",1:8)
    print("A matrix")
    for (i in 1:nrow(full.CVs)){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    # just use distance metric
    dist.vec <- vector(length = ncol(A)) #vector of distance metrics
    for (z in 1:ncol(A)){
        dist.vec[z] <- sum( (A[,z] - A[,test.elec.indx])^2 )
    }
    
    elec.target <- A[,test.elec.indx] # the "target" in the objective
    A <- A[,-test.elec.indx] # remove for optimization below
    
    w <- Int(bag.size)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    # Constrained Optimization (Integer Least Squares)
    problem <- Problem(objective, constraints = list(w >= 0, 
                                                     sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    min.dist <- sum((elec.target - X %*% w.solution )^2)
    
    dist.diff <- abs(max(dist.vec) - min.dist) # difference between maximum distance metric from closest feasible study strap to  test elec
    
    d <- 10^(seq(sprd, 0, length = straps)) # strap dist
    d[which.max(d)] <- 0 # make the furthest distance a 0 since the study strap would be just the furthest electrode
    seq.vec <-  min.dist + abs(dist.diff) * d
    
    return(seq.vec)
    
}


####################################
# Bag Generator Function Based on distances
####################################

bagGen.dist.crossVal <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, 
                                        lambda = 0, train.elecs, straps){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    full.CVs <- as.data.frame(full.CVs)
    print("bagGen.dist.crossVal test elec indx")
    test.elec.indx <- which(full.CVs$Electrode == test.elec) # index of test elec necessary for proper indexing
    bag.size <- length(train.elecs) # do this since bag.size is for the global ensemble not cross validation
    full.CVs <- as.matrix(full.CVs)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- train.elecs
    # make objective matrix A
    print("A matrix")
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    print("elec.target")
    elec.target <- A[,test.elec.indx] # the "target" in the objective
    
    A <- A[,-test.elec.indx] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # target electrode 
    print("for loop Optimization step")
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        # objective <- Minimize(  square(d - norm2(elec.target - X %*% w )   ) )  
        objective <- Minimize( pos( sum((elec.target - X %*% w)^2) - d) )
        #objective <- Minimize( sum((elec.target - X %*% w )^2))     
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        
        # make bag
        bag <- c()
        for (x in 1:length(w.solution)){
            bag <- c(bag, rep(x, w.solution[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        print("add to bag")
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w.solution  # add bag recipe to 
    }
    
    if(straps < nrow(bag.mat)){print("error: straps < nrow(bag.mat)")}
    return( list( unique(bag.mat)[1:min(straps, nrow(bag.mat)), ], 
                        unique(bag.recipe)[1:min(straps, nrow(bag.mat)),]  ) ) 
}



accCalc.strapsTuner <- function(full, error.type = 1, bag.mat, elecs, test.elec, bagRows.mat, 
                    sampProp = 1, test.prop = 1,  best.comp){
    # accuracy calculator for different numbers of straps
    # takes in objects from global environment and calculates RMSE (training), OOB error metrics
    # returns matrix of errors in cols and models in rows
    # study.strap.rows is a list of rows here
    # sampProp is proportion of full dataset to fit models 
    # test.prop is proportion of full dataset to test models on
    # rows are the subset of rows used to train the study strap on
    # best.comp is number of components used
    # bag is the study bag as usual
    # elecs are the elecs that can be OOB or trained on (not held out in anyway)
    # error.type is the type of error to calculate
    library(pls)
    
    pred.mat <- matrix(NA, ncol = nrow(bag.mat), nrow = nrow(full[full$Electrode == test.elec, ])) # matrix of predictions
    colnames(pred.mat) <- paste0("Strap", 1:ncol(pred.mat)) 
    
    error.vec <- vector(length = ncol(pred.mat)) # vector of errors of each progressive extra strap
    
    for (bag.num in 1:nrow(bag.mat)){
    
        # current strap
        print(paste0("accTune strap#_", bag.num))
        bag <- bag.mat[bag.num,]
        #study.strap.rows <- bagRows.mat[bag.num,]
        study.strap.rows <- bagRows.mat[[bag.num]]
        print(paste0("study.strap.rows strap#_", bag.num, "_nrows_", length(study.strap.rows)))
        print(paste0("test elec rows#_", bag.num, "_nrows_", 
                     nrow(full[full$Electrode == test.elec, ])))
        
        OOB <- setdiff(elecs, unique(bag))
        colmns <- (ncol(full)-1000):ncol(full)
        X.colmns <- (ncol(full)-999):ncol(full)
        
        # fit model
        print(" acc tune model")
        fit.tune <- pcr(DA ~., data = full[study.strap.rows , colmns], ncomp = best.comp, 
                        model = FALSE, y = FALSE)
        
        pred.mat[,bag.num] <- predict(fit.tune, full[full$Electrode == test.elec, X.colmns], ncomp = best.comp)
        print("head pred mat")
        print(head(pred.mat))
        if (bag.num == 1){
            # to avoid problems with rowMeans
            error.vec[bag.num]  <- sqrt( mean(  (  (pred.mat[,1]) - 
                                                       full$DA[full$Electrode == test.elec] )^2  )  ) # RMSE of adding each additional error vector
        }else{
            error.vec[bag.num]  <- sqrt( mean(  (  rowMeans(pred.mat[,1:bag.num]) - 
                                                       full$DA[full$Electrode == test.elec] )^2  )  ) # RMSE of adding each additional error vector
        }
        print(paste0("error.vec_bag.number_", bag.num))
        print(error.vec[bag.num])
        }
    
    return(error.vec)
}

randomBag <- function(full.CVs, distances = distances, test.elec = y, bag.size = bag.size, straps = 20){
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    elecs <- sort(study.num)
    K <- nrow(full.CVs)# number of studies including test elec
    #elecs <- elecs[-which(elecs == test.elec)] # remove test elec
    
    # bag matrix
    bag.recipe <- matrix(0, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", elecs)
    row.names(bag.recipe) <- paste0("Target", 1:straps)
    
    print("make bag.mat")
    elems <- bag.size * straps
    bag.mat <- matrix( sample(elecs, elems, replace = TRUE), ncol = bag.size, nrow = straps)
    
    print("bag.recipe iterate")
    for (x in 1:nrow(bag.mat)){
        bag <- bag.mat[x,]
        
        # iterate through bags
        for (i in 1:length(elecs)){
            # iterate through elecs in bags
            bag.recipe[x, i] <- sum(bag == elecs[i])
        }
    }
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) )
    
}


bagGen.DCA <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, 
                            lambda = 0, tol = 10, percent = 0.1){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix)
    # and average CV matrix without test elec and outputs a matrix of study bags
    # includes a penalization term to promote diversity that is a difference of convex functions
    # (l2 norm of diference in bags)
    # tol is the number of iterations before convergence until restart starting vector 
    # percent is the percentage of adjacent distances that the maximum promotion can achieve

    distances <- sort(as.vector(distances))

    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }

    elec.target <- A[,test.elec] # the "target" in the objective

    A <- A[,-test.elec] # remove for optimization below

    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))

    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    B <- matrix(NA, nrow = 1, ncol = bag.size)
    objective <- Minimize( sum((elec.target - X %*% w)^2)  )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    B[1,] <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    
    # B[1,] <- sample(1:bag.size, bag.size, replace = TRUE) #toy example
    
    # target electrode

    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]

        # function that stays constant for each iteration
        g1.fn <- function(X, w, elec.target, d){
           pos( sum((elec.target - X %*% w)^2) - d)
        }

        cnvg <- 1 # starting value
        w1 <- rep(1, K - 1) # starting value
        w0 <- w1 # start with w0 as  previous iteration value
        DCA <- 0 # added in
        counter <- 0 # number of iterations before convergence
            for ( y in 1:nrow(B)){
                # iterate through all previous solutions
                #DCA <-  sum( (B[y,] - w1)^2 ) +  2 * t( w1 - B[y,]) %*% (w - w1)  
                DCA <-  DCA + sum( (B[y,] - w1)^2 ) +  2 * t( w1 - B[y,]) %*% (w - w1)
            }
            lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
            fn.itr <- g1.fn(X, w, elec.target, d) - lambda * DCA # 
            while (cnvg > 0.1){
                
                counter <- counter + 1
                w0 <- w1
                
                if(counter > tol){
                    # if it is not convering try a new starting vector
                    counter <- 0
                    w0 <- rand_vect(length(study.num), bag.size) # random vector of starting value
                }

                objective <- Minimize( fn.itr  )
                #objective <- Minimize( sum((elec.target - X %*% w )^2))
                problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
                result <- solve(problem)
                w1 <- round(result$getValue(w)) # round since the 0s are not exactly 0
                print(w1)
                
                # update iteration
                DCA <- 0 # erase DCA so we can remake it
                for ( y in 1:nrow(B)){
                    # iterate through all previous solutions
                    DCA <- (DCA + sum( (B[y,] - w1)^2 ) +  2 * t( w1 - B[y,]) %*% (w - w1)   )  
                }
                
                lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
                fn.itr <- g1.fn(X, w, elec.target, d) - lambda * DCA # 
                
                 cnvg <- sum ( (w0 - w1)^2  ) # difference in l2 norm between previous iteration and current one
            }
            print("DCA complete")
            # after convergence, update B to include new solution
            B <- rbind(B, t(w1))
            
            # make bag
            bag <- c()
            for (x in 1:length(w1)){
                bag <- c(bag, rep(x, w1[x]))
            }
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            bag.mat[i,] <- study.num[bag] # add bag to matrix.
            bag.recipe[i, ] <- w1  # add bag recipe to
   
             } # double check placement of this bracket
    
    bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,]
    bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,]
    return( list( unique(bag.mat), unique(bag.recipe)  ) )
}

#  for(i in 1:nrow(bag.mat)){
#     # print( sum( (wiggle.fn(colMeans(full.CVs[bag.mat[i,],-c(1,2)])) - elec.target)^2  )  )
#    print( sum( (rowMeans( D[,bag.mat[i,] ]) - elec.target)^2)  )
#      
# }


lambda.gen <- function(distances, indx, percent, bag.size){
    # calculates lambda that is a percentage of surrounding distances
    
    if (indx == 1){
        # difference between current distance and next distance
        # square bag.size because at most vector can be bag.size^2
        # divide by index because we are summing across all the previous bags
        lambda <- ( distances[2] - distances[1] ) * percent / (bag.size^2)
    }else if (indx < length(distances)){
        # average difference between current distance and surrounding distances 
        lambda <-  mean ( c( abs(distances[indx + 1] - distances[indx]), 
                             abs(distances[indx - 1] - distances[indx]) )) * percent / (bag.size^2 * indx)
    }else if (indx == length(distances)){
        # difference between current distance and next distance
        lambda <- abs( distances[length(distances)] - distances[length(distances) - 1] ) * percent / (bag.size^2 * indx)
    }
    return(lambda)
}

# print("opt minimum")
# w <- Int(K - 1)
# X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
# W <- w
# W_k <- w1
# #prepare
# if(nrow(B) > 1){
#     for ( i in 1:(nrow(B)-1)){
#         W <- hstack(W,w)
#         W_k <- hstack(W_k, w1)
#         
#     }
# }
# 
# 
# dim(X)
# # optimize
# f2 <- function(X,W) {sum(cvxr_norm(W - t(B), 2, axis = 2))}
# f3 <- function(X,W,w,B,w1){-2 * sum( (B - t(W_k)) %*% (w - w1))}
# objective <- Minimize(g1.fn(X, w, elec.target, d)  - lambda * (f2(X,W_k) +  f3(X,W,w,B,w1) )     )          #   sum((elec.target - X %*% w)^2) )
# #objective <- Minimize(g1.fn(X, w, elec.target, d)  - lambda * (f2(X,W_k)   ))          #   sum((elec.target - X %*% w)^2) )
# 
# problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
# result <- solve(problem)
# w.opt <- round(result$getValue(w)) # round since the 0s are not exactly 0





pathGen <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, paths = 5, 
                    converge.lim = 50000, prox.straps = 0, write_dist = FALSE){
    # provides distance vector based upon fitting an exponential curve to random sample of paths
    # full.CVs is matrix of average features (not sure what else it includes)
    # test.elec is the study number of the target 
    # paths is the number of random paths to sample
    # converge.lim is the number of random samples to take before ending the path
    # prox.straps - proximal straps -- the number of extra straps (in addition to the optimum) very close to the target that we include PER PATH in the model fitting process
    # prox.straps also determines the total number of straps for each electrode-- is just the average straps/path + prox.straps is the total straps to be fit
    # write_dist indicator of whether to write a csv of the distances
    
    setwd("/n/home12/gloewinger")
    source("Study Strap Functions.R")
    total.proximal <- prox.straps * paths # total number of proximal straps to optimize to include in model fitting step
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    full.CVs <- full.CVs[,-c(1,2)] # remove all but covariates
    
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    
    target.elec <- full.CVs[test.elec,]
    full.CVs[test.elec,] <- NA
    CVs <- full.CVs[-test.elec,] # remove test elec
    elecs <- study.num
    
    ###############################################
    # Optimization Step prior to curve fitting
    ###############################################
    #### Prepare Constrained Optimization ####
    print("opt minimum")
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints

    # optimize
    objective <- Minimize( sum((elec.target - X %*% w)^2) )
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.opt <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # make bag
    bag <- c()
    for (x in 1:length(w.opt)){
        bag <- c(bag, rep(x, w.opt[x]))
    }
    bag <- study.num[bag]
    
    # calculate distance
    print(" minimum distance")
    feats <- wiggle.fn( colMeans( full.CVs[bag, ]) ) # calculate inflection points
    dist.opt <- sum(  (feats - elec.target)^2  ) # squared Euclidean Distance between inflection points
    
    # accept reject step
    counter <- 0 # number of samples before acceptance
    strap <- 1 # strap counter within path
    total.straps <- 1 # the total number of straps across all paths
    thresh <-  5e+14 # restart value: arbitraily large starting threshold
    threshold <- thresh # start at restart value
    distance.matrix <- matrix(NA, ncol = 5, nrow = 500) # arbitrarily large number of rows
    colnames(distance.matrix) <- c("Strap", "Distance", "Count", "Path", "Time")
    path.straps <- vector(length = paths) # total number of straps in each path
    for( path in 1:paths){
        
        timeStart <- Sys.time() 
        print(paste0("path ", path))
        # iterate through number of paths
        while (counter < converge.lim){ #while loop for restrapping: breaks for either acceptance (sim) or accpetance counter limit (converge.lim)
            
            # Choose electrodes to be in bag:

            bag <- sample(elecs, bag.size, replace = TRUE) # sample bag
            feats <- wiggle.fn( colMeans( full.CVs[bag, ]) ) # calculate inflection points
            current.dist <- sum(  (feats - elec.target)^2  ) # squared Euclidean Distance between inflection points
            
            if(current.dist < threshold ){
                # add distance to matrix if accepted

                distance.matrix[total.straps, ] <- c(strap, current.dist, counter, path, NA) # add to matrix

                counter <- 0 # reset counter
                strap <- strap + 1 # increase straps of current path
                total.straps <- total.straps + 1 # increase total straps
                threshold <- current.dist
            }else{
                # add counter
                counter <- counter + 1
            }
            
        }
        print(paste0("dist matrix test.elec ", test.elec))
        print(distance.matrix)
        
        # add extra values close to optimum
        if (prox.straps > 0){
            # check there are more proximal straps than just the optimum
            d <- seq(dist.opt, current.dist, length = prox.straps + 2) # make a sequence of evenly spaced distances between the optimum and the last distance of the current path
            d <- d[-c(1,length(d))] # remove the first and last elements so does not include limit values
            
            for (i in 1:length(d)){
                d.i <- d[i] # current dist
                # optimization step
                w <- Int(K - 1)
                objective <- Minimize( pos(  sum((elec.target - X %*% w)^2) - d.i) )
                problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
                result <- solve(problem)
                w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
                
                # make bag
                bag <- c()
                for (x in 1:length(w.solution)){
                    bag <- c(bag, rep(x, w.solution[x]))
                }
                bag <- study.num[bag]
                # calculate distance of optimized bag
                feats <- wiggle.fn( colMeans( full.CVs[bag, ]) ) # calculate inflection points
                dist.opt <- sum(  (feats - elec.target)^2  ) # squared Euclidean Distance between inflection points
                
                # add distances to matrix - counter = 0 so we know which is from optimization
                distance.matrix[total.straps, ] <- c(strap, current.dist, 0, path, NA) # add to matrix
                strap <- strap + 1
                total.straps <- total.straps + 1
            }
        }
        
        timeEnd <- Sys.time()      
        path.time <- as.numeric( difftime(timeEnd, timeStart, units='mins') )
        print(paste0("path ", path, "_elec_", test.elec, "_time_", path.time))
        # add optimum - counter = 0 so we know which is from optimization
        distance.matrix[total.straps, ] <- c(strap, dist.opt, 0, path, path.time) # for the last distance in every path, add the optimum distance
        
        path.straps[path] <- strap # save the total number of straps in this path
        strap <- 1 # reset strap counter for this bag
        threshold <- thresh # restart threshold
        counter <- 1
        total.straps <- total.straps + 1
    }
    
    if (write_dist == TRUE){
        # write distance matrix
        distance.matrix <- distance.matrix[rowSums(!is.na(distance.matrix)) > 0,] #eliminate NAs
        setwd("/n/home12/gloewinger/path_analysis")
        write.csv(distance.matrix, paste0("dist_mat_elec)_", test.elec))
        setwd("/n/home12/gloewinger")
    }
    
    # curve fitting
    distance.matrix <- as.data.frame( distance.matrix[rowSums(!is.na(distance.matrix)) > 0,] ) #eliminate NAs
    fit <- lm( log(Distance) ~ Strap, data = distance.matrix ) # fit log linear model 
    distances <- exp (  predict(fit, data.frame( Strap = 1:round(mean(path.straps))) )  ) # target distances are fitted values of curve for average straps
    
    return( distances ) 

}


# generates random vector that sums to M and is of length N
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
    vec <- rnorm(N, M/N, sd)
    if (abs(sum(vec)) < 0.01) vec <- vec + 1
    vec <- round(vec / sum(vec) * M)
    deviation <- M - sum(vec)
    for (. in seq_len(abs(deviation))) {
        vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
    }
    if (pos.only) while (any(vec < 0)) {
        negs <- vec < 0
        pos  <- vec > 0
        vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
        vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
    }
    vec
}





medianPath <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, paths = 5, 
                    converge.lim = 50000, prox.straps = 0, write_dist = FALSE, straps = 500){
    # provides bags from AR that is median of all those accepted
    # full.CVs is matrix of average features (not sure what else it includes)
    # test.elec is the study number of the target 
    # paths is the number of random paths to sample
    # converge.lim is the number of random samples to take before ending the path
    # prox.straps - proximal straps -- the number of extra straps (in addition to the optimum) very close to the target that we include PER PATH in the model fitting process
    # prox.straps also determines the total number of straps for each electrode-- is just the average straps/path + prox.straps is the total straps to be fit
    # write_dist indicator of whether to write a csv of the distances
    
    setwd("/n/home12/gloewinger")
    source("Study Strap Functions.R")
    total.proximal <- prox.straps * paths # total number of proximal straps to optimize to include in model fitting step
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    full.CVs <- full.CVs[,-c(1,2)] # remove all but covariates
    
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    target.elec <- full.CVs[test.elec,]
    full.CVs[test.elec,] <- NA
    CVs <- full.CVs[-test.elec,] # remove test elec
    elecs <- study.num
    
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size + 1, nrow = straps) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = straps) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))

    ###############################################
    # Optimization Step prior to curve fitting
    ###############################################
    #### Prepare Constrained Optimization ####
     # accept reject step
    counter <- 0 # number of samples before acceptance
    strap <- 1 # strap counter within path
    total.straps <- 1 # the total number of straps across all paths
    thresh <-  5e+14 # restart value: arbitraily large starting threshold
    threshold <- thresh # start at restart value
    distance.matrix <- matrix(NA, ncol = 5, nrow = 500) # arbitrarily large number of rows
    colnames(distance.matrix) <- c("Strap", "Distance", "Count", "Path", "Time")
    path.straps <- vector(length = paths) # total number of straps in each path
    for( path in 1:paths){
        
        timeStart <- Sys.time() 
        print(paste0("path ", path))
        # iterate through number of paths
        while (counter < converge.lim){ #while loop for restrapping: breaks for either acceptance (sim) or accpetance counter limit (converge.lim)
            
            # Choose electrodes to be in bag:
            
            bag <- sample(elecs, bag.size, replace = TRUE) # sample bag
            feats <- wiggle.fn( colMeans( full.CVs[bag, ]) ) # calculate inflection points
            current.dist <- sum(  (feats - elec.target)^2  ) # squared Euclidean Distance between inflection points
            
            if(current.dist < threshold ){
                # add distance to matrix if accepted
                
                distance.matrix[total.straps, ] <- c(strap, current.dist, counter, path, NA) # add to matrix
                bag.mat[total.straps,] <- c(path, bag)
                
                # make bag recpie
                for (study in elecs){
                    study.indx <- which(elecs == study) # study index
                    bag.recipe[total.straps, study.indx] <- sum(I(bag==study)) # count number
                }
                
                counter <- 0 # reset counter
                strap <- strap + 1 # increase straps of current path
                total.straps <- total.straps + 1 # increase total straps
                threshold <- current.dist
            }else{
                # add counter
                counter <- counter + 1
            }
            
        }
        
        timeEnd <- Sys.time()      
        path.time <- as.numeric( difftime(timeEnd, timeStart, units='mins') )
        print(paste0("path ", path, "_elec_", test.elec, "_time_", path.time))
        # add optimum - counter = 0 so we know which is from optimization
        distance.matrix[total.straps, ] <- c(strap, dist.opt, 0, path, path.time) # for the last distance in every path, add the optimum distance
        
        path.straps[path] <- strap # save the total number of straps in this path
        strap <- 1 # reset strap counter for this bag
        threshold <- thresh # restart threshold
        counter <- 1
        total.straps <- total.straps + 1
    }
    
    if (write_dist == TRUE){
        # write distance matrix
        distance.matrix <- distance.matrix[rowSums(!is.na(distance.matrix)) > 0,] #eliminate NAs
        setwd("/n/home12/gloewinger/path_analysis")
        write.csv(distance.matrix, paste0("dist_mat_elec)_", test.elec))
        setwd("/n/home12/gloewinger")
    }
    
    #### find median ###
    # find minimum number of straps so all can be compared
    min.straps <- min(path.straps) # find minimum number of straps so all can be compared
    dist.vec <- vector(length = paths)
    for(i in paths){
        dists <- distance.matrix$Distance[distance.matrix$Path == i] # vector of distances of each path
        dist.vec[i] <- sum(dists[1:min.straps]) # sum 
    }
    strap.median <- median(dist.vec) # median distances
    path.select <- which.min( sum( (strap.median - dist.vec)^2 )) # find the median path -- I use the l2 norm here in case paths is not an odd number then we wouldnt have equality to the median
    
    # just return rows of those from the median path
    bag.recipe <- bag.recipe[bag.mat[,1] == path.select,] # intentionally use the bag.mat rows since bag recipe doesnt have it
    bag.mat <- bag.mat[bag.mat[,1] == path.select, -1]
    
    return( list(unique(bag.mat), unique(bag.recipe)) ) 
    
}



bagGen.Dist_DCA <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, 
                       lambda = 0, tol = 10, percent = 0.1){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix)
    # and average CV matrix without test elec and outputs a matrix of study bags
    # includes a penalization term to promote diversity that is a difference of convex functions
    # (l2 norm of diference in bags)
    # tol is the number of iterations before convergence until restart starting vector 
    # percent is the percentage of adjacent distances that the maximum promotion can achieve
    
    distances <- sort(as.vector(distances))
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    #w <- Variable(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    B <- matrix(NA, nrow = 1, ncol = bag.size)
    objective <- Minimize( sum((elec.target - X %*% w)^2)  )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    B[1,] <- round(result$getValue(w)) # round since the 0s are not exactly 0
    tol <- 30
    
    # B[1,] <- sample(1:bag.size, bag.size, replace = TRUE) #toy example
    
    # target electrode
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        
        # function that stays constant for each iteration
        g.fn <- function(X, w, elec.target, d){
            sum((elec.target - X %*% w)^2)^2 + d^2
        }
        
        cnvg <- 1 # starting value
        w1 <- rep(1, K - 1) # starting value
        w0 <- w1 # start with w0 as  previous iteration value
        h.fn <- function(X, w, w1, elec.target, d){
            2 * d * sum((X %*% w1 - elec.target)^2) + 
            t(4 * d * t(X) %*% (X %*% w1 - elec.target) ) %*% (w - w1)
            }
        counter <- 0 # number of iterations before convergence
        #for ( y in 1:nrow(B)){
            # iterate through all previous solutions
            #DCA <-  sum( (B[y,] - w1)^2 ) +  2 * t( w1 - B[y,]) %*% (w - w1)  
            #DCA <-  DCA + sum( (B[y,] - w1)^2 ) +  2 * t( w1 - B[y,]) %*% (w - w1)
        #}
        #lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
        fn.itr <- g.fn(X, w, elec.target, d) -  h.fn(X, w, w1, elec.target, d) # 
        while (cnvg > 0.1){
            
            counter <- counter + 1
            w0 <- w1
            
            if(counter > tol){
                # if it is not convering try a new starting vector
                print("reset")
                counter <- 0
                w1 <- rand_vect(length(study.num), bag.size) # random vector of starting value
            }
            
            objective <- Minimize( fn.itr  )
            #objective <- Minimize( sum((elec.target - X %*% w )^2))
            problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
            result <- solve(problem)
            w1 <- round(result$getValue(w)) # round since the 0s are not exactly 0
            print(w1)
            
            # update iteration
            
            
            #lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
            fn.itr <- g.fn(X, w, elec.target, d) - h.fn(X, w, w1, elec.target, d) # 
            
            cnvg <- sum ( (w0 - w1)^2  ) # difference in l2 norm between previous iteration and current one
        }
        print("DCA complete")
        # after convergence, update B to include new solution
        B <- rbind(B, t(w1))
        
        # make bag
        bag <- c()
        for (x in 1:length(w1)){
            bag <- c(bag, rep(x, w1[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w1  # add bag recipe to
        
    } # double check placement of this bracket
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) )
}

# for(i in 1:nrow(bag.mat)){
#     # print( sum( (wiggle.fn(colMeans(full.CVs[bag.mat[i,],-c(1,2)])) - elec.target)^2  )  )
#     print( sum( (rowMeans( D[,bag.mat[i,] ]) - elec.target)^2)  )
#     
# }



bagGen.Cielo <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    distance_vector <- vector(length = K)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- as.matrix( A * (1/bag.size) )# Objective Matrix X reformulated so we can use integer constraints
    objective <- Minimize(sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    B <- matrix(w.solution, nrow = 1)
    w.min <- w.solution
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        objective <- Minimize(  sum(B %*% w) )    # maybe add a function that heavily penalizes when sum((X %*% w - elec.target )^2) >= d ) 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size, #,
                                                         sum((X %*% w - elec.target )^2) <= d ) ) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        print(w.solution)
        if (length(w.solution) > 1){
            B <- rbind(B, t(w.solution))
            # make bag
            bag <- c()
            for (x in 1:length(w.solution)){
                bag <- c(bag, rep(x, w.solution[x]))
            }
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            bag.mat[i,] <- study.num[bag] # add bag to matrix.
            bag.recipe[i, ] <- w.solution  # add bag recipe to 
        }
        
    }    
    # add minumum
    bag.recipe <- rbind(bag.recipe, as.vector(w.min))
    bag <- c()
    
    for (x in 1:length(w.min)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    bag.mat <- rbind(bag.mat, as.vector(bag))
    
    # eliminate redundant bags
    bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,]
    bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,]
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


### This function worked well -- show to Dr. Rahul Mazumder

bagGen.Cielo_norm <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    distance_vector <- vector(length = K)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- as.matrix( A * (1/bag.size) )# Objective Matrix X reformulated so we can use integer constraints
    objective <- Minimize(sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    B <- matrix(w.solution, nrow = 1)
    w.min <- w.solution
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        objective <- Minimize(  sum(B %*% w) / nrow(B) + sum(w^2) )    # maybe add a function that heavily penalizes when sum((X %*% w - elec.target )^2) >= d ) 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size, #,
                                                         sum((X %*% w - elec.target )^2) <= d ) ) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        print(w.solution)
        if (length(w.solution) > 1){
            B <- rbind(B, t(w.solution))
            # make bag
            bag <- c()
            for (x in 1:length(w.solution)){
                bag <- c(bag, rep(x, w.solution[x]))
            }
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            bag.mat[i,] <- study.num[bag] # add bag to matrix.
            bag.recipe[i, ] <- w.solution  # add bag recipe to 
        }
        
    }
    
    # add minumum
    bag.recipe <- rbind(bag.recipe, as.vector(w.min))
    bag <- c()
    
    for (x in 1:length(w.min)){
        bag <- c(bag, rep(x, w.solution[x]))
        }
    bag.mat <- rbind(bag.mat, as.vector(bag))
    
    # eliminate redundant bags
    bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,]
    bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,]
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}


bagGen.Minimize.norm <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, lambda = 0){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    
    distances <- as.vector(distances)
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    distance_vector <- vector(length = K)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_",study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- as.matrix( A * (1/bag.size) )# Objective Matrix X reformulated so we can use integer constraints
    objective <- Minimize(sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    B <- matrix(w.solution, nrow = 1)
    w.min <- w.solution
    # target electrode 
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        objective <- Minimize(  sum(w^2) )    # maybe add a function that heavily penalizes when sum((X %*% w - elec.target )^2) >= d ) 
        problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size, #,
                                                         sum((X %*% w - elec.target )^2) <= d ) ) #sum(w) == bag.size)) , sum(w^2) <= 100
        result <- solve(problem)
        w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
        print(w.solution)
        if (length(w.solution) > 1){
            B <- rbind(B, t(w.solution))
            # make bag
            bag <- c()
            for (x in 1:length(w.solution)){
                bag <- c(bag, rep(x, w.solution[x]))
            }
            #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
            # so need to use list of elecs (study num) and index them by bag
            bag.mat[i,] <- study.num[bag] # add bag to matrix.
            bag.recipe[i, ] <- w.solution  # add bag recipe to 
        }
        
    }
    
    # add minumum
    bag.recipe <- rbind(bag.recipe, as.vector(w.min))
    bag <- c()
    
    for (x in 1:length(w.min)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    bag.mat <- rbind(bag.mat, as.vector(bag))
    
    # eliminate redundant bags
    bag.mat <- bag.mat[rowSums(!is.na(bag.mat)) > 0,]
    bag.recipe <- bag.recipe[rowSums(!is.na(bag.recipe)) > 0,]
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) ) 
}



#### Rahul Mazumder method
bagGen.Dist_DCA2 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, 
                            lambda = 0, tol = 10, percent = 0.1){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix)
    # and average CV matrix without test elec and outputs a matrix of study bags
    # includes a penalization term to promote diversity that is a difference of convex functions
    # (l2 norm of diference in bags)
    # tol is the number of iterations before convergence until restart starting vector 
    # percent is the percentage of adjacent distances that the maximum promotion can achieve
    
    distances <- sort(as.vector(distances))
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    #w <- Variable(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    B <- matrix(NA, nrow = 1, ncol = bag.size)
    objective <- Minimize( sum((elec.target - X %*% w)^2)  )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    B[1,] <- round(result$getValue(w)) # round since the 0s are not exactly 0
    tol <- 30
    
    # B[1,] <- sample(1:bag.size, bag.size, replace = TRUE) #toy example
    
    # target electrode
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        
        # function that stays constant for each iteration
        g.fn <- function(X, w, elec.target, d){
            sum((elec.target - X %*% w)^2)^2 + d^2
        }
        
        cnvg <- 1 # starting value
        w1 <- rep(1, K - 1) # starting value
        w0 <- w1 # start with w0 as  previous iteration value
        h.fn <- function(X, w, w1, elec.target, d){
            2 * d * sum((X %*% w1 - elec.target)^2) + 
                t(4 * d * t(X) %*% (X %*% w1 - elec.target) ) %*% (w - w1)
        }
        counter <- 0 # number of iterations before convergence
        
        fn.itr <- sum_squares(elec.target - X %*% w)^2 + d^2 -  ( 2 * d * sum_squares(X %*% w1 - elec.target) + 
                                                      t(4 * d * t(X) %*% (X %*% w1 - elec.target) ) %*% (w - w1)) # 
        while (cnvg > 0.1){
            
            counter <- counter + 1
            w0 <- w1
            
            if(counter > tol){
                # if it is not convering try a new starting vector
                print("reset")
                counter <- 0
                w1 <- rand_vect(length(study.num), bag.size) # random vector of starting value
            }
            
            objective <- Minimize( fn.itr  )
            #objective <- Minimize( sum((elec.target - X %*% w )^2))
            problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
            result <- solve(problem)
            w1 <- round(result$getValue(w)) # round since the 0s are not exactly 0
            print(w1)
            
            # update iteration
            
            
            #lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
            fn.itr <- g.fn(X, w, elec.target, d) - h.fn(X, w, w1, elec.target, d) # 
            
            cnvg <- sum ( (w0 - w1)^2  ) # difference in l2 norm between previous iteration and current one
        }
        print("DCA complete")
        # after convergence, update B to include new solution
        B <- rbind(B, t(w1))
        
        # make bag
        bag <- c()
        for (x in 1:length(w1)){
            bag <- c(bag, rep(x, w1[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w1  # add bag recipe to
        
    } # double check placement of this bracket
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) )
}


bagGen.Dist_DCA3 <- function(full.CVs, distances, test.elec, bag.size = nrow(full.CVs) - 1, 
                             lambda = 0, tol = 10, percent = 0.1){
    # input is coordinates if inflection point coordinates are fed
    # input is distance if distances are fed
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix)
    # and average CV matrix without test elec and outputs a matrix of study bags
    # includes a penalization term to promote diversity that is a difference of convex functions
    # (l2 norm of diference in bags)
    # tol is the number of iterations before convergence until restart starting vector 
    # percent is the percentage of adjacent distances that the maximum promotion can achieve
    
    distances <- sort(as.vector(distances))
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    # bag matrix
    bag.mat <- matrix(NA, ncol = bag.size, nrow = length(distances)) # the actual bag
    bag.recipe <- matrix(NA, ncol = K - 1, nrow = length(distances)) # each element in vector is the number of studies from the corresponding study
    colnames(bag.recipe) <- paste0("P_", study.num)
    row.names(bag.recipe) <- paste0("Target",1:length(distances))
    #### Prepare Constrained Optimization ####
    # load package
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    #w <- Variable(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    B <- matrix(NA, nrow = 1, ncol = bag.size)
    objective <- Minimize( sum((elec.target - X %*% w)^2)  )
    #objective <- Minimize( sum((elec.target - X %*% w )^2))
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    B[1,] <- round(result$getValue(w)) # round since the 0s are not exactly 0
    tol <- 30
    
    # B[1,] <- sample(1:bag.size, bag.size, replace = TRUE) #toy example
    
    # target electrode
    
    for (i in 1:length(distances)){
        # integer program
        d <- distances[i]
        
        # function that stays constant for each iteration
        g.fn <- function(X, w, elec.target, d){
            sum((elec.target - X %*% w)^2)^2 + d^2
        }
        
        cnvg <- 1 # starting value
        w1 <- rep(1, K - 1) # starting value
        w0 <- w1 # start with w0 as  previous iteration value
        h.fn <- function(X, w, w1, elec.target, d){
            2 * d * sqrt( sum((X %*% w1 - elec.target)^2) ) + 
                t( (2 * d * t(X) %*% (X %*% w1 - elec.target) ) / (sqrt( sum((X %*% w1 - elec.target)^2))) ) %*% (w - w1)
        }
        counter <- 0 # number of iterations before convergence
        
        fn.itr <- sum_squares(elec.target - X %*% w)^2 + d^2 -  (2 * d * sqrt( sum((X %*% w1 - elec.target)^2) ) + 
                                                                     t( (2 * d * t(X) %*% (X %*% w1 - elec.target) ) / (sqrt( sum((X %*% w1 - elec.target)^2))) ) %*% (w - w1) )
        
        while (cnvg > 0.1){
            
            counter <- counter + 1
            w0 <- w1
            
            if(counter > tol){
                # if it is not convering try a new starting vector
                print("reset")
                counter <- 0
                w1 <- rand_vect(length(study.num), bag.size) # random vector of starting value
            }
            
            objective <- Minimize( fn.itr  )
            #objective <- Minimize( sum((elec.target - X %*% w )^2))
            problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
            result <- solve(problem)
            w1 <- round(result$getValue(w)) # round since the 0s are not exactly 0
            print(w1)
            
            # update iteration
            
            
            #lambda <- lambda.gen(distances, indx = i, percent = percent, bag.size = bag.size)
            fn.itr <- g.fn(X, w, elec.target, d) - h.fn(X, w, w1, elec.target, d) # 
            
            cnvg <- sum ( (w0 - w1)^2  ) # difference in l2 norm between previous iteration and current one
        }
        print("DCA complete")
        # after convergence, update B to include new solution
        B <- rbind(B, t(w1))
        
        # make bag
        bag <- c()
        for (x in 1:length(w1)){
            bag <- c(bag, rep(x, w1[x]))
        }
        #The "bag" is not the stud numbers but the element number index (i.e., it doesnt take out test elec)
        # so need to use list of elecs (study num) and index them by bag
        bag.mat[i,] <- study.num[bag] # add bag to matrix.
        bag.recipe[i, ] <- w1  # add bag recipe to
        
    } # double check placement of this bracket
    
    return( list( unique(bag.mat), unique(bag.recipe)  ) )
}



bagGen.base <- function(full.CVs, test.elec, bag.size = nrow(full.CVs) - 1){
    # input is coordinates if inflection point coordinates are fed
    # input is eta is the distance ot add to the minimum to be the circumfrence of the campfire
    # feed it coordinates of (r) target CVs (should be an r x 8 matrix) 
    # and average CV matrix without test elec and outputs a matrix of study bags
    # straps are number to go until stopping
    
    
    # import wiggle.fn
    source("Study Strap Functions.R")
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,-c(1,2)])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
   
    # minimum dist
    objective <- Minimize(sum((elec.target - X %*% w )^2))  # 1/K part is to 
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size))
    result <- solve(problem)
    w.solution <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # make bag
    bag <- c()
    for (x in 1:length(w.solution)){
        bag <- c(bag, rep(x, w.solution[x]))
    }
    
    bag <- study.num[bag] # remake with correct study indices
    return(bag)
}




AR.immitate <- function(full, full.CVs, test.elec, bag.size = nrow(full.CVs) - 1, paths = 5, 
                    converge.lim = 50000, prox.straps = 0, write_dist = FALSE, max.straps = 150){
    # provides distance vector based upon fitting an exponential curve to random sample of paths
    # full.CVs is matrix of average features (not sure what else it includes)
    # test.elec is the study number of the target 
    # paths is the number of random paths to sample
    # converge.lim is the number of random samples to take before ending the path
    # prox.straps - proximal straps -- the number of extra straps (in addition to the optimum) very close to the target that we include PER PATH in the model fitting process
    # prox.straps also determines the total number of straps for each electrode-- is just the average straps/path + prox.straps is the total straps to be fit
    # write_dist indicator of whether to write a csv of the distances
    
    setwd("/n/home12/gloewinger")
    source("Study Strap Functions.R")
    total.proximal <- prox.straps * paths # total number of proximal straps to optimize to include in model fitting step
    suppressWarnings(library(CVXR, warn.conflicts=FALSE))
    
    K <- nrow(full.CVs)# number of studies including test elec
    study.num <- seq(1, K)[-test.elec] # the elecs included
    full.CVs <- full.CVs[,-c(1,2)] # remove all but covariates
    
    # make objective matrix A
    A <- matrix(NA, nrow = 8, ncol = K)
    colnames(A) <- paste0("Elec_", 1:K)
    row.names(A) <- paste0("Coordinate_", 1:8)
    for (i in 1:K){
        A[,i] <- wiggle.fn(full.CVs[i,])
    }
    
    elec.target <- A[,test.elec] # the "target" in the objective
    
    A <- A[,-test.elec] # remove for optimization below
    
    
    target.elec <- full.CVs[test.elec,]
    full.CVs[test.elec,] <- NA
    CVs <- full.CVs[-test.elec,] # remove test elec
    elecs <- study.num
    
    ###############################################
    # Optimization Step prior to curve fitting
    ###############################################
    #### Prepare Constrained Optimization ####
    print("opt minimum")
    w <- Int(K - 1)
    X <- A * (1/bag.size) # Objective Matrix X reformulated so we can use integer constraints
    
    # optimize
    objective <- Minimize( sum((elec.target - X %*% w)^2) )
    problem <- Problem(objective, constraints = list(w >= 0, sum(w) == bag.size)) #sum(w) == bag.size)) , sum(w^2) <= 100
    result <- solve(problem)
    w.opt <- round(result$getValue(w)) # round since the 0s are not exactly 0
    
    # make bag
    bag <- c()
    for (x in 1:length(w.opt)){
        bag <- c(bag, rep(x, w.opt[x]))
    }
    bag <- study.num[bag]
    
    bagRows.list <- list(NA) # each element is a vector of row indices of full that are accepted straps in AR
    # accept reject step
    counter <- 0 # number of samples before acceptance
    strap <- 1 # strap counter within path
    total.straps <- 1 # the total number of straps across all paths
    thresh <-  5e+14 # restart value: arbitraily large starting threshold
    threshold <- thresh # start at restart value
    distance.matrix <- matrix(NA, ncol = 5, nrow = 500) # arbitrarily large number of rows
    colnames(distance.matrix) <- c("Strap", "Distance", "Count", "Path", "Time")
    path.straps <- vector(length = paths) # total number of straps in each path
    for( path in 1:paths){
        
        timeStart <- Sys.time() 
        print(paste0("path ", path))
        # iterate through number of paths
        while (counter < converge.lim){ #while loop for restrapping: breaks for either acceptance (sim) or accpetance counter limit (converge.lim)
            
            # Choose electrodes to be in bag:
            
            bag <- sample(elecs, bag.size, replace = TRUE) # sample bag
            feats <- wiggle.fn( colMeans( full.CVs[bag, ]) ) # calculate inflection points
            current.dist <- sum(  (feats - elec.target)^2  ) # squared Euclidean Distance between inflection points
            strap.table <- as.data.frame(table(bag)) #sample with replacement and put in data table format
            
            if(current.dist < threshold ){
                #if accepted: generate study strap row indices
                
                ########################################
                # remake pseudo study
                ########################################
                indx <- c() # vector of indices corresponding to rows being sub-sampled
                for(i in 1:nrow(strap.table)){
                    
                    sub.elec <- as.numeric(as.character(strap.table[i,1])) # current electrode being sampled from
                    elec.indx <- which(full$Electrode == sub.elec) # rows corresponding to electrode sampled in current bag
                    num.obs <- round( (length(elec.indx) / bag.size) * as.numeric(as.character(strap.table[i,2])) )
                    elec.indx <- elec.indx[sample(1:length(elec.indx), num.obs, replace = FALSE)] #sub-sample rows
                    indx <- c(indx, elec.indx ) # sample as many rows as indicated in elec.obs and add to indices vector
                   }
                
                bagRows.list[[strap]] <- indx
                rm(elec.indx)
                rm(indx)
                ###########################################
                
                counter <- 0 # reset counter
                strap <- strap + 1 # increase straps of current path
                total.straps <- total.straps + 1 # increase total straps
                threshold <- current.dist
            }else{
                # add counter
                counter <- counter + 1
            }
            
        }

        timeEnd <- Sys.time()      
        path.time <- as.numeric( difftime(timeEnd, timeStart, units='mins') )
        print(paste0("path ", path, "_elec_", test.elec, "_time_", path.time))
        # add optimum - counter = 0 so we know which is from optimization
        
        path.straps[path] <- strap # save the total number of straps in this path
        strap <- 1 # reset strap counter for this bag
        threshold <- thresh # restart threshold
        counter <- 1
        total.straps <- total.straps + 1
    }
    
     return( bagRows.list ) 
    
}




optFeatBag <- function(data,
                       mtry = NA,
                       R = 150,
                       corr = TRUE,
                       selfCorr = FALSE,
                       quad = FALSE,
                       lambda = 0.5,
                       initial = "all",
                       solver = "MOSEK",
                       thresh = 0,
                       rhoVec = NULL,
                       CDiters = 1){
    
    
    # mtry is number of feature bags to make
    
    # corr is an indicator of whether to include the
    #       feature correlation matrix in the objective
    
    # selfCorr is whether to include the S matrix in
    #       the objective with itself z_j^T S z_j
    # quad is whether to include the quadratic form for the 
    # self correlation
    
    # initial is how to initialize: has 3 values
    #       # "all" is for a model with all
    #       #  "top" is mtry most strongly correlated
    #       #  "random" randomly chooses mtry
    
    # data is expected to be in form Y   X_1    X_2 .... X_p
    
    # thresh is the threshold of correlations to be cut off--keep at 0 usually
    
    library(CVXR)
    p <- ncol(data) - 1 # Subtract 1 for y 
    C <- matrix(0, nrow = R, ncol = p)
    
    if(is.na(mtry)){
        # if not given then give heuristic of round p / 3 as number of features 
        mtry <- round(p / 3)
    }
    
    # pairwise correlations between y and each covariate--use if user does not provide
    if(is.null(rhoVec)){
        rhoVec <- abs( apply(data[,-1], 2, function(x) cor(data$Y, x)) )
    }
    
    #rhoVec <- scale(rhoVec, center = FALSE)
    
    # pairwise correlations between each covariate
    if(corr){
        
        # if not fast enough can change this to X^T X after scaling 
        # the data and then divide by the value on the diagonal. then equivalent
        
        S <- abs( cor(data[,-1], data[,-1]) ) 
        
        #S <- S * upper.tri(S) # upper triangle so we do not duplicate values
        
        Svec <- c( as.matrix( S * upper.tri(S, diag = FALSE) ) )
        Svec <- Svec[Svec != 0]
        S <- S * I(S >= thresh)
        Svec <- as.vector( Svec * I(Svec >= thresh) )
        
    }else{
        S <- diag(p)
    }
    
    
    # initialize first row with random features
    if(initial == "all"){
        
        # all features in the first feature bag
        C[1, ] <- 1 
        
    }else if(initial == "random"){
        # random one in feature bag
        indx <- sample.int(p, mtry) # randomly choose indices of features
        C[1, indx] <- 1 
        
    }else if(initial == "top"){
        # first feature bag is the top mtry features most correlated with outcome
        
        indx <- sort(rhoVec, index.return=TRUE, decreasing=TRUE) 
        indx <- indx$ix[1:mtry] # top most correlated features
        C[1, indx] <- 1 
    }
    
    
    if(selfCorr){
        if(quad){
            # if include quadratic form for z^T S z
            
            for(r in 2:R){
                print("quad")
                
                message(paste("Feature Bag", r))
                
                z <- CVXR::Variable(p, boolean = TRUE)
                
                obj <- CVXR::Maximize( 
                    lambda * sum_entries( t(z) %*% rhoVec ) -
                        (1 - lambda) / r * (sum_entries( C %*% S %*% z ) + quad_form(z, S) )
                ) # loss
                
                
                problem <- CVXR::Problem(obj, 
                                   constraints = 
                                       list(
                                           sum_entries(z) == mtry
                                       )
                ) 
                result <- CVXR::solve(problem, solver = solver, warm_start = TRUE) 
                
                C[r,] <- round(result$getValue(z), 6)
            }
            
        }else{
            # if not quadratic form then use vector representation
            print("self corr vector")
            #######################
            # Contrast Matrix
            #######################
            
            CMat <- matrix(0, ncol = p, nrow = choose(p, 2) ) # contrast matrix
            FMat <- matrix(0, ncol = p, nrow = choose(p, 2) ) # contrast matrix but with only first 1 in each row
            SMat <- matrix(0, ncol = p, nrow = choose(p, 2) ) # contrast matrix but with only second 1 in each row
            
            cr <- er <- 0 # keep track of current row
            for(cov in 1:(p-1) ){
                
                cr <- er + 1 # lower bound of rows
                er <- cr + p - cov - 1 # upper bound of rows
                
                CMat[cr:er, cov] <- 1
                FMat[cr:er, cov] <- 1
                
                ident <- diag(1, p - cov) # identity matrix of corresponding size
                
                CMat[cr:er, (1 + cov):p] <- ident
                SMat[cr:er, (1 + cov):p] <- ident
            }
            
            # remove duplicate rows if there are any (shouldnt be)
            indx <- apply(CMat, 1, function(x) length(unique(x)) == 1)
            CMat <- CMat[!indx,]
            
            
            for(r in 2:R){
                
                message(paste("Feature Bag", r))
                
                z <- CVXR::Variable(p, boolean = TRUE)
                alphaVec <- CVXR::Variable( choose(p, 2) )
                
                obj <- CVXR::Maximize( 
                    lambda * ( t(z) %*% rhoVec ) -
                        (1 - lambda) / r * (sum_entries( C %*% S %*% z ) + t(alphaVec) %*% Svec  )   
                ) # loss
                
                problem <- CVXR::Problem(obj, 
                                   constraints = 
                                       list(
                                           alphaVec <= FMat %*% ( z ), #, # ensures each individual alpha_i <= z_i 
                                           alphaVec <= SMat %*% ( z ), # # ensures each individual alpha_i <= z_j for all i, j \in R choose 2 
                                           alphaVec >= CMat %*% ( z ) - 1,
                                           alphaVec >= 0, # elementwise
                                           sum(z) == mtry # 
                                       ) ) 
                result <- CVXR::solve(problem, solver = solver, warm_start = TRUE) 
                
                C[r,] <- round(result$getValue(z), 6)
            }
            
        }
        
        
    }else{
        # if no self correlation
        print("no self corr")
        # if include quadratic form for z^T S z
        
        for(cD in 1:CDiters){
            # number of coordinate descent iterations
            for(r in 2:R){
                
                message(paste("Feature Bag", r))
                
                z <- CVXR::Variable(p, boolean = TRUE)
                
                obj <- CVXR::Maximize( 
                    lambda * sum_entries( t(z) %*% rhoVec ) -
                        (1 - lambda) / r * (sum_entries( C %*% S %*% z ) )
                ) # loss
                
                
                problem <- CVXR::Problem(obj, 
                                         constraints = 
                                             list(
                                                 sum_entries(z) == mtry
                                             )
                ) 
                
                result <- CVXR::solve(problem, solver = solver, warm_start = TRUE) 
                
                C[r,] <- round(result$getValue(z), 6)
            }
            
        }
            
        
    }
    
    
    
    # convert matrix into indices of selected instead of 0 and 1 indicators
    C_opt <- apply( C, 1, function(x) which(x == 1) )
    #C_opt <- do.call(rbind, C_opt) # doesnt work because of dimensions if first has all of them
    
    return(list( indic = C, indx = t(C_opt) ) )
    
    
    
}



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
            # indxList[[study]] <- which(data$Study != trainStudies[study]) # indices of studies to train on
            # HOOList[[study]] <- which(data$Study == trainStudies[study]) # indices of study to hold out
            # changed 9/11/20 to make indices match
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
            
            # indxList[[1]] <- which(data$Study == trainStudies[study]) # indices of studies to train on
            # HOOList[[1]] <- which(data$Study != trainStudies[study]) # indices of study to hold out
            
            # changed 9/11/20 to make indices match
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
            # indxList <- which(data$Study == trainStudies[study]) # indices of studies to train on
            # HOOList <- which(data$Study != trainStudies[study]) # indices of study to hold out
            # changed 9/11/20 to make indices match
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
            # indxList <- which(data$Study == trainStudies[study]) # indices of studies to train on
            # HOOList <- which(data$Study != trainStudies[study]) # indices of study to hold out
            # changed 9/11/20 to make indices match
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
                # oneVec <- rep(1, length(HOOList))
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
            
            # HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            # indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on 
            # changed 9/11/20 t0 make indices match
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
            
            # HOOList <- which(data$Study == trainStudies[study]) # indices of study to hold out
            # indxList <- which(data$Study != trainStudies[study]) # indices of studies to train on 
            
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                        # oneVec <- rep(1, length(HOOList))
                        # preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                        # rmseMat[f, tune]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
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
            #wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                        # oneVec <- rep(1, length(HOOList))
                        # preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                        # rmseMat[f, tune]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
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
            # wgt <- weights[indxList, indxList] * length(indxList) / sum(weights[indxList, indxList])
            
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
                
                # indx <- which(data$Study == trainStudies[study]) # indices of study to hold out
                # changed on 9/11/20 to make indices match
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
                        # oneVec <- rep(1, length(HOOList))
                        # preds <- as.vector(  as.matrix( cbind(1, data[HOOList, -c(1,2)]) ) %*% betaHat   )
                        # rmseMat[f, tune]  <- sqrt(mean( (preds - data$Y[HOOList] )^2  )) # rmse
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


