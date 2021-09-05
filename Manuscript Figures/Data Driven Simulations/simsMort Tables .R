######################################################
# compare different K and betaVar -- no regularization
# no W specialist training
######################################################
# zero
setwd("~/Desktop/Research/simsMort4")
library(latex2exp)
library(ggplot2)
library(kableExtra)

bVar <- c(0.25, 0.5, 1, 2, 4) # beta variances
tune <- c("zero")
clust <- 0 # no clusters
Kvec <- c(2, 5, 8, 25, 30) 

rmseMat2 <- rmseMat1 <- rmseMat <- matrix(nc = length(bVar), nr = length(Kvec))
colnames(rmseMat2) <- colnames(rmseMat1) <- colnames(rmseMat) <- bVar
rownames(rmseMat2) <- rownames(rmseMat1) <- rownames(rmseMat) <- Kvec

seMat2 <- seMat <- resMat <- matrix(nc = 6, nr = length(bVar) *  length(Kvec))
colnames(resMat) <- c("K", "bVar", "OEC", "Spec", "OEC_Zero", "Zero")
  
compMat <- array(data = NA, dim = c(  length(tune), length(bVar), 3  ),
                 dimnames = list(        )
)# first dimension is 3 because of 3 oec versions (generalist, zero out, standard)

# r.0 <- r.sse <- r.cv <- matrix(NA, nr = 3, nc = length(bVar))

compMat[1,,] # first dimension corresponds to tuning style


cnt <- 0
compMat2 <- array(data = NA, dim = c(  length(Kvec), length(bVar), 3  ) )
                  
for(t in 1:length(Kvec)){
    
    tn <- tune[1] # "zero"
    k <- Kvec[t]
    
    for(j in 1:length(bVar)){
        
        cnt <- cnt + 1
        setwd("~/Desktop/Research/simsMort4")
        x <- bVar[j]
        flNm <- paste0("oecMort_0.5.etTn_TRUE_cv_stTn_cv_oecTn_", tn, "_stTnIn_TRUE_stCV_cv_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_sclX_TRUE_glm_FALSE_sampSz_NA_smpSzW_1_tstTrn_52_minTrn_100_bVar_", x, "_reTn_FALSE_eta_TRUE_yr_2003_K_", k)
        
        if(file.exists(flNm)){
            # check to see if file exists
            d <- read.csv(flNm) 
            
            # number of simulation iterations for monte carlo error estimation
            numItrs <- nrow(d)
            
            rmseMat[t,j] <- mean( d$oec_test2 / d$stack2 )
            rmseMat1[t,j] <- mean( d$oec_country2 / d$stacking_country )
            rmseMat2[t,j] <- mean( d$oec_country0 / d$stacking_country_zeroOut )
            
            compMat2[t, j, ] <- colMeans(cbind(d$oec_test2 / d$stack2,
                                              d$oec_country2 / d$stacking_country,
                                              d$oec_country0 / d$stacking_country_zeroOut
            )
            )
            
            resMat[cnt, ] <- c(k, x, 
                               colMeans(
                                 cbind( d$oec_country2 / d$country,
                                        d$stacking_country / d$country,
                                        d$oec_country0 / d$country,
                                        d$stacking_country_zeroOut / d$country
                                      )
                                        )
                            )
            
            # monte carlo error
            seMat[cnt, ] <- c(k, x, 
                               apply(
                                 cbind( d$oec_country2 / d$country,
                                        d$stacking_country / d$country,
                                        d$oec_country0 / d$country,
                                        d$stacking_country_zeroOut / d$country
                                 ), 2, sd
                               ) / numItrs
            )
            
            # monte carlo error 2
            seMat2[cnt, ] <- c(k, x, 
                              apply(
                                cbind( d$oec_country2 / d$stacking_country,
                                       d$oec_country0 / d$stacking_country_zeroOut
                                ), 2, sd
                              ) / numItrs, NA, NA
                              
            )
        }
        
    }
    
}



##################
# tables
##################

library(kableExtra)

# specialist no regulariztion
kable( round( (rmseMat1), 2   ),  
      format = "latex", booktabs = T) %>% kable_styling(position = "center")

# zero-out no regularization
kable( round( (rmseMat2), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# combined specialist and zero out into one table
combinedMat <- cbind(rmseMat1, rmseMat2)
kable( round( cbind(rmseMat1, rmseMat2), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# monte carlo error 
apply(seMat[,3:6], 2, max)
apply(seMat2[,3:4], 2, max)

# table vs. country-specific model - no regularization
                    
kable( round( (resMat), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

rmseMat_reduced <- resMat %>% 
  as_tibble %>% 
  filter(K %in% c(2, 8, 30) ) %>% 
  filter(bVar %in% c(0.25, 1, 4) )

kable( round( (rmseMat_reduced), 3   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")


######################################################
# ridge term
######################################################
# cvCF
setwd("~/Desktop/Research/simsMort4")
library(latex2exp)
library(ggplot2)
library(kableExtra)

bVar <- c(0.25, 0.5, 1, 2, 4) # beta variances
tune <- c("cvCF")
Kvec <- c(2, 5, 8, 25, 30) 

rmseMat2 <- rmseMat1 <- rmseMat <- matrix(nc = length(bVar), nr = length(Kvec))
colnames(rmseMat2) <- colnames(rmseMat1) <- colnames(rmseMat) <- bVar
rownames(rmseMat2) <- rownames(rmseMat1) <- rownames(rmseMat) <- Kvec

seMat2 <- seMat <- resMat <- matrix(nc = 6, nr = length(bVar) *  length(Kvec))
colnames(resMat) <- c("K", "bVar", "OEC", "Spec", "OEC_Zero", "Zero")

compMat <- array(data = NA, dim = c(  length(tune), length(bVar), 3  ),
                 dimnames = list(        )
)# first dimension is 3 because of 3 oec versions (generalist, zero out, standard)

# r.0 <- r.sse <- r.cv <- matrix(NA, nr = 3, nc = length(bVar))

compMat[1,,] # first dimension corresponds to tuning style

cnt <- 0
compMat2 <- array(data = NA, dim = c(  length(Kvec), length(bVar), 3  ) )

for(t in 1:length(Kvec)){
  
  tn <- tune[1] # "zero"
  k <- Kvec[t]
  
  for(j in 1:length(bVar)){
    
    cnt <- cnt + 1
    setwd("~/Desktop/Research/simsMort4")
    x <- bVar[j]
    flNm <- paste0("oecMort_0.5.etTn_TRUE_cv_stTn_cv_oecTn_", tn, "_stTnIn_TRUE_stCV_cv_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_sclX_TRUE_glm_FALSE_sampSz_NA_smpSzW_1_tstTrn_52_minTrn_100_bVar_", x, "_reTn_FALSE_eta_TRUE_yr_2003_K_", k)
    
    if(file.exists(flNm)){
      # check to see if file exists
      d <- read.csv(flNm) 
      
      # number of simulation iterations for monte carlo error estimation
      numItrs <- nrow(d)
      
      rmseMat[t,j] <- mean( d$oec_test2 / d$stack2 )
      rmseMat1[t,j] <- mean( d$oec_country2 / d$stacking_country )
      rmseMat2[t,j] <- mean( d$oec_country0 / d$stacking_country_zeroOut )
      
      compMat2[t, j, ] <- colMeans(cbind(d$oec_test2 / d$stack2,
                                         d$oec_country2 / d$stacking_country,
                                         d$oec_country0 / d$stacking_country_zeroOut
      )
      )
      
      resMat[cnt, ] <- c(k, x, 
                         colMeans(
                           cbind( d$oec_country2 / d$country,
                                  d$stacking_country / d$country,
                                  d$oec_country0 / d$country,
                                  d$stacking_country_zeroOut / d$country
                           )
                         )
      )
      
      # monte carlo error
      seMat[cnt, ] <- c(k, x, 
                        apply(
                          cbind( d$oec_country2 / d$country,
                                 d$stacking_country / d$country,
                                 d$oec_country0 / d$country,
                                 d$stacking_country_zeroOut / d$country
                          ), 2, sd
                        ) / numItrs
      )
      
      # monte carlo error 2
      seMat2[cnt, ] <- c(k, x, 
                         apply(
                           cbind( d$oec_country2 / d$stacking_country,
                                  d$oec_country0 / d$stacking_country_zeroOut
                           ), 2, sd
                         ) / numItrs, NA, NA
                         
      )
    }
    
  }
  
}

##################
# tables with regularization
##################

# specialist with regulariztion
kable( round( (rmseMat1), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# zero-out with regularization
kable( round( (rmseMat2), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# combined specialist and zero out into one table
combinedMat <- cbind(rmseMat1, rmseMat2)
kable( round( cbind(rmseMat1, rmseMat2), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

# monte carlo error 
apply(seMat[,3:6], 2, max)
apply(seMat2[,3:4], 2, max)

# table vs. country-specific model - with regularization

kable( round( (resMat), 2   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")

rmseMat_reduced <- resMat %>% 
  as_tibble %>% 
  filter(K %in% c(2, 8, 30) ) %>% 
  filter(bVar %in% c(0.25, 1, 4) )

kable( round( (rmseMat_reduced), 3   ),  
       format = "latex", booktabs = T) %>% kable_styling(position = "center")


