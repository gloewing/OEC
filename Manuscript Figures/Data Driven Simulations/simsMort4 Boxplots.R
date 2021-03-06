# generalist specialist, zero out
######################################################
# compare different K and betaVar -- no regularization 
######################################################
setwd("~/Desktop/Research/simsMort4")
# setwd("~/Desktop/OEC/Manuscript Figures/Data Driven Simulations/simsMort4")
library(latex2exp)
library(tidyverse)
library(ggplot2)
thresh <- 0.5 
bVar <- c(0.25, 0.5, 1, 2, 4) # beta variances
tune <- c("zero","cvCF") 
minTime <- 100
Kvec <- c(2, 5, 8, 25, 30) 

year <- c(2003) 
ls_beta <- vector(length = length(bVar), "list")
ls_K <- vector(length = length(Kvec), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(tn in tune){
    for(yrs in year){
        for(b in 1:length(bVar)){
            for(k in 1:length(Kvec)){
                cnt <- cnt + 1 # counter
                K <- Kvec[k]
                x <- bVar[b]
                
                flNm <- paste0("oecMort_0.5.etTn_TRUE_cv_stTn_cv_oecTn_", tn, "_stTnIn_TRUE_stCV_cv_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_sclX_TRUE_glm_FALSE_sampSz_NA_smpSzW_1_tstTrn_52_minTrn_100_bVar_", x, "_reTn_FALSE_eta_TRUE_yr_",yrs,"_K_",K)
                
                if(file.exists(flNm)){
                    # check to see if file exists
                    d <- read.csv(flNm) 
                    
                    rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                          Specialist = d$oec_country2 / d$stacking_country,
                                          ZeroOut = d$oec_country0 / d$stacking_country_zeroOut
                    )
                    
                    ls_beta[[cnt]] <- cbind( 
                        gather(rmseMat), 
                        K, 
                        x,
                        yrs,
                        tn
                    )
                    
                }
            }
            
        }
    }
}
#############
# process data
#############
dat <- do.call(rbind, ls_beta)
dat$K <- as.factor(dat$K)
dat$x <- as.factor(dat$x)
dat$yrs <- as.factor(dat$yrs)
dat$tn <- as.factor(dat$tn)

#############
# figures
#############
meths <- c("Specialist", "ZeroOut")
tune <- c("zero", "cvCF") 

for(t in tune){
    for(m in meths){
        
        plt1 <- dat %>% 
            as_tibble() %>% 
            dplyr::filter(key == m)  %>% 
            dplyr::filter(tn == t) %>%
            ggplot(aes(y = value, x = x, fill = K, color = K)) + 
            geom_hline(yintercept = 1,
                       linetype   = 2,
                       color      = "darkgray") +
            geom_boxplot(size         = 0.20,
                         outlier.size = 0.50) +
            geom_boxplot(size          = 0.20,
                         color         = "black",
                         outlier.color = NA) +
            coord_cartesian(ylim = c(0.5, 1.5)) +
            scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Reds")) +
            scale_color_manual(values = RColorBrewer::brewer.pal(5, "Reds")) +
            labs(x = TeX('$\\mathbf{\\sigma^2_{\\theta}}}$'),
                 y = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{MSS}}$'),
                 color = "Studies (K)",
                 fill  = "Studies (K)") +
            theme_bw() + 
            theme(text = element_text(face = "bold"))

        # save figures
        setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Data-Driven Simulations/simsMort4 Figures")
        ggsave(paste0("mortSims_", m, "_", t, ".pdf"),
               plot   = plt1,
               dpi    = 300,
               width  = 6, 
               height = 4)
    }
}
