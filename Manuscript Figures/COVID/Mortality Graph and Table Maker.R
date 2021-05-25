#######################
# ggplot2 Specialist
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)
nhVec <- c(FALSE)
minYrs <- c(100) # actually in months
mnYr <- 100
    for(nH in nhVec){
    # cycle through minimum years and make a different figure for each
    for(mnYr in minYrs){
        
    if(mnYr == 150){
        yrs <- 2003:2019 # 1997:2019
    }else{
        yrs <- 2003:2019 # 1997:2019
    }
    tnVec <- c( "zero", "cvCF") # , "cvCF"
    proVec <- c(4) # linear term for time
    
    for(tn in tnVec){
        
        for(xPro in proVec){
            
            setwd("~/Desktop/Research/mort17")
            fileNm <- paste0("mrtNew_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nH",nH,"_mnSt_-01-01_mnTr_", mnYr, "_fld_5_smpSzWt_1")
            if(file.exists(fileNm)){
                
           
            a <- read.csv( fileNm )
            
            # remove southern hemisphere
            a <- a %>% tibble %>% dplyr::filter( !country.1 %in% c("New Zealand", "Chile", "Australia DCD") )

            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            par(mfrow = c(1,1))
            
            mat <- mat2 <- mat3 <- mat4 <- mat5 <- mat6 <- matrix(NA, nrow(a), ncol = length(yrs))
            mat9 <- mat10 <- mat7 <- mat8 <- mat
            for(j in 1:length(yrs)){
                mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
                mat2[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
                
                mat3[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country_noLinear[a$testYear == yrs[j] ]
                mat4[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$country_noLinear[a$testYear == yrs[j] ]
                
                # choose with samllest objective among different methods
                objSpec <- a %>% tibble() %>% 
                            dplyr::filter(testYear == yrs[j]) %>% 
                            dplyr::select(oec_specAnneal_leftObj, oec_spec_leftObj, oec_specAnneal_rightObj, oec_spec_rightObj )
                
                indx <- apply(objSpec, 1, which.min)
                
                rmseSpec <- a %>% tibble() %>% 
                        dplyr::filter(testYear == yrs[j]) %>% 
                        dplyr::select(oec_specAnneal_left, oec_spec_left, oec_specAnneal_right, oec_spec_right)
                
                # select ones with smallest objective (based on index)
                indx <- cbind(1:nrow(rmseSpec), indx)
        
                mat5[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmseSpec)[indx] / a$country_noLinear[a$testYear == yrs[j] ]
                mat7[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmseSpec)[indx] / a$oec_country2[a$testYear == yrs[j] ]
                mat9[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmseSpec)[indx] / a$stacking_country[a$testYear == yrs[j] ]
                
                # zero out
                
                # choose with samllest objective among different methods
                obj0 <- a %>% tibble() %>% 
                    dplyr::filter(testYear == yrs[j]) %>% 
                    dplyr::select(oec_0Anneal_leftObj, oec_0_leftObj, oec_0Anneal_rightObj, oec_0_rightObj )
                
                indx <- apply(obj0, 1, which.min)
                
                rmse0 <- a %>% tibble() %>% 
                    dplyr::filter(testYear == yrs[j]) %>% 
                    dplyr::select(oec_0Anneal_left, oec_0_left, oec_0Anneal_right, oec_0_right)
                
                # select ones with smallest objective (based on index)
                indx <- cbind(1:nrow(rmse0), indx)
                
                mat6[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmse0)[indx] / a$country_noLinear[a$testYear == yrs[j] ]
                mat8[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmse0)[indx] / a$oec_country0[a$testYear == yrs[j] ]
                mat10[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- as.matrix(rmse0)[indx] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
                
            }
          
            mat <- as.data.frame(mat)
            colnames(mat) <- yrs
            m <-cbind( gather(mat), "Specialist")
            
            mat2 <- as.data.frame(mat2)
            colnames(mat2) <- yrs
            m1 <- cbind(gather(mat2), "Zero Out")
            
            mat3 <- as.data.frame(mat3)
            colnames(mat3) <- yrs
            m3 <- cbind(gather(mat3), "Spec/noLin")
            
            mat4 <- as.data.frame(mat4)
            colnames(mat4) <- yrs
            m4 <- cbind(gather(mat4), "ZeroOut/noLin")
            
            
            mat5 <- as.data.frame(mat5)
            colnames(mat5) <- yrs
            m5 <- cbind(gather(mat5), "Spec_bestObj/noLin")
            
            mat6 <- as.data.frame(mat6)
            colnames(mat6) <- yrs
            m6 <- cbind(gather(mat6), "ZeroOut_bestObj/noLin")
            
            mat7 <- as.data.frame(mat7)
            colnames(mat7) <- yrs
            m7 <- cbind(gather(mat7), "Spec_bestObj/Spec")
            
            mat8 <- as.data.frame(mat8)
            colnames(mat8) <- yrs
            m8 <- cbind(gather(mat8), "ZeroOut_bestObj/OEC0")
            
            mat9 <- as.data.frame(mat9)
            colnames(mat9) <- yrs
            m9 <- cbind(gather(mat9), "Spec_bestObj/Stack")
            
            mat10 <- as.data.frame(mat10)
            colnames(mat10) <- yrs
            m10 <- cbind(gather(mat10), "ZeroOut_bestObj/Stack0")
            
            colnames(m)[3] <- colnames(m1)[3] <- colnames(m3)[3] <- colnames(m4)[3] <- colnames(m5)[3] <- colnames(m6)[3] <- colnames(m7)[3] <- colnames(m8)[3] <- "Method"
            colnames(m9)[3] <- colnames(m10)[3]<- "Method"
            m <- rbind(m, m1, m3, m4, m5, m6, m7, m8, m9, m10)
            
            colnames(m)[3] <- "Method"
            
            t <- ifelse(tn == "zero", "OLS", "Ridge")
            
            m <- cbind(m, factor(t) )
            colnames(m) <- c("Year", "RMSE", "Method", "Regularization")
            
            if(tn == "zero"){
                df <- m
            }else{
                # bind them together
                df <- rbind(df, m)
            }
            
            }
            
        }
    }
            df <- df[!is.na(df$RMSE),]
            df$Year <- as.factor(df$Year)
            
            plt = df %>% tibble %>%  filter(df$Method %in% c("Specialist", "Zero Out")) %>%
                # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
                # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
                # dplyr::filter(cl == 3) %>%
                # dplyr::filter(n == 300) %>%
                #dplyr::filter(Method == "OLS") %>%
                ggplot(aes( y = RMSE, x = Year, fill = Regularization )) + # 
               facet_wrap( ~ Method, nrow = 1) +
                geom_boxplot(
                    lwd = 1.5, 
                    fatten = 0.5, 
                    alpha = 0.5 
                ) + 
                geom_hline(yintercept=1, 
                           linetype="dashed", 
                           color = "black", 
                           size = rel(0.5),
                           alpha = 0.7) + #
                #ylim(0, 2) +
                ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
                xlab(TeX('Year')) + 
                # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
                scale_fill_manual(values = c("red", "blue")) + #"red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
                theme_classic(base_size = 12) +
                ylim(0.4, 1.1) + 
                theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                       axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                       axis.title = element_text(face="bold", color="black", size=rel(1.75)),
                       legend.key.size = unit(2, "line"), # added in to increase size
                       legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
                       legend.title = element_text(face="bold", color="black", size = rel(2)),
                       strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
                       #legend.position ="none"
                ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months"))  +  
                scale_x_discrete(breaks = seq(2003,2019)[c(1, 5,10,15, 17)]) 

            
            
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_", mnYr, ".png"),
                    plot = plt,
                    width = 22,
                    height = 7 
            )
            
            
            plt = df %>% tibble %>%  filter(df$Method %in% c("Spec_bestObj/Stack", "ZeroOut_bestObj/Stack0")) %>%
                # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
                # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
                # dplyr::filter(cl == 3) %>%
                # dplyr::filter(n == 300) %>%
                #dplyr::filter(Method == "OLS") %>%
                ggplot(aes( y = RMSE, x = Year, fill = Regularization )) + # 
                facet_wrap( ~ Method, nrow = 1) +
                geom_boxplot(
                    lwd = 1.5, 
                    fatten = 0.5, 
                    alpha = 0.5 
                ) + 
                geom_hline(yintercept=1, 
                           linetype="dashed", 
                           color = "black", 
                           size = rel(0.5),
                           alpha = 0.7) + #
                #ylim(0, 2) +
                ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
                xlab(TeX('Year')) + 
                # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
                scale_fill_manual(values = c("red", "blue")) + #"red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
                theme_classic(base_size = 12) +
                ylim(0.4, 1.1) + 
                theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                       axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                       axis.title = element_text(face="bold", color="black", size=rel(1.75)),
                       legend.key.size = unit(2, "line"), # added in to increase size
                       legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
                       legend.title = element_text(face="bold", color="black", size = rel(2)),
                       strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
                       #legend.position ="none"
                ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months")) +  
                scale_x_discrete(breaks = seq(2003,2019)[c(1, 5,10,15, 17)]) 
            
            
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_NewBestObj_", mnYr, ".png"),
                    plot = plt,
                    width = 22,
                    height = 7 
            )
            
            
            plt = df %>% tibble %>% filter(df$Method %in% c("Spec/noLin", "ZeroOut/noLin")) %>%
                ggplot(aes( y = RMSE, x = Year, fill = Regularization )) + # 
                facet_wrap( ~ Method, nrow = 1) +
                geom_boxplot(
                    lwd = 1.5, 
                    fatten = 0.5, 
                    alpha = 0.5 
                ) + 
                geom_hline(yintercept=1, 
                           linetype="dashed", 
                           color = "black", 
                           size = rel(0.5),
                           alpha = 0.7) + #
                #ylim(0, 2) +
                ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
                xlab(TeX('Year')) + 
                # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
                scale_fill_manual(values = c("red", "blue")) + #"red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
                theme_classic(base_size = 12) +
                ylim(0.4, 1.1) + 
                theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                       axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                       axis.title = element_text(face="bold", color="black", size=rel(1.75)),
                       legend.key.size = unit(2, "line"), # added in to increase size
                       legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
                       legend.title = element_text(face="bold", color="black", size = rel(2)),
                       strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
                       #legend.position ="none"
                ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months")) +  
                scale_x_discrete(breaks = seq(2003,2019)[c(1, 5,10,15, 17)]) 
            
            
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_NoLinear_", mnYr, ".png"),
                    plot = plt,
                    width = 22,
                    height = 7 
            )
            
            
            
            plt = df %>% tibble %>% filter(df$Method %in% c("Spec_bestObj/noLin", "ZeroOut_bestObj/noLin")) %>%
                ggplot(aes( y = RMSE, x = Year, fill = Regularization )) + # 
                facet_wrap( ~ Method, nrow = 1) +
                geom_boxplot(
                    lwd = 1.5, 
                    fatten = 0.5, 
                    alpha = 0.5 
                ) + 
                geom_hline(yintercept=1, 
                           linetype="dashed", 
                           color = "black", 
                           size = rel(0.5),
                           alpha = 0.7) + #
                #ylim(0, 2) +
                ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
                xlab(TeX('Year')) + 
                # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
                scale_fill_manual(values = c("red", "blue")) + #"red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
                theme_classic(base_size = 12) +
                ylim(0.4, 1.1) + 
                theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                       axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                       axis.title = element_text(face="bold", color="black", size=rel(1.75)),
                       legend.key.size = unit(2, "line"), # added in to increase size
                       legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
                       legend.title = element_text(face="bold", color="black", size = rel(2)),
                       strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
                       #legend.position ="none"
                ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months")) +  
                scale_x_discrete(breaks = seq(2003,2019)[c(1, 5,10,15, 17)]) 
            
            
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_NoLinearBest_", mnYr, ".png"),
                    plot = plt,
                    width = 22,
                    height = 7 
            )
            
            
            
            plt = df %>% tibble %>% filter(df$Method %in% c("Spec_bestObj/Spec", "ZeroOut_bestObj/OEC0")) %>%
                ggplot(aes( y = RMSE, x = Year, fill = Regularization )) + # 
                facet_wrap( ~ Method, nrow = 1) +
                geom_boxplot(
                    lwd = 1.5, 
                    fatten = 0.5, 
                    alpha = 0.5 
                ) + 
                geom_hline(yintercept=1, 
                           linetype="dashed", 
                           color = "black", 
                           size = rel(0.5),
                           alpha = 0.7) + #
                #ylim(0, 2) +
                ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
                xlab(TeX('Year')) + 
                # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
                scale_fill_manual(values = c("red", "blue")) + #"red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
                theme_classic(base_size = 12) +
                ylim(0.4, 1.1) + 
                theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                       axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                       axis.title = element_text(face="bold", color="black", size=rel(1.75)),
                       legend.key.size = unit(2, "line"), # added in to increase size
                       legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
                       legend.title = element_text(face="bold", color="black", size = rel(2)),
                       strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
                       #legend.position ="none"
                ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months")) +  
                scale_x_discrete(breaks = seq(2003,2019)[c(1, 5,10,15, 17)]) 
            
            
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_OldvsNew_", mnYr, ".png"),
                    plot = plt,
                    width = 22,
                    height = 7 
            )
    }
}
e <- cbind(a, a$oec_country0 / a$stacking_country_zeroOut)
colnames(e)[ncol(e)] <- "quant"

e <- cbind(e, a$oec_country0 / a$country_noLinear)
colnames(e)[ncol(e)] <- "noLinear"

e <- cbind(e, a$oec_country2 / a$stacking_country)
colnames(e)[ncol(e)] <- "spec"

# performance of zero out compared to stacking zero out method
e %>% tibble %>% 
    dplyr::group_by( country.1) %>% 
    #dplyr::filter( testYear > 2004) %>% 
    dplyr::summarize(my_mean = mean(quant) ) %>% print(n = Inf)

# cbind(e1$my_mean[-c(1,6,24)], e2$my_mean)
# performance of zero out compared to no country-specific model with no linear component
e %>% tibble %>% 
    dplyr::group_by( country.1, testYear) %>% 
    dplyr::summarize(my_mean = mean(noLinear) ) %>% print(n = Inf)

d <- e %>% tibble %>% 
    dplyr::group_by( country.1) %>% 
    dplyr::summarize(my_mean = mean(noLinear) ) 

boxplot(d$my_mean, ylim = c(0.93, 1.01))

d <- e %>% tibble %>% 
    dplyr::group_by( country.1) %>% 
    dplyr::summarize(my_mean = mean(quant) ) 

d %>% print(n = Inf)

boxplot(d$my_mean, ylim = c(0.93, 1.01))


# specialsit
d <- e %>% tibble %>% 
    dplyr::group_by( country.1) %>% 
    dplyr::summarize(my_mean = mean(spec) ) 

d %>% print(n = Inf)

d <- e %>% tibble %>% 
    dplyr::group_by( country.1) %>% 
    dplyr::summarize(my_mean = mean(quant) ) 

boxplot(d$my_mean, ylim = c(0.93, 1.01))

# looking at zero out on specific year
e %>% tibble %>% 
    filter(testYear == 2001) %>%
    dplyr::group_by( country.1) %>% 
    dplyr::summarize(my_mean = mean(quant) ) %>% print(n = Inf)


#######################
# All-together Figures AND Tables
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)

minYrs <- c(100) # actually in months
nH <- FALSE
# cycle through minimum years and make a different figure for each
for(mnYr in minYrs){
    
    if(mnYr == 150){
        yrs <- 2003:2019
    }else{
        yrs <- 2003:2019
    }
    tnVec <- c("zero", "cvCF") # 
    proVec <- c(4) # linear term for time
    
    # matrix for averages for tables
    avgMat <- matrix(nr = 4, nc = length(yrs))
    colnames(avgMat) <- yrs
    
    for(tn in tnVec){
        
        for(xPro in proVec){
            
            setwd("~/Desktop/Research/mort17")
            a <- read.csv( paste0("mrtNew_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nH",nH,"_mnSt_-01-01_mnTr_", mnYr, "_fld_5_smpSzWt_1") )
            
            #mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1
            
            setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            par(mfrow = c(1,3))
            
            mat <- mat1 <- mat2 <- mat3 <- matrix(NA, nrow(a), ncol = length(yrs))
            
            for(j in 1:length(yrs)){
                mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ] # oec specialist
                mat1[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$stacking_country_zeroOut[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ] # zero out specialist
                mat2[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ] # oec zero out
                mat3[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$stacking_country[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ] # # specialist
            }
            
            for(l in 1:4){
                num <- ifelse(l == 1, "", l - 1)
                avgMat[l,] <- colMeans( get(paste0("mat", num)), na.rm = TRUE )
            }
            avgMat <- rbind( avgMat[4, ], avgMat[-4, ] ) # put in right order
            rownames(avgMat) <- c("Specialist", "OEC", "Zero Out", "OEC Zero Out")
            avgMat <- round(avgMat, 3)
            
            print(paste0("Tune: ", tn, "; mnYr ", mnYr))
            print( kable(avgMat,  format = "latex", booktabs = T) %>% kable_styling(position = "center") )
            
            mat <- as.data.frame(mat)
            colnames(mat) <- yrs
            m <-cbind( gather(mat), "Specialist OEC")
            
            mat1 <- as.data.frame(mat1)
            colnames(mat1) <- yrs
            m1 <- cbind(gather(mat1), "Zero Out")
            
            mat2 <- as.data.frame(mat2)
            colnames(mat2) <- yrs
            m2 <- cbind(gather(mat2), "Zero Out OEC")
            
            mat3 <- as.data.frame(mat3)
            colnames(mat3) <- yrs
            m3 <- cbind(gather(mat3), "Specialist")
            
            colnames(m2)[3] <- colnames(m1)[3]<- colnames(m)[3] <- colnames(m3)[3] <- "Method"
            
            m <- rbind(m, m1, m2, m3)
            
            t <- ifelse(tn == "zero", "OLS", "Ridge")
            
            m <- cbind(m, factor(t) )
            colnames(m) <- c("Year", "RMSE", "Method", "Regularization")
            
            if(tn == "zero"){
                df <- m
            }else{
                # bind them together
                df <- rbind(df, m)
            }
            
        }
    }
    df <- df[!is.na(df$RMSE),]
    df$Year <- as.factor(df$Year)
    
    # ols
    plt = df %>% tibble %>% 
        # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
        # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
        # dplyr::filter(cl == 3) %>%
        # dplyr::filter(n == 300) %>%
        dplyr::filter(Regularization == "OLS") %>%
        ggplot(aes( y = RMSE, x = Year, fill = Method )) + # 
        # facet_wrap( ~ Method, nrow = 1) +
        geom_boxplot(
            lwd = 1.5, 
            fatten = 0.5, 
            alpha = 0.5 
        ) + 
        geom_hline(yintercept=1, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7) + #
        #ylim(0, 2) +
        ylab(TeX('$\\mathbf{RMSE/RMSE_{Country}}$') )+ 
        xlab(TeX('Year')) + 
        # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
        scale_fill_manual(values = c("red", "blue", "green", "#0868ac")) + # , "#E69F00" 
        theme_classic(base_size = 12) +
        ylim(0.1, 1.1) + 
        theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
               axis.text=element_text(face="bold",color="black", size=rel(1.75)),
               axis.title = element_text(face="bold", color="black", size=rel(1.75)),
               legend.key.size = unit(2, "line"), # added in to increase size
               legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
               legend.title = element_text(face="bold", color="black", size = rel(2)),
               strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
               #legend.position ="none"
        ) #+ ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months"))
    
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_OLS_Together", mnYr, ".pdf"),
            plot = plt,
            width = 16,
            height = 6
    )
    
    # ridge
    plt = df %>% tibble %>% 
        # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
        # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
        # dplyr::filter(cl == 3) %>%
        # dplyr::filter(n == 300) %>%
        dplyr::filter(Regularization == "Ridge") %>%
        ggplot(aes( y = RMSE, x = Year, fill = Method )) + # 
        # facet_wrap( ~ Method, nrow = 1) +
        geom_boxplot(
            lwd = 1.5, 
            fatten = 0.5, 
            alpha = 0.5 
        ) + 
        geom_hline(yintercept=1, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7) + #
        #ylim(0, 2) +
        ylab(TeX('$\\mathbf{RMSE/RMSE_{Country}}$') )+ 
        xlab(TeX('Year')) + 
        # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
        scale_fill_manual(values = c("red", "blue", "green", "#0868ac")) + # , "#E69F00" 
        theme_classic(base_size = 12) +
        ylim(0.1, 1.1) + 
        theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
               axis.text=element_text(face="bold",color="black", size=rel(1.75)),
               axis.title = element_text(face="bold", color="black", size=rel(1.75)),
               legend.key.size = unit(2, "line"), # added in to increase size
               legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
               legend.title = element_text(face="bold", color="black", size = rel(2)),
               strip.text.x = element_text(face="bold", color="black", size = rel(2)) #,
               #legend.position ="none"
        ) #+ ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months"))
    
    
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_Ridge_Together", mnYr, ".pdf"),
            plot = plt,
            width = 16,
            height = 6
    )
}
