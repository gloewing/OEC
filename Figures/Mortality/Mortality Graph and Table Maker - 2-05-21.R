# "mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_spec_alph0_FALSE_glm_FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1"
# good
# mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_zero_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1

# VERY GOOD BELOW
# mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_FALSE_eTnAvgW_FALSE_Wcv_cv_Wspec_FALSE_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE.pca.FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1
# mortTm_eta0.99_spInt_TRUE_xPro_3_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_FALSE_eTnAvgW_FALSE_Wcv_cv_Wspec_FALSE_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE.pca.FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1

# ----
# mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_5
# "mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_spec_alph0_FALSE_glm_FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1"
# good
# mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_zero_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1

###### Specialist-7_min50    550 500 Linear No regularization
#######################
# ggplot2 Specialist
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)

minYrs <- c(100, 150) # actually in months

# cycle through minimum years and make a different figure for each
for(mnYr in minYrs){
    
if(mnYr == 150){
    yrs <- 2014:2020
}else{
    yrs <- 2013:2020
}
tnVec <- c("zero", "cvCF") # 
proVec <- c(4) # linear term for time

for(tn in tnVec){
    
    for(xPro in proVec){
        
        setwd("~/Desktop/Research/mort14")
        a <- read.csv(paste0("mrtTmA_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_", mnYr, "_smpSzWt_1"))
        
        #mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1
        
        setwd("~/Desktop/Research Final/Mortality/Figures")
        par(mfrow = c(1,3))
        
        mat <- mat2 <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
            mat2[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
            
        }
        
        mat <- as.data.frame(mat)
        colnames(mat) <- yrs
        m <-cbind( gather(mat), "Specialist")
        
        mat2 <- as.data.frame(mat2)
        colnames(mat2) <- yrs
        m1 <- cbind(gather(mat2), "Zero Out")
        colnames(m)[3] <- colnames(m1)[3] <- "Method"
        
        m <- rbind(m, m1)
        
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
        
        plt = df %>% tibble %>% 
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
            ) + ggtitle(paste0("OEC: ", mnYr, " Minimum Training Months"))
        
        
        setwd("~/Desktop/Research Final/Mortality/Figures")
        ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge", mnYr, ".png"),
                plot = plt,
                width = 18,
                height = 7 
        )
}
######################################################################################################################
#######################
# All-together Figures AND Tables
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)

minYrs <- c(100, 150) # actually in months

# cycle through minimum years and make a different figure for each
for(mnYr in minYrs){
    
    if(mnYr == 150){
        yrs <- 2014:2020
    }else{
        yrs <- 2013:2020
    }
    tnVec <- c("zero", "cvCF") # 
    proVec <- c(4) # linear term for time
    
    # matrix for averages for tables
    avgMat <- matrix(nr = 4, nc = length(yrs))
    colnames(avgMat) <- yrs
    
    for(tn in tnVec){
        
        for(xPro in proVec){
            
            setwd("~/Desktop/Research/mort14")
            a <- read.csv(paste0("mrtTmA_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_", mnYr, "_smpSzWt_1"))
            
            #mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1
            
            setwd("~/Desktop/Research Final/Mortality/Figures")
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
    
    setwd("~/Desktop/Research Final/Mortality/Figures")
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
    
    
    setwd("~/Desktop/Research Final/Mortality/Figures")
    ggsave( paste("Mortality_Ridge_Together", mnYr, ".pdf"),
            plot = plt,
            width = 16,
            height = 6
    )
}

######################################################################################################################
        #######################
        # ggplot2 Zero Out
        #######################
        library(dplyr)
        library(tidyverse)
        library(ggplot2)
        library(latex2exp)
        minYrs <- 150 # actually in months
        
        if(minYrs == 150){
            yrs <- 2014:2020
        }else{
            yrs <- 2013:2020
        }
        tnVec <- c("zero", "cvCF") # 
        proVec <- c(4)
        
        for(tn in tnVec){
            
            for(xPro in proVec){
                
                setwd("~/Desktop/Research/mort14")
                a <- read.csv(paste0("mrtTmA_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_", minYrs, "_smpSzWt_1"))
                
                #mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1
                
                setwd("~/Desktop/Research Final/Mortality/Figures")
                par(mfrow = c(1,3))
                
                mat <- matrix(NA, nrow(a), ncol = length(yrs))
                
                for(j in 1:length(yrs)){
                    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
                    
                }
                
                mat <- as.data.frame(mat)
                colnames(mat) <- yrs
                m <- gather(mat)
                
                t <- ifelse(tn == "zero", "OLS", "Ridge")
                
                m <- cbind(m, factor(t) )
                colnames(m) <- c("Year", "RMSE", "Method")
                
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
        
        plt = df %>% tibble %>% 
            # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
            # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
            # dplyr::filter(cl == 3) %>%
            # dplyr::filter(n == 300) %>%
            #dplyr::filter(Method == "OLS") %>%
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
            ) + ggtitle("Zero Out")
        
        
        plt2 = m %>% tibble %>% 
            # dplyr::filter(key %in% c("Generalist", "Merged")) %>%
            # dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
            # dplyr::filter(cl == 3) %>%
            # dplyr::filter(n == 300) %>%
            #dplyr::filter(Method == "OLS") %>%
            ggplot(aes( y = RMSE, x = Year )) + # , fill = Year
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
            ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$') )+ 
            xlab(TeX('Year')) + 
            # scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) + # 
            # scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00", "#525252", "grey")) +
            theme_classic(base_size = 12) +
            ylim(0.4, 1.1) + 
            theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
                   axis.text=element_text(face="bold",color="black", size=rel(1.75)),
                   axis.title = element_text(face="bold", color="black", size=rel(1.5)),
                   strip.text.x = element_text(face="bold", color="black", size = rel(2)),
                   legend.position ="none"
            ) 
        
        
        setwd("~/Desktop/Research Final/Mortality/Figures")
        ggsave( "Mortality_Specialist_OLS_Ridge.png",
                plot = plt,
                width = 10,
                height = 7 
        )        
        
        
        
        
        
######################################################################################################################
        # 
        # boxplot(mat, yrs,na.rm = T,
        #         ylab = "RMSE / RMSE_Country",
        #         names = yrs,
        #         xlab = "Test Year",
        #         main = "Specialist By Year")
        # abline(h = 1)
        
        colMeans(mat, na.rm=T)
        
        
        mat <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
            
        }
        
        
        boxplot(mat, yrs,na.rm = T,
                ylab = "RMSE / RMSE_Country_0",
                names = yrs,
                xlab = "Test Year",
                main = "Zero Out Specialist By Year")
        abline(h = 1)
        colMeans(mat, na.rm=T)
        
        mat <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
            
        }
        
        
        boxplot(mat, yrs,na.rm = T,
                ylab = "RMSE / RMSE_Stacking",
                names = yrs,
                xlab = "Test Year",
                main = "OEC Generalist vs. Stacking")
        abline(h = 1)
        
        colMeans(mat, na.rm = T)
        dev.off()
        
    }
    
    
}



###### Specialist-7_min50    550 500 Linear No regularization
#######################
# Standard R Plot
#######################
minYrs <- 150 # actually in months

if(minYrs == 150){
    yrs <- 2014:2020
}else{
    yrs <- 2013:2020
}
tnVec <- c("zero", "cvCF")
proVec <- c(1, 4)

for(tn in tnVec){
    
    for(xPro in proVec){
        
        setwd("~/Desktop/Research/mort14")
        a <- read.csv(paste0("mrtTmA_eta0.5_xPro_", xPro, "_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_", minYrs, "_smpSzWt_1"))
        
        #mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1
        
        setwd("~/Desktop/Research Final/Mortality/Figures")
        jpeg(paste0("oec_mort14_xPro", xPro, "_min_", minYrs, "_tn_", tn, ".jpg"), width = 850, height = 300)
        
        par(mfrow = c(1,3))
        
        mat <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
            
        }
        boxplot(mat, yrs,na.rm = T,
                ylab = "RMSE / RMSE_Country",
                names = yrs,
                xlab = "Test Year",
                main = "Specialist By Year")
        abline(h = 1)
        
        colMeans(mat, na.rm=T)
        
        
        mat <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
            
        }
        
        
        boxplot(mat, yrs,na.rm = T,
                ylab = "RMSE / RMSE_Country_0",
                names = yrs,
                xlab = "Test Year",
                main = "Zero Out Specialist By Year")
        abline(h = 1)
        colMeans(mat, na.rm=T)
        
        mat <- matrix(NA, nrow(a), ncol = length(yrs))
        
        for(j in 1:length(yrs)){
            mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
            
        }
        
        
        boxplot(mat, yrs,na.rm = T,
                ylab = "RMSE / RMSE_Stacking",
                names = yrs,
                xlab = "Test Year",
                main = "OEC Generalist vs. Stacking")
        abline(h = 1)
        
        colMeans(mat, na.rm = T)
        dev.off()
        
    }
    
    
}


# yrs <- 2011:2020
# mrtTmA_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_zero_stdTn_zero_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-03-01_mnTr_150_smpSzWt_1

### Specialist-7_min50_byYear   550 500















mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear >= yrs[j] ]),j] <- a$oec_country2[a$testYear >= yrs[j] ] / a$stacking_country[a$testYear >= yrs[j] ]
    
    
}
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")
jpeg(paste0("Specialist_Linear_ZeroTn_Agg", ".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 1)")
abline(h = 1, col = "red")
dev.off()


colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2013:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
    
    
}
 ### Specialist-7_min50_byYear   550 500
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")

jpeg(paste0("Specialist_Linear_ZeroTn_Annual",".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
dev.off()

colMeans(mat, na.rm = T)

# ----------------------------------------------------------
# same as above but compare to country stacking
setwd("~/Desktop/Research/mort8")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_zero_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear >= yrs[j] ]),j] <- a$oec_country0[a$testYear >= yrs[j] ] / a$stacking_country_zeroOut[a$testYear >= yrs[j] ]
    
    
}
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")
jpeg(paste0("Specialist_Linear_ZeroTn_Agg_CountryStack", ".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 1)")
abline(h = 1, col = "red")
dev.off()


colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2014:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
    
    
}
### Specialist-7_min50_byYear   550 500
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")

jpeg(paste0("Specialist_Linear_ZeroTn_Annual_CountryStack",".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
dev.off()

colMeans(mat, na.rm = T)

#########################################################################



# ----------------------------------------------------------
# same as above but compare to country stacking
setwd("~/Desktop/Research/mort8")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_zero_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear >= yrs[j] ]),j] <- a$oec2[a$testYear >= yrs[j] ] / a$stacking[a$testYear >= yrs[j] ]
    
    
}
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")
jpeg(paste0("Specialist_Linear_ZeroTn_Agg_CountryStack", ".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 1)")
abline(h = 1, col = "red")
dev.off()


colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2014:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
    
    
}
### Specialist-7_min50_byYear   550 500
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")

jpeg(paste0("Specialist_Linear_ZeroTn_Annual_CountryStack",".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
dev.off()

colMeans(mat, na.rm = T)

#########################################################################

###### Specialist-7_min50    550 500 Linear Regularization
setwd("~/Desktop/Research/mort8")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_5")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear >= yrs[j] ]),j] <- a$oec_country2[a$testYear >= yrs[j] ] / a$country[a$testYear >= yrs[j] ]
    
    
}
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")
jpeg(paste0("Specialist_Linear_CVCTSTn_Agg", ".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 1)")
abline(h = 1, col = "red")
dev.off()


colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
    
    
}
### Specialist-7_min50_byYear   550 500
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")

jpeg(paste0("Specialist_Linear_CVCTSTn_Annual",".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
dev.off()

colMeans(mat, na.rm = T)

# ----------------------------------------------------------
# same as above but compare to country stacking
setwd("~/Desktop/Research/mort8")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_4_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_etaSpec_TRUE_eTnAvgW_FALSE_Wcv_cv_Wspec_cvSpecTS_TScv_TRUE_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCFTS_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_5")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear >= yrs[j] ]),j] <- a$oec_country2[a$testYear >= yrs[j] ] / a$stacking_country[a$testYear >= yrs[j] ]
    
    
}
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")
jpeg(paste0("Specialist_Linear_CVCTSTn_Agg_CountryStack", ".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 1)")
abline(h = 1, col = "red")
dev.off()


colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
    
    
}
### Specialist-7_min50_byYear   550 500
setwd("~/Desktop/Research Final/Mortality/Figures/mort8")

jpeg(paste0("Specialist_Linear_CVCTSTn_Annual_CountryStack",".jpg"), width = 550, height = 500)

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_CountryStack",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
dev.off()

colMeans(mat, na.rm = T)


######################################################



######################
# minimum of 100 observations

######     550 500
setwd("~/Desktop/Research/mort5")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_eTnAvgW_FALSE_Wcv_cv_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_100_smpSzWt_7")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear >= yrs[j] ]),j] <- a$oec_country2[a$testYear >= yrs[j] ] / a$country[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Specialist OEC Aggregated Across Years (Loss 7)")
abline(h = 1, col = "red")
colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2013:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
    
    
}
###   550 500
boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
colMeans(mat, na.rm = T)


######################


######################
# Generalist (min 50)

######     550 500
setwd("~/Desktop/Research/mort5")
a <- read.csv("mortTime_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_eTnAvgW_FALSE_Wcv_cv_sclX_TRUE_alpha0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMon_12_lowLmTn_FALSE_mnSt_-03-01_sampSzWt_7")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear >= yrs[j] ]),j] <- a$oec2[a$testYear >= yrs[j] ] / a$stacking[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Generalist OEC Aggregated Across Years (Loss 7)")
abline(h = 1, col = "red")
colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
    
    
}
### Generalist-7_min50_byYear   550 500
boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Test Year (Years since 2010)",
        main = "Generalist By Years (Loss 7)")
abline(h = 1)
colMeans(mat, na.rm = T)


######################

######################
# Generalist Average (min 50)

######     550 500
setwd("~/Desktop/Research/mort5")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_TRUE_etaTn_TRUE_cv_eTnAvgW_TRUE_Wcv_zero_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_7")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear >= yrs[j] ]),j] <- a$oec2[a$testYear >= yrs[j] ] / a$stacking[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Generalist Avg OEC Aggregated Across Years (Loss 7)")
abline(h = 1, col = "red")
colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
    
    
}
### Generalist-7_min50_byYear   550 500
boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Test Year (Years since 2010)",
        main = "Generalist Avg By Years (Loss 7)")
abline(h = 1)
colMeans(mat, na.rm = T)


######################


######################
# Generalist Average Compared to Standard Avg (with free intercept_ (min 50)

######     550 500
setwd("~/Desktop/Research/mort5")
a <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_TRUE_etaTn_TRUE_cv_eTnAvgW_TRUE_Wcv_zero_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_7")
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear >= yrs[j] ]),j] <- a$oec2[a$testYear >= yrs[j] ] / a$avg[a$testYear >= yrs[j] ]
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Generalist Avg OEC Aggregated Across Years (Loss 7)")
abline(h = 1, col = "red")
colMeans(a[,-24])
colMeans(mat, na.rm = T)

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec2[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$avg[a$testYear == yrs[j] ]
    
    
}
### Generalist-7_min50_byYear   550 500
boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Stacking",
        xlab = "Test Year (Years since 2010)",
        main = "Generalist Avg By Years (Loss 7)")
abline(h = 1)
colMeans(mat, na.rm = T)


############################
###### xGen comparing 1-3   550 500
setwd("~/Desktop/Research/mort5")
a1 <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_1_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_eTnAvgW_FALSE_Wcv_cv_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1")
a2 <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_2_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_eTnAvgW_FALSE_Wcv_cv_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1")
a3 <- read.csv("mortTm_eta0.99_spInt_TRUE_xPro_3_psi_0.Inf_Avg_FALSE_etaTn_TRUE_cv_eTnAvgW_FALSE_Wcv_cv_sclX_TRUE_alph0_FALSE_glm_FALSE_oecTn_cvCF_itr_20_tstCntMn_12_lwLmTn_FALSE_mnSt_-03-01_mnTr_50_smpSzWt_1")

mean(a1$country / a2$country)
mean(a1$country / a3$country)
mean(a2$country / a3$country)

###### Compare Country specific models ########
## approach 1 / 2
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a1$country[a1$testYear >= yrs[j] ]),j] <- a1$country[a$testYear >= yrs[j] ] / a2$country[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_1 / RMSE_2",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Country-Specific Model Approach 1 / 2")
abline(h = 1, col = "red")

## approach 1 / 3

mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a1$country[a1$testYear >= yrs[j] ]),j] <- a1$country[a$testYear >= yrs[j] ] / a3$country[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_1 / RMSE_3",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Country-Specific Model Approach 1 / 3")
abline(h = 1, col = "red")



## approach 2 / 3
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a2$country[a2$testYear >= yrs[j] ]),j] <- a2$country[a$testYear >= yrs[j] ] / a3$country[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_2 / RMSE_3",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "Country-Specific Model Approach 2 / 3",
        ylim = c(0,2))
abline(h = 1, col = "red")




###### Compare Specialist OEC specific models ########
## approach 1 / 2
yrs <- 2011:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a1$oec_country2[a1$testYear >= yrs[j] ]),j] <- a1$oec_country2[a$testYear >= yrs[j] ] / a2$oec_country2[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_1 / RMSE_2",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "OEC Specialist Model Approach 1 / 2")
abline(h = 1, col = "red")

## approach 1 / 3

mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a1$oec_country2[a1$testYear >= yrs[j] ]),j] <- a1$oec_country2[a$testYear >= yrs[j] ] / a3$oec_country2[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_1 / RMSE_3",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "OEC Specialist Model Approach 1 / 3")
abline(h = 1, col = "red")



## approach 2 / 3
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a2$testYear >= yrs[j] ]),j] <- a2$oec_country2[a$testYear >= yrs[j] ] / a3$oec_country2[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE_2 / RMSE_3",
        xlab = "Starting From Test Year (Years since 2010)",
        main = "OEC Specialist Model Approach 2 / 3",
        ylim = c(0,2))
abline(h = 1, col = "red")


colMeans(a[,-24])
colMeans(mat, na.rm = T)
# ------------------------

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
    
    
}
### Specialist-7_min50_byYear   550 500
boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)",
        main = "Specialist By Years (Loss 7)")
abline(h = 1)
colMeans(mat, na.rm = T)


######################

######################
mat <- matrix(NA, nrow(a), ncol = length(yrs))
for(j in 1:length(yrs)){
    mat[1:length(a$oec[a$testYear >= yrs[j] ]),j] <- a$oec2[a$testYear >= yrs[j] ] / a$stacking[a$testYear >= yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = " >= Test Year (Years since 2010)")
abline(h = 1)
colMeans(a[,-24])

yrs <- 2012:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)")
abline(h = 1)
colMeans(mat,na.rm = T)





for(j in 1:length(yrs)){
    mat[1:length(a$oec[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$avg[a$testYear == yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)")
abline(h = 1)
colMeans(mat,na.rm = T)



yrs <- 2013:2020
mat <- matrix(NA, nrow(a), ncol = length(yrs))

for(j in 1:length(yrs)){
    mat[1:length(a$oec[a$testYear == yrs[j] ]),j] <- a$oec2[a$testYear == yrs[j] ] / a$stacking[a$testYear == yrs[j] ]
    
    
}

boxplot(mat, yrs,na.rm = T,
        ylab = "RMSE / RMSE_Country",
        xlab = "Test Year (Years since 2010)")
abline(h = 1)
colMeans(mat,na.rm = T)

