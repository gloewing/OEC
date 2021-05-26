#######################
# ggplot2 Specialist
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)
nhVec <- c(FALSE)
mnYr <- 100
yrs <- 2003:2019
tnVec <- c( "zero", "cvCF") 
for(tn in tnVec){
    
    setwd("~/Desktop/Research/mort17")
    # setwd("~/Desktop/OEC/Manuscript Figures/COVID/mort17")
    fileNm <- paste0("mrtNew_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-01-01_mnTr_100_fld_5_smpSzWt_1")
    if(file.exists(fileNm)){
        
        a <- read.csv( fileNm )
        
        # remove southern hemisphere
        a <- a %>% tibble %>% dplyr::filter( !country.1 %in% c("New Zealand", "Chile", "Australia DCD") )
        
        # setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
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
df <- df[!is.na(df$RMSE),]
df$Year <- as.factor(df$Year)

plt <- df %>% 
    as_tibble() %>%
    filter(df$Method %in% c("Specialist", "Zero Out")) %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::filter(Year >= 2010) %>%
    mutate(Year = factor(Year)) %>%
    ggplot(aes(y = RMSE, x = Year, fill = Regularization, color = Regularization)) + 
    facet_wrap(~Method) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size          = 0.20,
                 outlier.size  = 0.50) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_color_manual(values = c("#ca0020", "#0868ac")) +
    scale_fill_manual(values = c("#ca0020", "#0868ac")) +
    labs(y        = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$'),
         x        = TeX('$\\mathbf{Year}$')) +
    coord_cartesian(ylim = c(0.3, 1.3)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal",
          text = element_text(face = "bold"))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_", mnYr, ".pdf"),
        plot = plt,
        width  = 8,
        height = 4)

# rename labels
supp.labs <- c("Specialist", "Zero Out")
names(supp.labs) <- c("Spec/noLin", "ZeroOut/noLin")

plt <- df %>% 
    as_tibble() %>% 
    filter(df$Method %in% c("Spec/noLin", "ZeroOut/noLin")) %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::filter(Year >= 2010) %>%
    mutate(Year = factor(Year)) %>%
    ggplot(aes(y = RMSE, x = Year, fill = Regularization, color = Regularization)) + 
    facet_wrap(~ Method,
               labeller = labeller( Method = supp.labs )
               ) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size          = 0.20,
                 outlier.size  = 0.50) +
    geom_boxplot(size          = 0.20,
                 outlier.size  = 0.50) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_color_manual(values = c("#ca0020", "#0868ac")) +
    scale_fill_manual(values = c("#ca0020", "#0868ac")) +
    labs(y        = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{SSM}}$'),
         x        = TeX('$\\mathbf{Year}$')) +
    coord_cartesian(ylim = c(0.75, 1.45)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal",
          text = element_text(face = "bold"))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_NoLinear_", mnYr, ".pdf"),
        plot   = plt,
        width  = 8,
        height = 4)

#######################
# All-together Figures AND Tables
#######################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(kableExtra)

minYrs <- mnYr <- 100 # in months
nH <- FALSE
yrs <- 2003:2019
tnVec <- c("zero", "cvCF") 
proVec <- 4 

# matrix for averages for tables
avgMat <- matrix(nr = 4, nc = length(yrs))
colnames(avgMat) <- yrs

for(tn in tnVec){
    setwd("~/Desktop/Research/mort17")
    #setwd("~/Desktop/OEC/Manuscript Figures/COVID/mort17")
    a <- read.csv( paste0("mrtNew_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nH",nH,"_mnSt_-01-01_mnTr_100_fld_5_smpSzWt_1") )
    # setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    par(mfrow = c(1,3))
    
    # full data
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
    
    # subset from 2010 - 2019
    print(paste0("Subset Tune: ", tn, "; mnYr ", mnYr))
    print( kable(avgMat[,-seq(1, 7)],  format = "latex", booktabs = T) %>% kable_styling(position = "center") )
    
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

df <- df[!is.na(df$RMSE),]
df$Year <- as.factor(df$Year)

# ols
plt <- df %>%
    as_tibble() %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::filter(Year >= 2010) %>%
    mutate(Year = factor(Year)) %>%
    dplyr::filter(Regularization == "OLS") %>%
    ggplot(aes(y = RMSE, x = Year, fill = Method, color = Method)) + 
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size          = 0.20,
                 outlier.size  = 0.50) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    labs(y        = TeX('$\\mathbf{RMSE/RMSE_{SSM}}$'),
         x        = TeX('$\\mathbf{Year}$')
         ) + #,
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                      labels = c( expression(bold(Specialist)), expression(bold(OEC[Spec])), expression(bold("Zero Out")),  expression(bold(OEC[Zero])) )
    ) +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                       labels = c( expression(bold(Specialist)), expression(bold(OEC[Spec])), expression(bold("Zero Out")),  expression(bold(OEC[Zero])) )
    ) +
    coord_cartesian(ylim = c(0.1, 1.1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(face = "bold")) 
    
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_OLS_Together", mnYr, ".pdf"),
            plot   = plt,
            width  = 6,
            height = 3)

# ridge
plt <- df %>%
        as_tibble() %>%
        mutate(Year = as.numeric(as.character(Year))) %>%
        dplyr::filter(Year >= 2010) %>%
        mutate(Year = factor(Year)) %>%
        dplyr::filter(Regularization == "Ridge") %>%
        ggplot(aes(y = RMSE, x = Year, fill = Method, color = Method)) + 
        geom_hline(yintercept = 1,
                   linetype   = 2,
                   color      = "darkgray") +
        geom_boxplot(size          = 0.20,
                     outlier.size  = 0.50) +
        geom_boxplot(size          = 0.20,
                     color         = "black",
                     outlier.color = NA) +
        labs(y        = TeX('$\\mathbf{RMSE/RMSE_{SSM}}$'),
             x        = TeX('$\\mathbf{Year}$')
             ) + #,
        scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                          labels = c( expression(bold(Specialist)), expression(bold(OEC[Spec])), expression(bold("Zero Out")),  expression(bold(OEC[Zero])) )
        ) +
        scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                           labels = c( expression(bold(Specialist)), expression(bold(OEC[Spec])), expression(bold("Zero Out")),  expression(bold(OEC[Zero])) )
        ) +
        coord_cartesian(ylim = c(0.1, 1.1)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(face = "bold")) 
    
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_Ridge_Together", mnYr, ".pdf"),
            plot   = plt,
            width  = 6,
            height = 3)
    
########################################
# OEC vs. SSM BOTH without linear terms
########################################

    library(dplyr)
    library(tidyverse)
    library(ggplot2)
    library(latex2exp)
    library(kableExtra)
    nhVec <- c(FALSE)
    mnYr <- 100
    yrs <- 2003:2019
    tnVec <- c( "zero", "cvCF") 
    for(tn in tnVec){
        
        setwd("~/Desktop/Research/mort17")
        # setwd("~/Desktop/OEC/Manuscript Figures/COVID/mort17")
        fileNm <- paste0("mrtNew_eta0.5_xPro_5_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_", tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_12.nHFALSE_mnSt_-01-01_mnTr_100_fld_5_smpSzWt_1")
        if(file.exists(fileNm)){
            
            a <- read.csv( fileNm )
            
            # remove southern hemisphere
            a <- a %>% tibble %>% dplyr::filter( !country.1 %in% c("New Zealand", "Chile", "Australia DCD") )
            
            # setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            par(mfrow = c(1,1))
            
            mat <- mat2 <- mat3 <- mat4 <- mat5 <- mat6 <- matrix(NA, nrow(a), ncol = length(yrs))
            mat9 <- mat10 <- mat7 <- mat8 <- mat
            for(j in 1:length(yrs)){
                mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
                mat2[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
                
                mat3[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
                mat4[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
                
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
    df <- df[!is.na(df$RMSE),]
    df$Year <- as.factor(df$Year)
    
##################################################################
# OEC vs. stacking when both do not have linear terms
    
    plt <- df %>% 
        as_tibble() %>%
        filter(df$Method %in% c("Specialist", "Zero Out")) %>%
        mutate(Year = as.numeric(as.character(Year))) %>%
        dplyr::filter(Year >= 2010) %>%
        mutate(Year = factor(Year)) %>%
        ggplot(aes(y = RMSE, x = Year, fill = Regularization, color = Regularization)) + 
        facet_wrap(~Method) +
        geom_hline(yintercept = 1,
                   linetype   = 2,
                   color      = "darkgray") +
        geom_boxplot(size          = 0.20,
                     outlier.size  = 0.50) +
        geom_boxplot(size          = 0.20,
                     color         = "black",
                     outlier.color = NA) +
        scale_color_manual(values = c("#ca0020", "#0868ac")) +
        scale_fill_manual(values = c("#ca0020", "#0868ac")) +
        labs(y        = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stacking}}$'),
             x        = TeX('$\\mathbf{Year}$')) +
        coord_cartesian(ylim = c(0.7, 1.25)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "top",
              legend.direction = "horizontal",
              text = element_text(face = "bold"))
    
    # both stacking ond oec have no linear term for time
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_NoLinOEC_vs_NoLinStack.pdf"),
            plot = plt,
            width  = 8,
            height = 4)
    
    
    # OEC with no linear term and country-specific model with no linear term
    # rename
    supp.labs <- c("Specialist", "Zero Out")
    names(supp.labs) <- c("Spec/noLin", "ZeroOut/noLin")
    
    
    plt1 <- df %>% 
        as_tibble() %>% 
        filter(df$Method %in% c("Spec/noLin", "ZeroOut/noLin")) %>%
        mutate(Year = as.numeric(as.character(Year))) %>%
        dplyr::filter(Year >= 2010) %>%
        mutate(Year = factor(Year)) %>%
        ggplot(aes(y = RMSE, x = Year, fill = Regularization, color = Regularization)) + 
        facet_wrap(~ Method,
                   labeller = labeller( Method = supp.labs )
        ) +
        geom_hline(yintercept = 1,
                   linetype   = 2,
                   color      = "darkgray") +
        geom_boxplot(size          = 0.20,
                     outlier.size  = 0.50) +
        geom_boxplot(size          = 0.20,
                     outlier.size  = 0.50) +
        geom_boxplot(size          = 0.20,
                     color         = "black",
                     outlier.color = NA) +
        scale_color_manual(values = c("#ca0020", "#0868ac")) +
        scale_fill_manual(values = c("#ca0020", "#0868ac")) +
        labs(y        = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{SSM}}$'),
             x        = TeX('$\\mathbf{Year}$')) +
        coord_cartesian(ylim = c(0.7, 1.25)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "top",
              legend.direction = "horizontal",
              text = element_text(face = "bold"))
    
    # both stacking ond country have no linear term for time
    setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
    ggsave( paste("Mortality_NoLinOEC_vs_NoLinCountry.pdf"),
            plot = plt1,
            width  = 8,
            height = 4)
    
