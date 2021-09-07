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
monthsVec <- c(12, 24, 36, 48, 60)
tnVec <- c( "zero", "cvCF") 
cnt = 0
for(months in monthsVec){
    for(tn in tnVec){
        cnt <- cnt + 1
        
        setwd("~/Desktop/Research/mort17")
        # setwd("~/Desktop/OEC/Manuscript Figures/COVID/mort17")
        fileNm <- paste0("mrtNew_eta0.5_xPro_4_etaTn_cv_etaSpec_cvSpec_Wcv_cv_Wspec_FALSE_Wspec0_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE_hr_10_oecTn_",tn, "_stdTn_", tn, "_tnOr_FALSE_tstCntMn_", months, ".nHFALSE_mnSt_-01-01_mnTr_100_fld_5_smpSzWt_1")
        if(file.exists(fileNm)){
            
            a <- read.csv( fileNm )
            
            # remove southern hemisphere
            a <- a %>% tibble %>% dplyr::filter( !country.1 %in% c("New Zealand", "Chile", "Australia DCD") )
            
            # setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
            par(mfrow = c(1,1))
            
            mat <- mat2 <- mat3 <- mat4 <- matrix(NA, nrow(a), ncol = length(yrs))
            for(j in 1:length(yrs)){
                mat[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$stacking_country[a$testYear == yrs[j] ]
                mat2[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$stacking_country_zeroOut[a$testYear == yrs[j] ]
                
                mat3[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country2[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
                mat4[1:length(a$oec_country2[a$testYear == yrs[j] ]),j] <- a$oec_country0[a$testYear == yrs[j] ] / a$country[a$testYear == yrs[j] ]
                
            }
            
            mat <- as.data.frame(mat)
            colnames(mat) <- yrs
            m <-cbind( gather(mat), "Specialist")
            
            mat2 <- as.data.frame(mat2)
            colnames(mat2) <- yrs
            m1 <- cbind(gather(mat2), "Zero Out")
            
            mat3 <- as.data.frame(mat3)
            colnames(mat3) <- yrs
            m3 <- cbind(gather(mat3), "Spec/country")
            
            mat4 <- as.data.frame(mat4)
            colnames(mat4) <- yrs
            m4 <- cbind(gather(mat4), "ZeroOut/country")
            
            colnames(m)[3] <- colnames(m1)[3] <- colnames(m3)[3] <- colnames(m4)[3] <- "Method"
            m <- rbind(m, m1, m3, m4)
            
            colnames(m)[3] <- "Method"
            
            t <- ifelse(tn == "zero", "OLS", "Ridge")
            print(months)
            m <- cbind(m, factor(t), months )
            colnames(m) <- c("Year", "RMSE", "Method", "Regularization", "Months")
            
            if(cnt == 1){
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

# rename labels
supp.labs <- c(expression(bold(paste("OEC"^"S", " vs. MSS"^"S"))), 
               expression(bold(paste("OEC"^"S", " vs. SSM")))
) 
names(supp.labs) <- c("Specialist", "Spec/country")

plt <- 
    df %>% 
    as_tibble() %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    dplyr::filter(Year >= 2010) %>%
    filter(Method %in% c("Specialist", "Spec/country")) %>%
    mutate( Months = factor(Months ) ) %>%
    mutate( Method = factor(Method, labels= supp.labs ) ) %>%
    ggplot(aes(y = RMSE, x = Months, fill = Regularization, color = Regularization)) + 
    facet_wrap(~ Method,
               labeller=label_parsed
    ) +
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
    labs(y        = TeX('$\\mathbf{RMSE_{OEC^S}/RMSE_{Method}}$'),
         x        = "Months of Training Data") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.direction = "horizontal",
          text = element_text(face = "bold"))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/Final Covid Figures/Mort17")
ggsave( paste("Mortality_Specialist_ZeroOut_OLS_Ridge_New_", mnYr, "Months_Spec.pdf"),
        plot = plt,
        width  = 8,
        height = 4)

