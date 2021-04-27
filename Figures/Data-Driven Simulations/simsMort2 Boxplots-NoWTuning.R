# No W tuning for specialist and 100 min training for training countriesa

# generalist specialist, zero out
######################################################
# compare different K and betaVar -- no regularization 
######################################################
setwd("~/Desktop/Research/simsMort2")
library(latex2exp)
bVar <- c(0.25, 0.5, 1, 1.5, 3) # beta variances
tune <- c("zero")
minTime <- 100
clust <- 0 # no clusters
Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")
ls_K <- vector(length = length(Kvec), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(b in 1:length(bVar)){
    for(k in 1:length(Kvec)){
        cnt <- cnt + 1 # counter
        K <- Kvec[k]
        tn <- tune[1]
        x <- bVar[b]
        flNm <- paste0("oecSpc.Mort0.5.etTn_TRUE_cv_stTn_cv_oecTn_", tn, "_stTnIn_TRUE_stCV_cv_minTm_",minTime,"_AvgW_FALSE_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_sampSz_NA_smpSzW_1_tstTrn_52_minTrn_", minTime, "_bVar_", x, "_reTn_FALSE_eta_TRUE_K_", K, "_s_511")
        
        if(file.exists(flNm)){
            # check to see if file exists
            d <- read.csv(flNm) 
            rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                  Specialist = d$oec_country2 / d$stacking_country,
                                  ZeroOut = d$oec_country0 / d$stacking_country_zeroOut
            )
            
            # rmseMat2 <- data.frame(Generalistoec = d$oec_test2 ,
            #                       Specialistoec = d$oec_country2 ,
            #                       ZeroOutoec = d$oec_country0,
            #                       stacking = d$stack2,
            #                       spec = d$stacking_country,
            #                       d$stacking_country_zeroOut
            #                       
            # )
            # print( colMeans(rmseMat2) )
            
            ls_beta[[cnt]] <- cbind( 
                gather(rmseMat), 
                K, 
                x
            )
            
        }
        
        
        
    }
}

dat <- do.call(rbind, ls_beta)
dat$K <- as.factor(dat$K)
dat$x <- as.factor(dat$x)
# all three of gen/spec/spec0
ggplot(dat, aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",


plt1 <- ggplot(dat[dat$key == "Specialist",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Specialist')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.4, 1.5) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

# zero out
plt0 <- ggplot(dat[dat$key == "ZeroOut",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Zero Out')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.55, 1.2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

# generalist
plt2 <- ggplot(dat[dat$key == "Generalist",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Generalist')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.55, 1.2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",


# save figures
setwd("~/Desktop/Research Final/Mortality/Figures/simsMort2 Figures")
ggsave("mortSims_spec_NoWTrain100.pdf", 
       plot = plt1,
       width = 10, height = 8)

ggsave("mortSims_spec0_NoWTrain100.pdf", 
       plot = plt0,
       width = 10, height = 8)


ggsave("mortSims_generalist_NoWTrain100.pdf", 
       plot = plt2,
       width = 10, height = 8)


##########################################################################
# No W tuning for specialist and 150 min training for training countriesa

# generalist specialist, zero out
######################################################
# compare different K and betaVar -- WITH regularization
######################################################
setwd("~/Desktop/Research/simsMort2")
library(latex2exp)
bVar <- c(0.25, 0.5, 1, 1.5, 3) # beta variances
tune <- c("cvCF")
minTime <- 100
clust <- 0 # no clusters
Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")
ls_K <- vector(length = length(Kvec), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(b in 1:length(bVar)){
    for(k in 1:length(Kvec)){
        cnt <- cnt + 1 # counter
        K <- Kvec[k]
        tn <- tune[1]
        x <- bVar[b]
        flNm <- paste0("oecSpc.Mort0.5.etTn_TRUE_cv_stTn_cv_oecTn_", tn, "_stTnIn_TRUE_stCV_cv_minTm_",minTime,"_AvgW_FALSE_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_sampSz_NA_smpSzW_1_tstTrn_52_minTrn_", minTime, "_bVar_", x, "_reTn_FALSE_eta_TRUE_K_", K, "_s_511")
        
        if(file.exists(flNm)){
            # check to see if file exists
            d <- read.csv(flNm) 
            rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                  Specialist = d$oec_country2 / d$stacking_country,
                                  ZeroOut = d$oec_country0 / d$stacking_country_zeroOut
            )
            
            # rmseMat2 <- data.frame(Generalistoec = d$oec_test2 ,
            #                       Specialistoec = d$oec_country2 ,
            #                       ZeroOutoec = d$oec_country0,
            #                       stacking = d$stack2,
            #                       spec = d$stacking_country,
            #                       d$stacking_country_zeroOut
            #                       
            # )
            # print( colMeans(rmseMat2) )
            
            ls_beta[[cnt]] <- cbind( 
                gather(rmseMat), 
                K, 
                x
            )
            
        }
        
        
        
    }
}

dat <- do.call(rbind, ls_beta)
dat$K <- as.factor(dat$K)
dat$x <- as.factor(dat$x)
# all three of gen/spec/spec0
ggplot(dat, aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",


plt1 <- ggplot(dat[dat$key == "Specialist",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Specialist')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.4, 1.5) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

# zero out
plt0 <- ggplot(dat[dat$key == "ZeroOut",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Zero Out')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.55, 1.2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

# generalist
plt2 <- ggplot(dat[dat$key == "Generalist",], aes( y = value, x = x, fill = K )) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    ) +
    # facet_grid( .~  key) +
    ggtitle( TeX(sprintf(paste0('Generalist')))) + 
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.55, 1.2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",



# save figures
setwd("~/Desktop/Research Final/Mortality/Figures/simsMort2 Figures")
ggsave("mortSims_spec_NoWTrain100cvCF.pdf", 
       plot = plt1,
       width = 10, height = 8)

ggsave("mortSims_spec0_NoWTrain100cvCF.pdf", 
       plot = plt0,
       width = 10, height = 8)


ggsave("mortSims_generalist_NoWTrain100cvCF.pdf", 
       plot = plt2,
       width = 10, height = 8)



