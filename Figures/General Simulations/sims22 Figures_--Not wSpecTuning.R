

######################################################
# compare different K and betaVar -- no regularization
######################################################
setwd("~/Desktop/Research/sims22")
library(dplyr)
library(latex2exp)
bVar <- c( 0.0001, 0.05, 0.25, 0.5) # beta variances #  0.001, 0.005, 0.01,
xVar <- c(0, round( sqrt( c(0.01, 0.5, 1.5 ) ), 2) )  # beta variances 
tune <- c("cvCF", "zero", "cv") # , "zero", "cv"
etaSt <- c(TRUE)
sclX <- c(TRUE)
clusts <- c(3,6) # no clusters
e = TRUE
sc = TRUE
nVec <- c(50, 150, 300) # 50,150,
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(cl in clusts){
    for(tn in tune){
        for(n in nVec){
            for(bl in 1:length(bVar)){
                for(xv in 1:length(xVar)){
                    cnt <- cnt + 1 # counter
                    b <- bVar[bl]
                    x <- xVar[xv]
                    
                    flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                    
                    if(x == 0)                    flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", 0, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                    
                    if(file.exists(flNm)){
                        # check to see if file exists
                        d <- read.csv(flNm) 
                        rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                              Specialist = d$oec_country2 / d$stacking_country,
                                              ZeroOut = d$oec_country0 / d$stacking_country_zeroOut,
                                              Merged = d$oec_test2 / d$merge2,
                                              MS = d$stack2 / d$merge2,
                                              SpecOecMerged = d$oec_country2 / d$merge, # the merged tested on specialist
                                              ZeroOecMerged = d$oec_country0 / d$merge,
                                              SpecMerged = d$stacking_country / d$merge, # the merged tested on specialist
                                              ZeroMerged = d$stacking_country_zeroOut / d$merge
                                              
                        )
                        ls_beta[[cnt]] <- cbind( 
                            gather(rmseMat), 
                            x, # put on variance scale
                            cl,
                            b,
                            e,
                            sc,
                            n,
                            tn
                        )
                        
                    }
                    
                    
                    
                }
            }
        }
    }
}

dat <- do.call(rbind, ls_beta)
dat$cl <- as.factor(dat$cl)
dat$x <- round(dat$x^2, 2)
dat$x <- ifelse(dat$x == 1.49, 1.5, dat$x)
dat$x <- as.factor(dat$x)
dat$b <- as.factor(dat$b)
dat$e <- as.factor(dat$e)
dat$n <- as.factor(dat$n)

dat %>% group_by(key, x, b, n, cl, tn) %>% summarize(my_mean = mean(value) ) %>% print(n = Inf)

# Generalist / Merged

# clusters
plt = dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_generalist_merged_CLust_oec.png",
        plot = plt,
        width = 15,
        height = 10 
)

# clusters - plot together all standardized by the merged tested on the specialist country
plt = dat %>% tibble %>% 
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key )) +
    #facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
                      values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title="Method"))

setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_specialist_clusts_oecTogether.png",
        plot = plt,
        width = 15,
        height = 10 
)

# no clusters
plt1 = dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c( "red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_generalist_merged_noCLust_oec.png",
        plot = plt1,
        width = 15,
        height = 10 
)

# clusters - plot together all standardized by the merged tested on the specialist country
plt = dat %>% tibble %>% 
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key )) +
    #facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
                      values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title="Method"))

setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_specialist_Noclusts_oecTogether.png",
        plot = plt,
        width = 15,
        height = 10 
)



plt2 = dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))


setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_specialist_clusts_oec.png",
        plot = plt2,
        width = 15,
        height = 10 
)


plt2 = dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    ylim(0.5, 1.2) + 
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2)),
           strip.text.x = element_text(face="bold", color="black", size = rel(2))
    ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))


setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
ggsave( "sims22_specialist_Noclusts_oec.png",
        plot = plt2,
        width = 15,
        height = 10 
)

# generalist specialist, zero out
######################################################
# compare different K and betaVar -- no regularization
######################################################
setwd("~/Desktop/Research/sims22")
library(latex2exp)
bVar <- c(0.001, 0.05, 0.25) # beta variances
xVar <- c(0, 0.71) # beta variances
tune <- c("cvCF")
etaSt <- c(TRUE, FALSE)
sclX <- c(TRUE, FALSE)
cl <- 3 # no clusters
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(sc in sclX){
    
    for(e in etaSt){
    for(bl in 1:length(bVar)){
    for(xv in 1:length(xVar)){
        cnt <- cnt + 1 # counter
        b <- bVar[bl]
        x <- xVar[xv]
        flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
        
        if(file.exists(flNm)){
            # check to see if file exists
            d <- read.csv(flNm) 
            rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                  Specialist = d$oec_country2 / d$stacking_country,
                                  ZeroOut = d$oec_country0 / d$stacking_country_zeroOut,
                                  GenMerged = d$oec_test2 / d$merge2
            )
            ls_beta[[cnt]] <- cbind( 
                gather(rmseMat), 
                x, 
                cl,
                b,
                e,
                sc
            )
            
        }
        
        
        
    }
}
}
}

dat <- do.call(rbind, ls_beta)
dat$cl <- as.factor(dat$cl)
dat$x <- as.factor(dat$x)
dat$b <- as.factor(dat$b)
dat$e <- as.factor(dat$e)

dat2 <- dat[dat$x == 0.71,]

# all three of gen/spec/spec0
ggplot(dat2[dat2$key == "Specialist",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ e, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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


# all three of gen/spec/spec0
ggplot(dat2[dat2$key == "GenMerged",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ e, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    #ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

ggplot(dat2[dat2$key == "Specialist",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ e, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    #ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

dat %>% dplyr::filter(key == "ZeroOut") %>%
ggplot( aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ e, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    #ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

##### neither scale X or eta but covariate shift is important, clustering is important

# generalist specialist, zero out
######################################################
# compare different K and betaVar -- no regularization
######################################################
setwd("~/Desktop/Research/sims22")
library(dplyr)
library(latex2exp)
bVar <- c(0.001, 0.05, 0.25, 0.5, 1, 1.5) # beta variances
xVar <- c( round( sqrt(c( 0.0001, 0.001, 0.005, 0.01)), 2), 0.71, 1) # beta variances
tune <- c("cv", "zero")
etaSt <- c(TRUE)
sclX <- c(TRUE)
clusts <- c(3,6) # no clusters
e = TRUE
sc = TRUE
nVec <- c(50,150, 300)
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(cl in clusts){
for(tn in tune){
    for(n in nVec){
        for(bl in 1:length(bVar)){
            for(xv in 1:length(xVar)){
                cnt <- cnt + 1 # counter
                b <- bVar[bl]
                x <- xVar[xv]
                flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                
                if(file.exists(flNm)){
                    # check to see if file exists
                    d <- read.csv(flNm) 
                    rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                          Specialist = d$oec_country2 / d$stacking_country,
                                          ZeroOut = d$oec_country0 / d$stacking_country_zeroOut,
                                          GenMerged = d$oec_test2 / d$merge2,
                                          MS = d$stack2 / d$merge2
                    )
                    ls_beta[[cnt]] <- cbind( 
                        gather(rmseMat), 
                        x, 
                        cl,
                        b,
                        e,
                        sc,
                        n,
                        tn
                    )
                    
                }
                
                
                
            }
        }
    }
}
}

dat <- do.call(rbind, ls_beta)
dat$cl <- as.factor(dat$cl)
dat$x <- as.factor(dat$x)
dat$b <- as.factor(dat$b)
dat$e <- as.factor(dat$e)
dat$n <- as.factor(dat$n)

dat %>% dplyr::filter(key == "GenMerged") %>%
    #dplyr::filter(tn == "cv") %>%
    dplyr::filter(cl == 3) %>%
    ggplot( aes( y = value, x = b, fill = n )) +
    facet_wrap( ~ tn, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    )

dat %>% dplyr::filter(key == "Specialist") %>%
    ggplot( aes( y = value, x = b, fill = n )) +
    facet_wrap( ~ cl, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.7, 1.25) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

dat %>% dplyr::filter(key == "ZeroOut") %>%
    ggplot( aes( y = value, x = b, fill = n )) +
    facet_wrap( ~ cl, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.7, 1.25) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",


dat %>% dplyr::filter(key == "MS") %>%
    ggplot( aes( y = value, x = b, fill = n )) +
    facet_wrap( ~ cl, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    #ylim(0, 2) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",



dat %>% group_by(key, x, b, n, cl, tn) %>% 
    summarize(my_mean = mean(value) ) %>% print(n = Inf)
# mutate is like cbind
# spread is like inverse of gather
##########################################################################

# 
# plt2 = dat %>% tibble %>% 
#     #dplyr::filter(n == 6) %>%
#     dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "cvCF") %>%
#     ggplot(aes( y = value, x = b, fill = x )) +
#     facet_wrap( ~ key, nrow = 1) +
#     geom_boxplot(
#         lwd = 1.5, 
#         fatten = 0.5, 
#         alpha = 0.5 
#     ) + 
#     geom_hline(yintercept=1, 
#                linetype="dashed", 
#                color = "black", 
#                size = rel(0.5),
#                alpha = 0.7) + #
#     #ylim(0, 2) +
#     ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$') )+ 
#     xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
#     scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
#     scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     )
# 
# 
# setwd("~/Desktop/Research Final/OEC/OEC Ridge Plots/OEC Sims")
# ggsave( "sims22_specialist_clusts_oec.png",
#         plot = plt,
#         width = 15,
#         height = 10 
# )

######################################################
# compare different K and betaVar -- no regularization
######################################################
setwd("~/Desktop/Research/sims22")
library(dplyr)
library(latex2exp)
bVar <- c( 0.0001, 0.001, 0.005, 0.01, 0.05, 0.25, 0.5, 1, 1.5) # beta variances
xVar <- round( sqrt( c( 0, 0.5, 1, 1.5 ) ), 2)  # beta variances
tune <- c("cvCF") # , "zero", "cv"
etaSt <- c(TRUE)
sclX <- c(TRUE)
clusts <- c(3,6) # no clusters
e = TRUE
sc = TRUE
nVec <- c( 300) # 50,150,
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(cl in clusts){
    for(tn in tune){
        for(n in nVec){
            for(bl in 1:length(bVar)){
                for(xv in 1:length(xVar)){
                    cnt <- cnt + 1 # counter
                    b <- bVar[bl]
                    x <- xVar[xv]
                    flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                    
                    if(file.exists(flNm)){
                        # check to see if file exists
                        d <- read.csv(flNm) 
                        rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                              Specialist = d$oec_country2 / d$stacking_country,
                                              ZeroOut = d$oec_country0 / d$stacking_country_zeroOut,
                                              Merged = d$oec_test2 / d$merge2,
                                              MS = d$stack2 / d$merge2
                        )
                        ls_beta[[cnt]] <- cbind( 
                            gather(rmseMat), 
                            x, 
                            cl,
                            b,
                            e,
                            sc,
                            n,
                            tn
                        )
                        
                    }
                    
                    
                    
                }
            }
        }
    }
}

dat <- do.call(rbind, ls_beta)
dat$cl <- as.factor(dat$cl)
dat$x <- as.factor(dat$x)
dat$b <- as.factor(dat$b)
dat$e <- as.factor(dat$e)
dat$n <- as.factor(dat$n)

dat %>% group_by(key, x, b, n, cl, tn) %>% summarize(my_mean = mean(value) ) %>% print(n = Inf)

# Generalist / Merged

ggplot(dat[dat$key == "Merged",], aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ cl, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
    geom_hline(yintercept=1, 
               linetype="dashed", 
               color = "black", 
               size = rel(0.5),
               alpha = 0.7) + #
    ylim(0.5, 1.5) +
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00"))  # "#525252",

dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(cl == 3) %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ key, nrow = 1) +
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
    ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Stack}}$') )+ 
    xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
    scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
    scale_fill_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) +
    theme_classic(base_size = 12) +
    theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
           axis.text=element_text(face="bold",color="black", size=rel(1.75)),
           axis.title = element_text(face="bold", color="black", size=rel(1.5)),
           legend.key.size = unit(2, "line"), # added in to increase size
           legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
           legend.title = element_text(face="bold", color="black", size = rel(2))
    )

dat %>% tibble %>% dplyr::filter(key == "Generalist") %>%
    ggplot(aes( y = value, x = b, fill = x )) +
    facet_wrap( ~ cl, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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

#############################################################
# tuning style
setwd("~/Desktop/Research/sims22")
library(latex2exp)
bVar <- c(0.001, 0.05, 0.25) # beta variances
xVar <- c(0.71) # beta variances
tune <- c("zero", "cvCF", "cv")
etaSt <- c(TRUE)
sclX <- c(TRUE)
cl <- 3 # no clusters
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

itrs <- 100 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = length(bVar), nr = itrs)
for(sc in sclX){
    
    for(tn in tune){
        for(bl in 1:length(bVar)){
            for(xv in 1:length(xVar)){
                cnt <- cnt + 1 # counter
                b <- bVar[bl]
                x <- xVar[xv]
                flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn, "_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                
                if(file.exists(flNm)){
                    # check to see if file exists
                    d <- read.csv(flNm) 
                    rmseMat <- data.frame(Generalist = d$oec_test2 / d$stack2,
                                          Specialist = d$oec_country2 / d$stacking_country,
                                          ZeroOut = d$oec_country0 / d$stacking_country_zeroOut,
                                          GenMerged = d$oec_test2 / d$merge2
                    )
                    ls_beta[[cnt]] <- cbind( 
                        gather(rmseMat), 
                        x, 
                        cl,
                        b,
                        e,
                        tn
                    )
                    
                }
                
                
                
            }
        }
    }
}

dat <- do.call(rbind, ls_beta)
dat$cl <- as.factor(dat$cl)
dat$x <- as.factor(dat$x)
dat$b <- as.factor(dat$b)
dat$e <- as.factor(dat$e)
dat$tn <- as.factor(dat$tn)

#dat2 <- dat[dat$x == 0.71,]
dat2 = dat

# Generalist / Merged
ggplot(dat2[dat2$key == "GenMerged",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ tn, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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

ggplot(dat2[dat2$key == "Generalist",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ tn, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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


ggplot(dat2[dat2$key == "Specialist",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ tn, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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


ggplot(dat2[dat2$key == "ZeroOut",], aes( y = value, x = b, fill = sc )) +
    facet_wrap( ~ tn, nrow = 1) +
    geom_boxplot(
        lwd = 1.5, 
        fatten = 0.5, 
        alpha = 0.5 
    ) + 
    theme_classic(base_size = 12) +
    #facet_grid( .~  key) +
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


##########################################################################
a <- read.csv("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_TRUE_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_0.05_xVar_0.71_clst_3_reTn_FALSE_K_5_s_513")
boxplot(a$oec_test2 / a$stack2)
boxplot(a$oec_test2 / a$merge2)

mean(a$oec_test2 / a$stack2)
mean(a$oec_test2 / a$merge2)

# works
a <- read.csv("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_TRUE_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_0.05_xVar_0.71_clst_3_reTn_FALSE_K_5_s_513")
boxplot(a$oec_test2 / a$stack2)
boxplot(a$oec_test2 / a$merge2)

mean(a$oec_test2 / a$stack2)
mean(a$oec_test2 / a$merge2)

a <- read.csv("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_TRUE_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_FALSE_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_0.05_xVar_0.71_clst_3_reTn_FALSE_K_5_s_513")
boxplot(a$oec_test2 / a$stack2)
boxplot(a$oec_test2 / a$merge2)

mean(a$oec_test2 / a$stack2)
mean(a$oec_test2 / a$merge2)

a <- read.csv("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_TRUE_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_0.25_xVar_0_clst_3_reTn_FALSE_K_5_s_513")
colMeans(a)

a <- read.csv("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_cvCF_stTnIn_TRUE_stCV_cv_etaSt_TRUE_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_cvSpec_TScv_FALSE_sclX_TRUE_glm_FALSE.pca.FALSE_p20_n_150.150_eps_1.2_smpSzW_1_bVar_0.25_xVar_0.71_clst_3_reTn_FALSE_K_5_s_513")

boxplot(a$oec_test2 / a$stack2)
boxplot(a$oec_test2 / a$merge2)

mean(a$oec_test2 / a$stack2)
mean(a$oec_test2 / a$merge2)



