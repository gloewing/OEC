# Figures and Tables for sims23 
######################################################
# compare different K and betaVar -- no regularization
######################################################
# setwd("~/Desktop/Research/sims23")
setwd("~/Desktop/OEC/Manuscript Figures/General Simulations/sims23")
library(dplyr)
library(latex2exp)
library(kableExtra)
bVar <- c( 0.0001, 0.01, 1) # beta variances #  0.001, 0.005, 0.01,
xVar <- c(0, round( sqrt( c(0.01, 0.5, 1.5 ) ), 2) )  # beta variances 
tune <- c("cvCF", "zero") # , "zero", "cv"
etaSt <- c(TRUE)
sclX <- c(TRUE)
clusts <- c(3,6) # no clusters
e = TRUE
sc = TRUE
nVec <- c(50, 150, 300) # 50,150,
#Kvec <- c(2,4,6,8,12)
ls_beta <- vector(length = length(bVar), "list")

# save average results for table
tableMat <- matrix( nc = 10, nr = length(bVar) * length(xVar) * length(clusts) * length(nVec) * length(tune))
tableMat2 <- matrix( nc = 11, nr = nrow(tableMat))
colnames(tableMat) <- c("Tune", "n", "Clusters", "XVar", "BVar", 
                                                "Generalist OEC", "Merged", 
                                                "Specialist OEC", "Study-Specific",
                                                "Zero Out OEC")
colnames(tableMat2) <- c("Tune", "n", "Clusters", "XVar", "BVar", 
                         "Generalist OEC", "Generalist", 
                         "Specialist OEC", "Specialist",
                         "Zero Out OEC",
                         "Zero Out")

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
                    
                    flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", x, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                    
                    if(x == 0)                    flNm <- paste0("oecSpc2.0.5.etTn_TRUE_cv_stTn_cv_oecTn_",tn,"_stTnIn_TRUE_stCV_cv_etaSt_",e,"_psi_0.Inf_etaSpc_TRUE_Wcv_cv_Wspc_FALSE_TScv_FALSE_sclX_", sc, "_glm_FALSE.pca.FALSE_p20_n_",n,".150_eps_1.2_smpSzW_1_bVar_",b,"_xVar_", 0, "_clst_",cl, "_reTn_FALSE_K_5_s_513")
                    
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
                        
                        # for table -- each relative to their stacking counterpart
                        tableMat[cnt,] <- c(tn, n, cl, x, b, 
                                            colMeans(
                                                cbind(d$oec_test2 / d$stack2,
                                                    d$merge2 / d$stack2,
                                                    d$oec_country2 / d$stacking_country,
                                                    d$country / d$stacking_country,
                                                    d$oec_country0 / d$stacking_country_zeroOut
                                                    )
                                                    )
                                            )
                        
                        # for table - all relative to merging
                        tableMat2[cnt,] <- c(tn, n, cl, x, b, 
                                            colMeans(
                                                cbind(d$oec_test2 / d$merge2,
                                                      d$stack2 / d$merge2,
                                                      d$oec_country2 / d$country,
                                                      d$stacking_country / d$country,
                                                      d$oec_country0 / d$country,
                                                      d$stacking_country_zeroOut / d$country
                                                      )
                                                    )
                                            )
                        # tableMat2[cnt,] <- c(tn, n, cl, x, b, 
                        #                     colMeans(
                        #                             cbind(d$stack2 / d$merge2,
                        #                                     d$oec_test2 / d$merge2,
                        #                                     d$stacking_country / d$merge2,
                        #                                     d$oec_country2 / d$merge2,
                        #                                     d$country / d$merge2,
                        #                                     d$stacking_d$merge2,
                        #                                     d$oec_country0 / d$merge2
                        #                                     )
                        #                             )
                        #                     )
                        
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

dat %>% dplyr::group_by(key, x, b, n, cl, tn) %>% dplyr::summarize(my_mean = mean(value) ) %>% print(n = Inf)
# setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
#########################
# Generalist / Merged
#########################
# clusters - cvCF
plt_cvCF <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap(~key) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.5, 1.2)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("Generalist", "Merged")) %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

# clusters - zero
plt_zero <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap(~key) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.5, 1.2)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("Generalist", "Merged")) %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_generalist_merged_CLust_oec_NoWspecTn_cvCF.png",
        plot   = plt_cvCF,
        width  = 3,
        height = 4)

ggsave( "sims23_generalist_merged_CLust_oec_NoWspecTn_zero.png",
        plot   = plt_zero,
        width  = 3,
        height = 4)

rm(plt_cvCF, plt_zero)

# clusters - plot together all standardized by the merged tested on the specialist country - cvCF
plt_cvCF <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key, color = key )) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                      labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                       labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = "Method") +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "cvCF") %>%
#     dplyr::filter(x == 0.01) %>%
#     ggplot(aes( y = value, x = b, fill = key )) +
#     #facet_wrap( ~ key, nrow = 1) +
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
#     ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
#     xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
#     scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
#     scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
#                       values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title="Method"))

# clusters - plot together all standardized by the merged tested on the specialist country - zero
plt_zero <- dat %>% tibble %>%
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key, color = key )) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                      labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                       labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = "Method") +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
#     dplyr::filter(x == 0.01) %>%
#     ggplot(aes( y = value, x = b, fill = key )) +
#     #facet_wrap( ~ key, nrow = 1) +
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
#     ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
#     xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
#     scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
#     scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
#                       values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title="Method"))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_specialist_clusts_oecTogether_NoWspecTn_cvCF.png",
        plot   = plt_cvCF,
        width  = 5,
        height = 3)

ggsave( "sims23_specialist_clusts_oecTogether_NoWspecTn_zero.png",
        plot   = plt_zero,
        width  = 5,
        height = 3)

rm(plt_cvCF, plt_zero)

# no clusters - cvCF
plt1_cvCF <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap(~key) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.5, 1.2)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("Generalist", "Merged")) %>%
#     dplyr::filter(cl == 6) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "cvCF") %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
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
#     scale_fill_manual(values = c( "red", "blue", "green", "#0868ac", "#E69F00")) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

# no clusters - zero
plt1_zero <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("Generalist", "Merged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap(~key) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.5, 1.2)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("Generalist", "Merged")) %>%
#     dplyr::filter(cl == 6) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
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
#     scale_fill_manual(values = c( "red", "blue", "green", "#0868ac", "#E69F00")) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_generalist_merged_noCLust_oec_NoWspecTn_cvCF.png",
        plot   = plt1_cvCF,
        width  = 8,
        height = 3)

ggsave( "sims23_generalist_merged_noCLust_oec_NoWspecTn_zero.png",
        plot   = plt1_zero,
        width  = 8,
        height = 3)
rm(plt1_cvCF, plt1_zero)

# clusters - plot together all standardized by the merged tested on the specialist country - cvCF
plt_cvCF <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key, color = key )) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                      labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                       labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = "Method") +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
#     dplyr::filter(cl == 6) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "cvCF") %>%
#     dplyr::filter(x == 0.01) %>%
#     ggplot(aes( y = value, x = b, fill = key )) +
#     #facet_wrap( ~ key, nrow = 1) +
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
#     ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
#     xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
#     scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
#     scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
#                       values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title="Method"))

# clusters - plot together all standardized by the merged tested on the specialist country - zero
plt_zero <- dat %>% tibble %>% 
    dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    dplyr::filter(x == 0.01) %>%
    ggplot(aes( y = value, x = b, fill = key, color = key )) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    scale_fill_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                      labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                       labels = c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = "Method") +
    theme_bw()
# dat %>% tibble %>% 
#     dplyr::filter(key %in% c("ZeroMerged", "SpecMerged", "ZeroOecMerged", "SpecOecMerged")) %>%
#     dplyr::filter(cl == 6) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
#     dplyr::filter(x == 0.01) %>%
#     ggplot(aes( y = value, x = b, fill = key )) +
#     #facet_wrap( ~ key, nrow = 1) +
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
#     ylab(TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Merged}}$') )+ 
#     xlab(TeX('$\\mathbf{\\sigma^2_{\\beta}}}$')) + 
#     scale_color_manual(values = c("red", "blue", "green", "#0868ac", "#E69F00")) + # "#525252",
#     scale_fill_manual(labels=c("Specialist","Specialist-OEC","Zero Out", "Zero Out OEC"),
#                       values = c("red", "blue", "green", "#0868ac", "#E69F00") ) +
#     theme_classic(base_size = 12) +
#     ylim(0.5, 1.2) + 
#     theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2), face="bold"),
#            axis.text=element_text(face="bold",color="black", size=rel(1.75)),
#            axis.title = element_text(face="bold", color="black", size=rel(1.5)),
#            legend.key.size = unit(2, "line"), # added in to increase size
#            legend.text = element_text(face="bold", color="black", size = rel(1.75)), # 3 GCL
#            legend.title = element_text(face="bold", color="black", size = rel(2)),
#            strip.text.x = element_text(face="bold", color="black", size = rel(2))
#     ) + guides(fill=guide_legend(title="Method"))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_specialist_Noclusts_oecTogether_NoWspecTn_cvCF.png",
        plot   = plt_cvCF,
        width  = 5,
        height = 3)

ggsave( "sims23_specialist_Noclusts_oecTogether_NoWspecTn_zero.png",
        plot   = plt_zero,
        width  = 5,
        height = 3)
rm(plt_cvCF, plt_zero)

# cvCF clusters
plt2_cvCF <- dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap( ~ key, nrow = 1) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.50, 1.20)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     #dplyr::filter(n == 6) %>%
#     dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "cvCF") %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

# zero clusters
plt2_zero <- dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 3) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap( ~ key, nrow = 1) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.50, 1.20)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     #dplyr::filter(n == 6) %>%
#     dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
#     dplyr::filter(cl == 3) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
#     dplyr::filter(x %in% c(0.01, 0.5, 1.5) ) %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_specialist_clusts_oec_NoWspecTn_cvCF.png",
        plot   = plt2_cvCF,
        width  = 8,
        height = 3)

ggsave( "sims23_specialist_clusts_oec_NoWspecTn_zero.png",
        plot   = plt2_zero,
        width  = 8,
        height = 3)
rm(plt2_cvCF, plt2_zero)

# clusters cvCF
plt2_cvCF <- dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "cvCF") %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap( ~ key, nrow = 1) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.50, 1.20)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     #dplyr::filter(n == 6) %>%
#     dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
#     dplyr::filter(cl == 6) %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))

# clusters zero
plt2_zero <- dat %>% tibble %>% 
    #dplyr::filter(n == 6) %>%
    dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
    dplyr::filter(cl == 6) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter(tn == "zero") %>%
    ggplot(aes( y = value, x = b, fill = x, color = x )) +
    facet_wrap( ~ key, nrow = 1) +
    geom_hline(yintercept = 1,
               linetype   = 2,
               color      = "darkgray") +
    geom_boxplot(size         = 0.20,
                 outlier.size = 0.50,
                 show.legend = FALSE) +
    geom_boxplot(size          = 0.20,
                 color         = "black",
                 outlier.color = NA) +
    coord_cartesian(ylim = c(0.50, 1.20)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Reds")) +
    labs(x    = TeX('$\\mathbf{\\sigma^2_{\\beta}}}$'),
         y    = TeX('$\\mathbf{RMSE_{OEC}/RMSE_{Method}}$'),
         fill = TeX('$\\mathbf{\\sigma^2_{x}}$')) +
    theme_bw()
# dat %>% tibble %>% 
#     #dplyr::filter(n == 6) %>%
#     dplyr::filter(key %in% c("Specialist", "ZeroOut")) %>%
#     dplyr::filter(cl == 6) %>%
#     dplyr::filter(n == 300) %>%
#     dplyr::filter(tn == "zero") %>%
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
#     ) + guides(fill=guide_legend(title=TeX('$\\mathbf{\\sigma^2_{x}}$')))


setwd("~/Desktop/Research Final/Mortality/Figures/Final Figures/General Simulations/sims23")
ggsave( "sims23_specialist_Noclusts_oec_NoWspecTn_cvCF.png",
        plot   = plt2_cvCF,
        width  = 8,
        height = 3)

ggsave( "sims23_specialist_Noclusts_oec_NoWspecTn_zero.png",
        plot   = plt2_zero,
        width  = 8,
        height = 3)
rm(plt2_cvCF, plt2_zero)

###############################################################################################
######################################################
# Tables - Stacking comparison (tableMat) - compared to relevant stacking version
######################################################
xV <- c(0, round( sqrt( c(0.01, 0.5, 1.5 ) ), 2) )  # x variances 
tm <- tableMat
tuned <- "zero"

############
# full table
############
tm <- as.data.frame( tableMat[ rowSums(is.na(tableMat)) < ncol(tableMat), ] ) %>% # remove rows with all NAs   -c(1,2,3)
    tibble %>%
    dplyr::filter(Tune %in% tuned ) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter( Clusters %in%  clusts ) %>%
    dplyr::filter(XVar %in% xV )

tm <- tm[ , -c(1,2)] %>%
    mutate_all(function(x) as.numeric(as.character(x))) # convert to numeric

tm$XVar <- round(tm$XVar^2, 2) # turn from sd into variance
tm$XVar <- ifelse(tm$XVar == 1.49, 1.50, tm$XVar) # just to get the re-rounding right

##### full table #####
print(paste0("Full table Tuning: ", tuned))
# table for generalist
tm %>% 
    round(3) %>%
    kable(format = "latex", booktabs = T) %>% 
    kable_styling(position = "center") %>% 
    print

############
# main table
############
tm <- tableMat
tuned <- "zero"

# includes a subset for space
xV <- round( sqrt( c(0.01, 1.5 ) ), 2)
bV <- round(c(1e-4, 1.00), 4)
tm <- as.data.frame( tableMat[ rowSums(is.na(tableMat)) < ncol(tableMat), ] ) %>% # remove rows with all NAs   -c(1,2,3)
    tibble %>%
    dplyr::filter(Tune %in% tuned ) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter( Clusters %in%  clusts ) %>%
    dplyr::filter(XVar %in% xV ) %>%
    dplyr::filter(BVar %in% bV )

tm <- tm[ , -c(1,2)] %>%
    mutate_all(function(x) as.numeric(as.character(x))) # convert to numeric

tm$XVar <- round(tm$XVar^2, 2) # turn from sd into variance
tm$XVar <- ifelse(tm$XVar == 1.49, 1.50, tm$XVar) # just to get the re-rounding right

##### full table #####
print(paste0("Full table Tuning: ", tuned))
# table for generalist
tm %>% 
    round(3) %>%
    kable(format = "latex", booktabs = T) %>% 
    kable_styling(position = "center") %>% 
    print

#######################################################################
# Tables - All relative to study-specific model or merged (tableMat2)
#######################################################################
xV <- c(0, round( sqrt( c(0.01, 0.5, 1.5 ) ), 2) )  # x variances 
tm <- tableMat2
tuned <- "zero"

############
# full table
############
tm <- as.data.frame( tableMat2[ rowSums(is.na(tableMat2)) < ncol(tableMat2), ] ) %>% # remove rows with all NAs   -c(1,2,3)
    tibble %>%
    dplyr::filter(Tune %in% tuned ) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter( Clusters %in%  clusts ) %>%
    dplyr::filter(XVar %in% xV )

tm <- tm[ , -c(1,2)] %>%
    mutate_all(function(x) as.numeric(as.character(x))) # convert to numeric

tm$XVar <- round(tm$XVar^2, 2) # turn from sd into variance
tm$XVar <- ifelse(tm$XVar == 1.49, 1.50, tm$XVar) # just to get the re-rounding right

##### full table #####
print(paste0("Full table Tuning: ", tuned))
# table for generalist
tm %>% 
    round(3) %>%
    kable(format = "latex", booktabs = T) %>% 
    kable_styling(position = "center") %>% 
    print

############
# main table
############
tm <- tableMat2
tuned <- "zero"

# includes a subset for space
xV <- round( sqrt( c(0.01, 1.5 ) ), 2)
bV <- round(c(1e-4, 0.01, 1.00), 3)
tm <- as.data.frame( tableMat2[ rowSums(is.na(tableMat2)) < ncol(tableMat2), ] ) %>% # remove rows with all NAs   -c(1,2,3)
    tibble %>%
    dplyr::filter(Tune %in% tuned ) %>%
    dplyr::filter(n == 300) %>%
    dplyr::filter( Clusters %in%  clusts ) %>%
    dplyr::filter(XVar %in% xV ) %>%
    dplyr::filter(BVar %in% bV )

tm <- tm[ , -c(1,2)] %>%
    mutate_all(function(x) as.numeric(as.character(x))) # convert to numeric

tm$XVar <- round(tm$XVar^2, 2) # turn from sd into variance
tm$XVar <- ifelse(tm$XVar == 1.49, 1.50, tm$XVar) # just to get the re-rounding right

##### full table #####
print(paste0("Main table Tuning: ", tuned))
# table for generalist
tm %>% 
    round(3) %>%
    kable(format = "latex", booktabs = T) %>% 
    kable_styling(position = "center") %>% 
    print


