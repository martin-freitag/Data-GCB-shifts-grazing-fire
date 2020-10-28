



### get the spatial fire and abandonment data ----------------------

spatial_all <- read.csv("04_Settlements_4km_GrazFire.csv", header=T, sep=",")

spatial <- subset(spatial_all,  Grassl_M >= 0.5)           # restrict to settlements with at least 50% grassland in the 4 km buffer around settlemtns
spatial$Status2015[which(spatial$Status2015 > 100)] <- 100 # change growing settlements to 100% intactness
spatial$Status2015[which(spatial$Status2015 == 6)] <- 10   # round to 10% steps 
spatial$Status2015[which(spatial$Status2015 == 1)] <- 0    # round to zero



spatial$BA_increase <- spatial$BA1990_M - spatial$BA2015_M
spatial$BA_increase <- spatial$BA_increase * 100

spatial$BA1990_M    <- spatial$BA1990_M * 100 # change to percent burned area in 4km buffer
spatial$BA2015_M    <- spatial$BA2015_M * 100


spatial1990 <- spatial                        # reorganize data for plotting 1990 and 2015 burned area simultaneously
spatial1990$Period <- rep("1990", nrow(spatial1990))
spatial1990$BA_proportion <- spatial$BA1990_M

spatial2015 <- spatial
spatial2015$Period <- rep("2015", nrow(spatial2015))
spatial2015$BA_proportion <- spatial$BA2015_M

spatial_both <- rbind(spatial1990, spatial2015)
spatial_both$Period <- factor(spatial_both$Period)





# Model the 90% quantiles using splines with smooth dimension k=4. 
# k set to this value by subject choice, as it is a trade-off between smoothness (high k) and overall patterns
# Bayesian framework used as package 'quantreg' was not flexible enough and models partly did not converge
# default priors

library(brms) #version 2.13


# smooth k
k <- 4

# no of iterations
iter <- 10000

MODIS_90 <- brm(bf(FireFr_M ~ s(Status2015, k=k), quantile = 0.90), data = spatial,
              family = asym_laplace(), control = list(adapt_delta = 0.9, max_treedepth = 15), 
              iter = iter, chains = 4, cores=4)

BA1990_90 <- brm(bf(BA1990_M ~ s(Status2015, k=k), quantile = 0.90), data = spatial,
                 family = asym_laplace(), control = list(adapt_delta = 0.9, max_treedepth = 15), 
                 iter = iter, chains = 4, cores=4)

BA2015_90 <- brm(bf(BA2015_M ~ s(Status2015, k=k), quantile = 0.90), data = spatial,
                family = asym_laplace(), control = list(adapt_delta = 0.9, max_treedepth = 15), 
                iter = iter, chains = 4, cores=4)


# generate newdata to create predicted regression lines of the mean +- credible interval

newdata_MODIS_90 <- data.frame(Status2015 = seq(0, 100, by=1) )
fitted_MODIS_90 <- fitted(MODIS_90, newdata = newdata_MODIS_90, dpar="mu")
# rename the fitted value to name of the original variable and generate fitted values
colnames(fitted_MODIS_90)[1] <- "FireFr_M"
fitted_MODIS_90 <- cbind(fitted_MODIS_90, newdata_MODIS_90)


newdata_BA1990_90 <- data.frame(Status2015 = seq(0, 100, by=1) )
fitted_BA1990_90 <- fitted(BA1990_90, newdata = newdata_BA1990_90, dpar="mu")
# rename the fitted value to name of the original variable and generate fitted values
colnames(fitted_BA1990_90)[1] <- "BA_proportion"
fitted_BA1990_90 <- cbind(fitted_BA1990_90, newdata_BA1990_90, Period = rep("1990", nrow(newdata_BA1990_90)))


newdata_BA2015_90 <- data.frame(Status2015 = seq(0, 100, by=1) )
fitted_BA2015_90 <- fitted(BA2015_90, newdata = newdata_BA2015_90, dpar="mu")
# rename the fitted value to name of the original variable and generate fitted values
colnames(fitted_BA2015_90)[1] <- "BA_proportion"
fitted_BA2015_90 <- cbind(fitted_BA2015_90, newdata_BA2015_90, Period = rep("2015", nrow(newdata_BA2015_90)))


# arrange predicted data of 90% quantiles into one data frame  
fitted_BA_90 <- rbind(fitted_BA1990_90, fitted_BA2015_90)






# plot data and 90% quantiles

library(ggplot2)
library(ggpubr) # to ggarrange plots


# some aesthetics
my_theme <-   theme_bw() + theme(panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 plot.title = element_text(size=8),
                                 axis.title.y=element_text(colour="black", size=8) ,
                                 axis.title.x=element_text(colour="black", size=8) ,
                                 axis.ticks.x=element_line(colour="black", size=0.5) ,
                                 axis.text=element_text(size=8, colour="black"),
                                 panel.border=element_rect(colour="black", fill=NA),
                                 legend.title=element_text(size=8),
                                 legend.text=element_text(size=8) ) 



# manual legend items (lines) for MODIS fire frequency plot
freq_line <- data.frame(x = c(56, 68), y = c(3.5, 3.5))


p_freq <- ggplot(data=fitted_MODIS_90, aes(x=Status2015, FireFr_M)) +
  geom_boxplot(data = spatial, aes(group=Status2015), 
               fill="gray90", outlier.color="gray50", outlier.size=0.15, lwd=0.25,
               width=5) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.5, show.legend = FALSE ) +
  geom_line(size=0.5, linetype="dotted", show.legend = FALSE ) +
  geom_line(data = freq_line, aes(x=x, y=y, group=1), 
            size=0.5, linetype="dotted") +
  annotate(geom="text", x=88, y=3.5, label="90% quantile", size=2.75) +
  scale_x_continuous(breaks = seq(0, 100, by=10), labels = c(0, "", "20", "", "40", "", "60", "", "80", "", "100")) +
  labs(title = NULL, x = "Livestock infrastructure intactness ca. 2013", y = "Mean fire frequency 2000-2015 in 4 km radius") + 
  my_theme
p_freq




# manual legend items (lines) for burned area plot
area_line <- data.frame(x = c(70, 78), y = c(60, 60))


p_area <- ggplot(data=fitted_BA_90, aes(x=Status2015, y=BA_proportion  )) +
  geom_boxplot(data = spatial_both, aes(group=interaction(Status2015, Period), colour = Period ), 
               fill="gray90", outlier.size=0.15, lwd=0.25,
               position=position_dodge(7.5), width=6.5) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill=Period), alpha = 0.5,
              position=position_dodge(7.5), show.legend = FALSE ) +
  geom_line( aes(colour=Period), size=0.5, linetype="dotted", 
             position=position_dodge(7.5), show.legend = FALSE ) +
  scale_x_continuous(breaks = seq(0, 100, by=10), labels = c(0, "", "20", "", "40", "", "60", "", "80", "", "100")) +
  geom_line(data = area_line, aes(x=x, y=y, group=1), 
            size=0.5, linetype="dotted") +
  annotate(geom="text", x=92.5, y=60, label="90% quantiles", size = 2.75) +
  labs(title = NULL, x = "Livestock infrastructure intactness ca. 2013", y = "Area burned in 4 km radius (%)") + 
  my_theme +
  theme (legend.position=c(.75, .8))
p_area



ggarrange(p_freq, p_area, labels = c("a)", "b)"), 
          ncol=2, nrow=1,
          font.label = list(size = 16),
          vjust=1.1, hjust=-0.1,
          widths = c(1, 1.5)
)

ggsave("04-fire_spatial_k4.pdf", plot=last_plot(), device="pdf",
       width=16, height=8, units="cm", dpi=1200)

ggsave("04-fire_spatial_k4.png", plot=last_plot(), device="png",
       width=16, height=8, units="cm", dpi=1200)







### histogram of fire frequency and burned area ---------


max(spatial$FireFr_M)
min(spatial$FireFr_M)


hist_freq <- ggplot(data=spatial, aes(x=FireFr_M)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, 
                 color="gray40", fill="gray80", binwidth=0.1)+
  labs(title = NULL, x = "Mean fire frequency 2000-2015 in 4 km radius", y = "Density") + 
  my_theme
hist_freq


max(spatial_both$BA_proportion)

hist_ba <- ggplot(data=spatial_both, aes(x=BA_proportion, color=Period, fill=Period )) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, 
                 binwidth=2)+
  scale_x_continuous(breaks = seq(0, 100, by=10), labels = c(0, "", "20", "", "40", "", "60", "", "80", "", "100")) +
  labs(title = NULL, x = "Area burned in 4 km radius (%)", y = "Density") + 
  my_theme +
  theme (legend.position=c(.75, .8))
hist_ba



spatial_both$log_BA_proportion <- spatial_both$BA_proportion + 1 


hist_ba_log <- ggplot(data=spatial_both, aes(x=log_BA_proportion, color=Period, fill=Period )) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, 
                 binwidth=0.025)+
  scale_x_log10(breaks = c(1, 11, 101), labels = c("0", "10", "100") ) +
  labs(title = NULL, x = "Area burned in 4 km radius (%)", y = "Density") + 
  my_theme +
  theme (legend.position=c(.75, .8))
hist_ba_log



ggarrange(hist_freq, hist_ba_log, 
          labels = c("a)", "b)"), 
          ncol=2, nrow=1,
          font.label = list(size = 16),
          vjust=1.1, hjust=-0.1
)

ggsave("Figures/04-fire_spatial_hist.pdf", plot=last_plot(), device="pdf",
       width=16, height=8, units="cm", dpi=1200)

ggsave("Figures/04-fire_spatial_hist.png", plot=last_plot(), device="png",
       width=16, height=8, units="cm", dpi=1200)






# ### segmented regression with brms --------------
# 
# # this Bayesian version of segmented regression was used at some point to test 
# # whether there is evidence for break-points in the response of fire freq / burned area
# # to settlement abandonment status. This analysis was not used in the end, but I'll keep it here as 
# # I might re-use this model code
# 
# 
# library(brms)
# 
# # The model
# bform <- bf(
#   BA_proportion_2015 ~ Intercept + slope1 * Status2015 * step(change - Status2015) +  # Section 1
#                        (slope1 * change + slope2 * (Status2015 - change)) * step(Status2015 - change),  # Section 2
#   Intercept + slope1 + slope2 + alpha ~ 1,  # Fixed intercept and slopes
#   #change ~ 1 ,  # no variation in change-points per groups or else
#   nl = TRUE,
#   nlf(change ~ inv_logit(alpha) * 10)
#   #quantile = 0.5 # specify quantile to model
# )
# 
# # Priors
# bprior <- prior(normal(0, 100), nlpar = "Intercept") +
#   prior(normal(0, 10), nlpar = "slope1") +
#   prior(normal(0, 10), nlpar = "slope2") +
#   prior(normal(0, 10), nlpar = "alpha")  # Within observed range of Status2015: [0,100]
# 
# # # Initial values
# # inits = list(list(
# #  slope1 = 0.5,
# #  slope2 = 0.5,
# #  Intercept = 0.5
# # ))
# 
# # Fit it!
# fit3 <- brm(bform, prior = bprior, chains = 4, # inits = inits,  # maybe no need for initial values
#             iter=10000,
#             family=Beta(link = "logit", link_phi = "log"),
#             #family = asym_laplace(link="inverse"),  # asym laplace family needed for quantile regression
#             data = spatial)
# 
# summary(fit3)
# conditional_effects(fit3)


