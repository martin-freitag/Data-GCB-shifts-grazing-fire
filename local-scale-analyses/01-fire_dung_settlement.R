
rm(list=ls())


library(glmmTMB)
library(MASS) #for negbinomial family in plain glm
library(MuMIn)
library(performance)
library(DHARMa)
library(segmented)
library(ggplot2)
library(ggeffects)
library(ggrepel)
library(ggpubr) # to ggarrange plots

# Retrieving Raw Data and Cleaning ----------------------------------------

envST <- read.table("02_envST.csv", sep=",", header=T, row.names=1)
envST$Dung <- log(envST$dung_livestock+1)




### model grazing intensity as a function of distance to nearest settlement -----------

dungmod <- glm(dung_livestock ~ dist_settl_km, data = envST, family = negative.binomial(theta=2, link = "log")) #negbinomial family from MASS package
summary(dungmod)
r2(dungmod)

res <- simulateResiduals(dungmod)
plot(res)

segdung <- segmented(dungmod, seg.Z= ~dist_settl_km)
summary(segdung)
plot(segdung)

# segmented model is no better (AIC) than the plain GLM, so no threshold


# predict regression line 

dung_pred <- ggpredict(dungmod, "dist_settl_km")

newdata <- data.frame(dist_settl_km = seq(0,30, by=1))

pred_dung <- predict(dungmod, 
                        newdata = newdata,
                        type = "response", se.fit = TRUE
)
str(dung_pred)  

dung_pred <- rbind(dung_pred      , dung_pred, dung_pred, dung_pred, dung_pred[c(1:3),] )

dung_pred$x         <- newdata$dist_settl_km
dung_pred$predicted <- pred_dung$fit  
dung_pred$std.error <- pred_dung$se.fit  
dung_pred$conf.low  <- pred_dung$fit - 1.96*pred_dung$se.fit  
dung_pred$conf.high <- pred_dung$fit + 1.96*pred_dung$se.fit  

plot(dung_pred)





### model fire frequency as a function of grazing intensity -----------

firemod <- glm(fire_freq_check ~ dung_livestock, data = envST, family = poisson)
summary(firemod)
r2(firemod)

res <- simulateResiduals(firemod)
plot(res)
plot(residuals(res) ~ envST$dung_livestock)

# predict simple model only to have a template, which will be filled with the predictions from the segmented model later
fire_pred <- ggpredict(firemod, "dung_livestock")


# check if there is a threshold with the segmented package, i.e. is the model fit better (AIC) with a broken-stick regression
segfire <- segmented(firemod, seg.Z= ~ dung_livestock, psi = c(5))
summary(segfire)
r2(segfire) 
# the segmented fire model has an AIC lower than without the segments, delta AIC is 1.97 - so, there is a threshold
# no R2 available for this kind of model - crude correlation
cor(envST$fire_freq_check, predict(segfire, type="link"))^2



res <- simulateResiduals(segfire)
plot(res, asFactor = T)
plotResiduals(res, form=envST$dung_livestock)
# model fit is ok


# predict regression line for plotting
newdata <- data.frame(dung_livestock = seq(0,40, by=1))

pred_segfire <- predict(segfire, 
                        newdata = newdata,
                        type = "response", se.fit = TRUE
                        )
str(pred_segfire)  

fire_pred <- rbind(fire_pred      , fire_pred[-1, ] )

fire_pred$x         <- newdata$dung_livestock
fire_pred$predicted <- pred_segfire$fit  
fire_pred$std.error <- pred_segfire$se.fit  
fire_pred$conf.low  <- pred_segfire$fit - 1.96*pred_segfire$se.fit  
fire_pred$conf.high <- pred_segfire$fit + 1.96*pred_segfire$se.fit  

plot(fire_pred)




### save the residuals of fire frequency, after correcting for grazing intensity -----------------
resid_fire <- residuals(segfire, type="response") #before I used the 'deviance' residuals, now on response scale - correlation between both is 0.96



### plot dung and fire regressions ----------------------

my_theme <-   theme_bw() + theme(panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 plot.title = element_text(size=11),
                                 axis.title.y=element_text(colour="black", size=9) ,
                                 axis.title.x=element_text(colour="black", size=9) ,
                                 axis.ticks.x=element_line(colour="black", size=0.5) ,
                                 axis.text=element_text(size=8, colour="black"),
                                 panel.border=element_rect(colour="black", fill=NA) ) 

p_fire_dung <- plot(fire_pred)      +  
  geom_jitter(data = envST, aes(x=dung_livestock, y=fire_freq_check), alpha = 0.3, shape = 16, size = 1,
              width=0, height=0.1) +
  labs(title = NULL, x = "No. dung piles (per 200m²)", y = "Fire frequency (2000-2015)") + 
  annotate("text", x=30, y=4.5, label= "Pseudo R²=0.60\nN=204", size = 2.75) +
  my_theme

p_fire_dung


p_dung_dist <- plot(dung_pred)      +  
  geom_point(data = envST, aes(x=dist_settl_km, y=dung_livestock), alpha = 0.3, shape = 16, size = 1) +
  labs(title = NULL, x = "Distance to nearest settlement (km)", y = "No. dung piles (per 200m²)") + 
  annotate("text", x=25, y=35, label= "Pseudo R²=0.87\nN=204", size = 2.75) +
  my_theme

p_dung_dist

ggarrange(p_dung_dist    , 
          p_fire_dung    , 
          labels = c("(a)", "(b)"),
          font.label = list(size = 12),
          hjust = 0.05,
          ncol = 2, nrow = 1
) 




ggsave("01-FireFreq-dung-segmented.pdf", plot=last_plot(), device="pdf",
       width=16, height=7, units="cm", dpi=1200)

ggsave("01-FireFreq-dung-segmented.png", plot=last_plot(), device="png",
       width=16, height=7, units="cm", dpi=1200)


