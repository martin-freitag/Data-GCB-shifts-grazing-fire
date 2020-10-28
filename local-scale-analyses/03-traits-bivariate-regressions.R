


# updating to R v. 4.0.1 changed the default data.frame behaviour. it sets stringsAsFactors = FALSE and breaks my form of subsetting
# therfore; i need to change from indexing to select()


library(vegan)
library(FD)
library(glmmTMB)
library(performance)
library(DHARMa)
library(ggplot2)
library(ggeffects)
library(ggpubr) # to ggarrange plots

# Retrieving Raw Data and Cleaning ----------------------------------------

envST <- read.table("02_envST.csv", sep=",", header=T, row.names=1)
speST <- read.table("01_speST.csv", sep=",", header=T, row.names=1)
traST <- read.table("03_traST_cat.csv", sep=",", header=T, row.names=1, na.strings= c("", "NA"))

# add the residuals of the firefrequency ~ dungpiles regression to the data!
#source: r script 'revision_fire_dung_settlement.R'
envST$resid_fire <- resid_fire
  
  

colnames(speST) <- make.cepnames(colnames(speST)) # abbreviate species names
rownames(traST) <- make.cepnames(rownames(traST))


##prepare data frame containing constraining variables
cons <- envST[, -c(1 : ncol(envST))]
cons$FireFreq      <- log(envST$fire_freq_check+1)
cons$resid_fire    <- envST$resid_fire
#cons$TimeSinceFire <- log(envST$timesincefire_check+1) 
cons$Dung          <- log(envST$dung_livestock + 1)   # think about logging Dung
cons$eC            <- log( (envST$LF_530 + envST$LF_3060)/2) #mean eC. think about log-transforming
#cons$pH            <- envST$pH_530
#cons$Sand          <- envST$Sand_530
#cons$BM_yield_m2   <- envST$BM_yield_m2


##check, if trait information of at least 80% cover is available. SLA is available for the same species as LDMC is
covLDMC <- speST[, -which(is.na(traST$LDMC))]
rel_LDMC <- (rowSums(covLDMC) / rowSums(speST))*100
del_LDMC <- rel_LDMC[which(rel_LDMC < 80)]

##remove plots with traits for less than 80% cover for LDMC, SLA and Height
species <- speST[ - which(row.names(speST) %in% names(del_LDMC)) ,]
cons <- cons[ - which(row.names(cons) %in% names(del_LDMC)) ,]
rm(covLDMC, rel_LDMC, del_LDMC)

## aim for as many complete cases as possible! Workflow including LDMC and SLA

cons <- cons[ which(complete.cases(cons)) ,] #one plots is missing soil data
species <- species[ row.names(species)%in%row.names(cons) ,]
selection<-colSums(species) > 0
species <- species[,selection] # take only the species that still occur in the species data frame

rm(selection)

## if needed: remove species that do not have a full trait set
#traits <- traST[complete.cases(traST) , ]
#species <- species[,names(species) %in% row.names(traits)] #keep species which are kept in traST
#traits <- traits[row.names(traits)%in%names(species) , ]

# for now, use the full set with NAs
traits <- traST
traits <- traits[row.names(traits)%in%names(species) , ]
# 92 species with full set, i.e. with SLA and LDMC


### calculate Community weighted means CWMs -----------------




######################################################################
traits <- subset(traits, select = c("LDMC", "SLA", "growthform", "life_history") )
######################################################################




# species data must be matrix for functcomp function

growthforms <- functcomp(subset(traits, select="growthform"), as.matrix(species), CWM.type="all")
colnames(growthforms) <- levels(traits$growthform)
growthforms <- data.frame(apply(growthforms, MARGIN = 2, function(x) (x * rowSums(species) / 100)  ), 
                          stringsAsFactors = TRUE)
# now on cover scale 0 to 1

# linear transformation to replace zeros, because beta regression does not allow for zeros. Recommended by Douma & Weedon 2020 Methods Ecol Evol
#growthforms["dwarfshrub"] <- (growthforms["dwarfshrub"]* (nrow(growthforms)-1)  + (1/ncol(growthforms)) ) / nrow(growthforms)
growthforms[,2] <- (growthforms[,2]* (nrow(growthforms)-1)  +  0.5 ) / nrow(growthforms)
growthforms[,3] <- (growthforms[,3]* (nrow(growthforms)-1)  + 0.5 ) / nrow(growthforms)
growthforms[,4] <- (growthforms[,4]* (nrow(growthforms)-1)  + 0.5 ) / nrow(growthforms)


# append predictors to growthforms DF, scale the predictor variables first
cons_scl <- apply(cons, MARGIN = 2, function(x) scale(x)[,1] )
growthforms <- cbind(growthforms, cons_scl)


# calculate CWMs of leaf traits to fit the same models
CWMleaf <- functcomp(subset(traits, select = c("SLA", "LDMC") ), as.matrix(log(species+1)))

# or with original cover
### CWMleaf <- functcomp(traits[c("SLA", "LDMC")], as.matrix(species))
CWMleaf <- cbind(CWMleaf, cons_scl)


# CWM for height
# we do not use global height values for species, but have measured height ourselves
dimcwm <- read.table("03_PlantDimST.csv", sep=",", head=T)
dimcwm$meanH <- apply(dimcwm[,c("HeightIndividual1", "HeightIndividual2", "HeightIndividual3", "HeightIndividual4")], 
                      MARGIN = 1, function(x) mean(x, na.rm = TRUE))

dimcwm$log_cover <- log(dimcwm$cover + 1)
dimcwm$covxmH <- dimcwm$log_cover * dimcwm$meanH #weight height by cover (log-transformed)
CWMheight <- with(dimcwm, by(covxmH, plot, sum) / by(log_cover, plot, sum))

### # weighted by original cover
### dimcwm$covxmH <- dimcwm$cover * dimcwm$meanH #weight height by original cover
### CWMheight <- with(dimcwm, by(covxmH, plot, sum) / by(cover, plot, sum))

#procedure was tested by hand in Excel, correct
CWMheight <- data.frame(CWMheight = as.vector(CWMheight))
rownames(CWMheight) <- rownames(envST) #unnecessary step, but I'll keep this the way it is...
CWMheight <- CWMheight[ row.names(CWMheight)%in%row.names(cons) , ] #subset
CWMheight <- data.frame(cbind(CWMheight, cons_scl))
max(CWMheight$CWMheight)

# if log species cover used
CWMheight <- CWMheight[ -which(CWMheight$CWMheight > 50) , ]



### start modelling the growthform covers with beta distribution ------------------

# first define a pseudo R squared function for glmmTMB model objects, as the glmmTMB package doesn't provide any R2 (nor does performance package)
pseudoR2 <- function(x, link = NULL) {
  if(class(x) != "glmmTMB") {"No glmmTMB object"}
  if(link %in% c("logit", "identity", "log")) {
    eta <- predict(x, type="link")
    }
  if(link == "logit") {
    y_link <- log(x$frame[,1]/(1-x$frame[,1]))
    }
  if(link == "log") {
    y_link <- log(x$frame[,1])
  }
  if(link == "identity") {
    y_link <- x$frame[,1]
    }
  cor(eta, y_link)^2
}



# model with the glmmTMB package
betagrass <- glmmTMB(grass ~ Dung + resid_fire + eC  , data = growthforms, # pH excluded due to collinearity with soil eC
                     family = beta_family(link = "logit")
                     ) 
summary(betagrass) #no fire: AIC -333.5, no dung: AIC -351.8, both= 350.1
pseudoR2(betagrass, link ="logit")
round(p.adjust(summary(betagrass)$coefficients$cond[,4], method = "fdr"),4)
res <- simulateResiduals(betagrass)
plot(res) 
plotResiduals(res, form = growthforms$eC) 
plotResiduals(res, form = growthforms$Dung) 
plotResiduals(res, form = growthforms$resid_fire)


### library(betareg)
###  betagrass_betareg <- betareg(forb ~ eC  + Dung + resid_fire, data = growthforms) 
###  summary(betagrass_betareg) 
###  # conclusion: my pseudo-R2 is identical to the betareg pdeudo R2 given above


betadwarfshrub <- glmmTMB(dwarfshrub ~ Dung + resid_fire + eC , data = growthforms, # pH excluded due to collinearity with soil eC
                     family = beta_family(link = "logit"),
                     ziformula = ~.
                     ) 
summary(betadwarfshrub) #no fire: AIC 0.3, no dung: AIC -17.4, both= -13.5
cor(growthforms$dwarfshrub, predict(betadwarfshrub, type="response"))^2
round(p.adjust(summary(betadwarfshrub)$coefficients$cond[,4], method = "fdr"),4)
round(p.adjust(summary(betadwarfshrub)$coefficients$zi[,4],   method = "fdr"),4)
res <- simulateResiduals(betadwarfshrub)
plot(res) 
plotResiduals(res, form = growthforms$eC) 
plotResiduals(res, form = growthforms$Dung) 
plotResiduals(res, form = growthforms$resid_fire)


betaforb <- glmmTMB(forb ~ Dung + resid_fire + eC , data = growthforms, # pH excluded due to collinearity with soil eC
                          family = beta_family(link = "logit")
                    ) 
summary(betaforb) #no fire: AIC -832.8, no dung: AIC -829.6, both: -838.2
pseudoR2(betaforb, link ="logit")
round(p.adjust(summary(betaforb)$coefficients$cond[,4], method = "fdr"),4)
res <- simulateResiduals(betaforb)
plot(res) 
plotResiduals(res, form = growthforms$eC) 
plotResiduals(res, form = growthforms$Dung) 
plotResiduals(res, form = growthforms$resid_fire)


betawoodyforb <- glmmTMB(forb_woody ~ Dung + resid_fire + eC , data = growthforms, # pH excluded due to collinearity with soil eC
                    family = beta_family(link = "logit")
                    ) 
summary(betawoodyforb) #no fire: AIC -612.0, no dung: AIC -658.0, both : -657.7
pseudoR2(betawoodyforb, link ="logit")
round(p.adjust(summary(betawoodyforb)$coefficients$cond[,4], method = "fdr"),4)
res <- simulateResiduals(betawoodyforb)
plot(res) 
plotResiduals(res, form = growthforms$eC) 
plotResiduals(res, form = growthforms$Dung) 
plotResiduals(res, form = growthforms$resid_fire)



### model the leaf trait CWMs ----------------------
modSLA <- glmmTMB(SLA ~  Dung + resid_fire + eC , data = CWMleaf, # pH excluded due to collinearity with soil eC
                  family = gaussian(link="identity")
                  )
summary(modSLA) #no fire: AIC 830.0, no dung: AIC 825.7 both : 826.5
pseudoR2(modSLA, link="identity")
#round(p.adjust(summary(modSLA)$coefficients$cond[,4], method = "fdr", n = 4),4)
res <- simulateResiduals(modSLA)
plot(res)
plotResiduals(res, form = CWMleaf$eC)
plotResiduals(res, form = CWMleaf$Dung)
plotResiduals(res, form = CWMleaf$resid_fire)



modLDMC <- glmmTMB(LDMC ~ Dung + resid_fire + eC + I(eC^2) , data = CWMleaf, # pH excluded due to collinearity with soil eC
                   family = gaussian(link="identity")
                   )
summary(modLDMC) #no fire: AIC 2106.5, no dung: AIC 2071.0, both : 2078.8
pseudoR2(modLDMC, link="identity")
#round(p.adjust(summary(modLDMC)$coefficients$cond[,4], method = "fdr", n = 4),4)
res <- simulateResiduals(modLDMC)
plot(res)
plotResiduals(res, form = CWMleaf$eC) 
plotResiduals(res, form = CWMleaf$Dung) 
plotResiduals(res, form = CWMleaf$resid_fire)


modheight <- glmmTMB(CWMheight ~ Dung + resid_fire + eC  , data = CWMheight, # pH excluded due to collinearity with soil eC
                     family = Gamma(link="log")
)
summary(modheight) #no fire: AIC 2106.5, no dung: AIC 2071.0, both : 2078.8
pseudoR2(modheight, link = "log")
#round(p.adjust(summary(modheight)$coefficients$cond[,4], method = "fdr", n = 4),4)
res <- simulateResiduals(modheight)
plot(res) 
plotResiduals(res, form = CWMheight$eC) 
plotResiduals(res, form = CWMheight$Dung) 
plotResiduals(res, form = CWMheight$resid_fire)




### make predictions with the ggeffects package -----------------------------

# label scaled predicotr on original scale
xlabel_fire_at <- c( ( -2 - mean(cons$resid_fire)) / sd(cons$resid_fire), ( 0 - mean(cons$resid_fire)) / sd(cons$resid_fire), ( 2 - mean(cons$resid_fire)) / sd(cons$resid_fire) ) 
xlabel_fire <- c(-2, 0, 2)

xlabel_eC_at <- c( (log(20) - mean(cons$eC)) / sd(cons$eC), (log(200) - mean(cons$eC)) / sd(cons$eC), (log(2000) - mean(cons$eC)) / sd(cons$eC) )
xlabel_eC <- c(20, 200, 2000)

xlabel_dung_at <- c( (log(0+1) - mean(cons$Dung)) / sd(cons$Dung), (log(5+1) - mean(cons$Dung)) / sd(cons$Dung), (log(25+1) - mean(cons$Dung)) / sd(cons$Dung) ) 
xlabel_dung <- c(0, 5, 25)



ylim_dwarf <- range(growthforms$dwarfshrub)
ybreak_dwarf <- c(0, 0.2, 0.4)
ylab_dwarf   <- c(0,20,40)

ylim_grass <- range(growthforms$grass)
ybreak_grass <- c(0, 0.3, 0.6)
ylab_grass   <- c(0,30,60)

ylim_forb <- range(growthforms$forb)
ybreak_forb <- c(0, 0.1, 0.2)
ylab_forb   <- c(0,10,20)

ylim_woodyforb <- range(growthforms$forb_woody)
ybreak_woodyforb <- c(0, 0.2, 0.4)
ylab_woodyforb   <- c(0,20,40)

ylim_SLA <- range(CWMleaf$SLA)
ybreak_SLA <- c(5, 10, 15)

ylim_LDMC <- range(CWMleaf$LDMC)
ybreak_LDMC <- c(300, 400)

ylim_height <- range(CWMheight$CWMheight)
ybreak_height <- c(10, 20, 30)



xlab_p_fire <- max(growthforms$resid_fire)*0.85
xlab_p_dung <- max(growthforms$Dung)*0.85


# set general plotting theme
my_theme <-   theme_bw() + theme(panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 plot.title = element_text(size=8),
                                 axis.title.y=element_text(colour="black", size=8) ,
                                 axis.title.x=element_text(colour="black", size=8) ,
                                 axis.ticks.x=element_line(colour="black", size=0.5) ,
                                 axis.text=element_text(size=8, colour="black"),
                                 panel.border=element_rect(colour="black", fill=NA) ) 


p_grass_fire <- plot(ggpredict(betagrass, "resid_fire [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=resid_fire, y=grass), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_grass)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_grass, breaks = ybreak_grass, labels = ylab_grass) +
  labs(title = NULL, x = NULL, y = NULL) + 
  my_theme
  
p_grass_eC       <- plot(ggpredict(betagrass, "eC [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=eC, y=grass), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_grass)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_grass, breaks = ybreak_grass, labels = ylab_grass) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_grass_dung       <- plot(ggpredict(betagrass, "Dung [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=Dung, y=grass), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_grass)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_grass, breaks = ybreak_grass, labels = ylab_grass) +
  labs(title = NULL,  x = NULL, y = "Grass cover (%)") +
  my_theme

p_dwarf_fire     <- plot(ggpredict(betadwarfshrub, "resid_fire [all]", type = "fe"))       +  
  geom_point(data = growthforms, aes(x=resid_fire, y=dwarfshrub), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_dwarf)*0.975, label= "p=0.02", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_dwarf, breaks = ybreak_dwarf, labels = ylab_dwarf) +
  labs(title = NULL,  x = NULL, y = NULL) +
  my_theme

p_dwarf_eC       <- plot(ggpredict(betadwarfshrub, "eC [all]", type = "fe"))            +  
  geom_point(data = growthforms, aes(x=eC, y=dwarfshrub), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_dwarf)*0.975, label= "p=0.01", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_dwarf, breaks = ybreak_dwarf, labels = ylab_dwarf) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_dwarf_dung       <- plot(ggpredict(betadwarfshrub, "Dung [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=Dung, y=dwarfshrub), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_dwarf)*0.975, label= "p=0.04", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_dwarf, breaks = ybreak_dwarf, labels = ylab_dwarf) +
  labs(title = NULL,  x = NULL, y = "Dwarfshrub cover (%)") +
  my_theme

p_forb_fire      <- plot(ggpredict(betaforb, "resid_fire [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=resid_fire, y=forb), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_forb)*0.975, label= "p=0.004", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_forb, breaks = ybreak_forb, labels = ylab_forb) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_forb_eC        <- plot(ggpredict(betaforb, "eC [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=eC, y=forb), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_forb)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_forb, breaks = ybreak_forb, labels = ylab_forb) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_forb_dung       <- plot(ggpredict(betaforb, "Dung [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=Dung, y=forb), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_forb)*0.975, label= "p=0.04", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_forb, breaks = ybreak_forb, labels = ylab_forb) +
  labs(title = NULL, x = NULL, y = "Forb cover (%)") +
  my_theme

p_woodyforb_fire <- plot(ggpredict(betawoodyforb, "resid_fire [all]", type = "fe"))        +  
  geom_point(data = growthforms, aes(x=resid_fire, y=forb_woody), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_woodyforb)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_woodyforb, breaks = ybreak_woodyforb, labels = ylab_woodyforb) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_woodyforb_eC   <- plot(ggpredict(betawoodyforb, "eC [all]", type = "fe"))       +  
  geom_point(data = growthforms, aes(x=eC, y=forb_woody), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_woodyforb)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_woodyforb, breaks = ybreak_woodyforb, labels = ylab_woodyforb) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_woodyforb_dung       <- plot(ggpredict(betawoodyforb, "Dung [all]", type = "fe"))      +  
  geom_point(data = growthforms, aes(x=Dung, y=forb_woody), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_woodyforb)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_woodyforb, breaks = ybreak_woodyforb, labels = ylab_woodyforb) +
  labs(title = NULL, x = NULL, y = "Woody forb cover (%)") +
  my_theme

p_SLA_eC <- plot(ggpredict(modSLA, "eC [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=eC, y=SLA), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_SLA)*0.975, label= "p=0.009", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_SLA, breaks = ybreak_SLA) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_SLA_fire <- plot(ggpredict(modSLA, "resid_fire [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=resid_fire, y=SLA), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_SLA)*0.975, label= "p=0.02", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_SLA, breaks = ybreak_SLA) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_SLA_dung <- plot(ggpredict(modSLA, "Dung [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=Dung, y=SLA), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_SLA)*0.975, label= "p=0.01", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_SLA, breaks = ybreak_SLA) +
  labs(title = NULL, x = NULL, y = "SLA (mm²/mg)") +
  my_theme

p_LDMC_eC <- plot(ggpredict(modLDMC, "eC [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=eC, y=LDMC), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_LDMC)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_LDMC, breaks = ybreak_LDMC) +
  labs(title = NULL, x = "Soil eC (µS)", y = NULL) +
  my_theme

p_LDMC_fire <- plot(ggpredict(modLDMC, "resid_fire [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=resid_fire, y=LDMC), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_LDMC)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_LDMC, breaks = ybreak_LDMC) +
  labs(title = NULL, x = "Residual of fire frequency", y = NULL) +
  my_theme

p_LDMC_dung <- plot(ggpredict(modLDMC, "Dung [all]", type = "fe"))          + 
  geom_point(data = CWMleaf, aes(x=Dung, y=LDMC), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_LDMC)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_LDMC, breaks = ybreak_LDMC) +
  labs(title = NULL, x = "No. dung piles (per 200m²)", y = "LDMC (mg/g)") +
  my_theme

p_height_eC <- plot(ggpredict(modheight, "eC [all]", type = "fe"))          + 
  geom_point(data = CWMheight, aes(x=eC, y=CWMheight), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=max(xlabel_eC_at), y=max(ylim_height)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_eC_at, labels = xlabel_eC) + 
  scale_y_continuous(limits = ylim_height, breaks = ybreak_height) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_height_fire <- plot(ggpredict(modheight, "resid_fire [all]", type = "fe"))          + 
  geom_point(data = CWMheight, aes(x=resid_fire, y=CWMheight), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_fire, y=max(ylim_height)*0.975, label= "p=0.68", size = 2.75) +
  scale_x_continuous(breaks = xlabel_fire_at, labels = xlabel_fire) + 
  scale_y_continuous(limits = ylim_height, breaks = ybreak_height) +
  labs(title = NULL, x = NULL, y = NULL) +
  my_theme

p_height_dung <- plot(ggpredict(modheight, "Dung [all]", type = "fe"))          + 
  geom_point(data = CWMheight, aes(x=Dung, y=CWMheight), alpha = 0.3, shape = 16, size = 1) +
  annotate("text", x=xlab_p_dung, y=max(ylim_height)*0.975, label= "p<0.001", size = 2.75) +
  scale_x_continuous(breaks = xlabel_dung_at, labels = xlabel_dung) + 
  scale_y_continuous(limits = ylim_height, breaks = ybreak_height) +
  labs(title = NULL, x = NULL, y = "Plant height (cm)") +
  my_theme






  
ggarrange(p_grass_dung    , 
          p_grass_fire    , 
          p_grass_eC      , 
          p_forb_dung     , 
          p_forb_fire     , 
          p_forb_eC       , 
          p_woodyforb_dung,
          p_woodyforb_fire,
          p_woodyforb_eC  ,
          p_dwarf_dung    , 
          p_dwarf_fire    , 
          p_dwarf_eC      , 
          p_height_dung   , 
          p_height_fire   , 
          p_height_eC     , 
          p_SLA_dung      , 
          p_SLA_fire      , 
          p_SLA_eC        , 
          p_LDMC_dung     , 
          p_LDMC_fire     , 
          p_LDMC_eC       , 
        #  labels = c("(a)", "", "",  
        #             "(b)", "", "",  
        #             "(c)", "", "",  
        #             "(d)", "", "",
        #             "(e)", "", "",
        #             "(f)", "", "",
        #             "(g)", "", "" 
        #  ),
        #  font.label = list(size = 12),
          hjust = -2,
          ncol = 3, nrow = 7
          ) 

ggsave("03-traits_bivariate.pdf", plot=last_plot(), device="pdf",
       width=18, height=24, units="cm", dpi=1200)

ggsave("03-traits_bivariate.png", plot=last_plot(), device="png",
       width=18, height=24, units="cm", dpi=1200)


