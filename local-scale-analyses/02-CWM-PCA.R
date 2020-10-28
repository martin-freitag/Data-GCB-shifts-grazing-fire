


# Retrieving Raw Data and Cleaning ----------------------------------------

envST <- read.table("02_envST.csv", sep=",", header=T, row.names=1)
speST <- read.table("01_speST.csv", sep=",", header=T, row.names=1)
traST <- read.table("03_traST_cat.csv", sep=",", header=T, row.names=1, na.strings= c("", "NA"))

# add the residuals of the firefrequency ~ dungpiles regression to the data!
envST$resid_fire <- resid_fire


library(vegan)
colnames(speST) <- make.cepnames(colnames(speST)) # abbreviate species names
rownames(traST) <- make.cepnames(rownames(traST))

library(ade4)

##prepare data frame containing constraining variables
cons               <- envST[, -c(1 : ncol(envST))]
cons$FireFreq      <- log(envST$fire_freq_check + 1)
cons$resid_fire    <- envST$resid_fire
cons$TimeSinceFire <- log(envST$timesincefire_check + 1)
cons$Dung          <- log(envST$dung_livestock + 1)  # think about logging Dung
cons$eC            <- log( (envST$LF_530 + envST$LF_3060)/2) #mean eC. think about log-transforming
cons$pH            <- envST$pH_530
cons$Sand          <- envST$Sand_530

##check, if trait information of at least 80% cover is available. SLA is available for the same species as LDMC
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


# for now, use the full set with NAs
traits <- traST
traits <- traits[row.names(traits)%in%names(species) , ]
# 92 species with full set, i.e. with SLA and LDMC


### calculate Community weighted means CWMs -----------------

library(FD)
library(ggplot2)
library(ggrepel)
library(ggvegan) #for fortify.cca function to extract predictor scores 


traits <- traits[c("LDMC", "SLA", "growthform", "life_history")] #excluded clonal strategy

# species data must be matrix for functcomp function

CWM <- functcomp(traits, as.matrix(log(species+1)), CWM.type="all")

# we do not use global height values for species, but have measured height ourselves
dimcwm <- read.table("03_PlantDimST.csv", sep=",", head=T)
summary(dimcwm)
dimcwm$meanH <- apply(dimcwm[,c("HeightIndividual1", "HeightIndividual2", "HeightIndividual3", "HeightIndividual4")], 
                      MARGIN = 1, function(x) mean(x, na.rm = TRUE))
head(dimcwm) #correct

dimcwm$log_cover <- log(dimcwm$cover + 1)

dimcwm$covxmH <- dimcwm$log_cover * dimcwm$meanH #weight height by cover (log-transformed)

height <- with(dimcwm, by(covxmH, plot, sum) / by(log_cover, plot, sum))
#procedure was tested by hand in Excel, correct
height <- data.frame(height = as.vector(height))
rownames(height) <- rownames(envST) #unnecessary step, but I'll keep this the way it is...
height <- height[ row.names(height)%in%row.names(CWM) , ] #subset


# and add height to the CWM data frame
CWM$height <- height
CWM <- apply(CWM, MARGIN = 2, scale )





### make a simple PCA with overlay of environmental variables -------------

CWMpca <- rda(CWM)

# importance (eigenvalue) first and second axis
summary(CWMpca)

# overlay of environmental variables
pca_env <- envfit(CWMpca ~  FireFreq + TimeSinceFire + Dung + eC + pH + Sand, data = cons, perm = 4999)
pca_env

CWM_scores  <- fortify(CWMpca, axes=1:2, display="sp")
site_scores <- fortify(CWMpca, axes=1:2, display="wa")
env_scores <- fortify(pca_env)
names(env_scores)[2] <- "Score"

# axis scores are exactly the opposite compared to later RDA - multiply b< -1 here to change direction
CWM_scores[c("PC1")]  <- -1 * CWM_scores[c("PC1")]
site_scores[c("PC1")] <- -1 * site_scores[c("PC1")]
env_scores[c("PC1")]  <- -1 * env_scores[c("PC1")]


pd <- position_dodge(0.1) #non-overlapping points
L_env <- 1 # scale arrow length
L_CWM <- 0.4

ggplot(site_scores, aes(x = PC1, y = PC2, colour = Score)) +
 # coord_fixed(ratio = 0.75) +
  geom_point(position=pd,size=1.25) +
  geom_segment(data=CWM_scores, mapping=aes(x=0, y=0, xend=PC1*L_CWM, yend=PC2*L_CWM), 
               arrow=arrow(angle = 15, length = unit(0.1, "inches"), type = "closed"), size = 0.75) +
  geom_text_repel(data=CWM_scores, aes(label = Label, x = PC1*L_CWM , y = PC2*L_CWM) ) +
  geom_segment(data=env_scores, mapping=aes(x=0, y=0, xend=PC1*L_env, yend=PC2*L_env), 
               arrow=arrow(angle = 15, length = unit(0.1, "inches"), type = "closed"), size = 0.75) +
  geom_text_repel(data=env_scores, aes(label = Label, x = PC1*L_env  , y = PC2*L_env )) +
  scale_colour_manual(values=c("gray50", "black", "darkgreen") ) + 
  scale_x_continuous(limits = c(-1.25, 1.25), breaks = c(-1, 0, 1), labels = c(-1, 0, 1)) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(), 
                     legend.position="none",
                     axis.title.y=element_text(colour="black", size=12) ,
                     axis.title.x=element_text(colour="black", size=12) ,
                     axis.ticks.x=element_line(colour="black", size=0.5) ,
                     axis.text=element_text(size=11, colour="black"),
                     panel.border=element_rect(colour="black", fill=NA) )

#ggsave("02-traits_PCA.pdf", plot=last_plot(), device="pdf",
#       width=18, height=14, units="cm", dpi=1200)





### calculate a RDA with environmental variables as predictors ------------
### this is only to check how much variation in the CWM is explained by these three predictors

CWMrda <- rda(CWM[,c(1:8)] ~ FireFreq + Dung + eC, data = cons)
CWMrda


### calculate a DCA with environmental variables as overlay ------------



DCA <- decorana(log(species + 1), iweigh=1)
DCA_env <- envfit(DCA ~  FireFreq + TimeSinceFire + Dung + eC + pH + Sand, data = cons, perm = 4999)
DCA_env

species_scores  <- fortify(DCA, axes=1:2, display="sp")
site_scores     <- fortify(DCA, axes=1:2, display="sites")
env_scores      <- fortify(DCA_env)
names(env_scores)[2] <- "Score"

# label only the 30 species with highest weight, those values (from weightes averages) stored in DCA$adotj (thanks to Rike!)
#freq <- colSums(species != 0)
species_scores <- species_scores[order(DCA$adotj, decreasing = TRUE)[1:30] , ]




### # same with NMDS, see if there is a difference
### NMDS <- metaMDS(log(species + 1), distance = "bray", trymax = 200, previous.best = NMDS)
### NMDS_env <- envfit(NMDS ~  FireFreq + TimeSinceFire + Dung + eC + pH + Sand, data = cons, perm = 4999)
### 
### species_scores  <- data.frame(scores(NMDS, display="sp"))
### species_scores$Score <- rep("species", nrow(species_scores))
### site_scores     <- data.frame(scores(NMDS, display="sites"))
### site_scores$Score <- rep("site", nrow(site_scores))
### env_scores      <- fortify(NMDS_env)
### names(env_scores)[2] <- "Score"




# add space in abbreviated species names
species_scores$Label <- paste(substr(species_scores$Label, 1, 4),
                                  substr(species_scores$Label, 5, 8),
                                  sep=" ")


pd <- position_dodge(0.1) #non-overlapping points
L <- 2 # scale arrow length

ggplot(site_scores, aes(x = DCA1, y = DCA2, colour = Score)) +
  coord_fixed(ratio = 0.75) +
  geom_point(position=pd,size=0.75) +
  geom_point(data=species_scores, position=pd,size=1, shape = 4) +
  geom_text_repel(data=species_scores, aes(label = Label, x = DCA1 , y = DCA2), size = 2.5) +
  geom_segment(data=env_scores, mapping=aes(x=0, y=0, xend=DCA1*L, yend=DCA2*L), 
               arrow=arrow(angle = 15, length = unit(0.1, "inches"), type = "closed")) +
  geom_text_repel(data=env_scores, aes(label = Label, x = DCA1*L  , y = DCA2*L )) +
  scale_colour_manual(values=c("gray50", "gray25", "darkgreen") ) + 
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(), 
                     legend.position="none",
                     axis.title.y=element_text(colour="black", size=10) ,
                     axis.title.x=element_text(colour="black", size=10) ,
                     axis.ticks.x=element_line(colour="black", size=0.5) ,
                     axis.text=element_text(size=10, colour="black"),
                     panel.border=element_rect(colour="black", fill=NA) )

#ggsave("02-species_DCA.pdf", plot=last_plot(), device="pdf",
#       width=18, height=14, units="cm", dpi=1200)








### # Fourth-Corner Statistics ------------------------------------------------
### 
### #Combine Fire Frequency and Time-Since_Fire to one Fire variable using PCA
### 
### env4th <- cbind(cons[c("Dung", "eC")]) #dataframe to be used in the fourthcorner statistic
### 
### # if needed: remove species that do not have a full trait set
### tra4th <- traits[complete.cases(traits) , ]
### spe4th <- species[,names(species) %in% row.names(tra4th)] #keep species which are kept in traST
### tra4th <- tra4th[row.names(tra4th)%in%names(spe4th) , ]
### 
### 
### library(ade4)
### 
### four.comb <- fourthcorner(env4th, spe4th,
###                           tra4th, modeltype = 6, #modeltype 6: permuting samples and species to get a correct alpha-error
###                           p.adjust.method.G = "none",  #Holm's correction accounts for linearly dependent factor levels
###                           p.adjust.method.D = "none", nrepet = 4999, 
###                           p.adjust.D = "levels") #p-values adjusted for number of levels, not globally for all tests
### 
### 
### plot(four.comb, alpha = 0.05, stat = "D2")
### print(four.comb, stat = "D2") #significance codes







