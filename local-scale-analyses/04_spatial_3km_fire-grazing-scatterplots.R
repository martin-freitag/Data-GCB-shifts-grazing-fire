

spatial_3km <- read.table("04_Grazing_BA_3km_means.csv", sep=",", dec=".", header=T)
str(spatial_3km)

# bring to % scale
spatial_3km$GRAZ1990_mean <- spatial_3km$GRAZ1990_mean/100
spatial_3km$GRAZ2000_mean <- spatial_3km$GRAZ2000_mean/100
spatial_3km$GRAZ2015_mean <- spatial_3km$GRAZ2015_mean/100
spatial_3km$GRAZ_CHANGE   <- spatial_3km$GRAZ_CHANGE  /100

spatial_3km$BA1990_mean <- spatial_3km$BA1990_mean * 100
spatial_3km$BA2015_mean <- spatial_3km$BA2015_mean * 100


# remove 3km pixels that have less than 50% grassland or negative grazing probability (at the edges of grazing maps)
spatial_3km <- subset(spatial_3km, GRASSLAND_mean >= 0.5 & GRAZ1990_mean >= 0 & GRAZ2015_mean >= 0)



library(ggplot2)
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






# grazing probability and burned area ca. 1990


p1 <- ggplot(data=spatial_3km, aes(x=GRAZ1990_mean, BA1990_mean)) +
  geom_hex( bins = 100) +   #aes(fill = log(..count..)) 
  scale_fill_viridis_c( ) +
  labs(title = NULL, x = "Grazing probability 1990 (%)", y = "Burned area 1990 (%)") + 
  my_theme
p1


# grazing probability and burned area ca. 2015


p2 <- ggplot(data=spatial_3km, aes(x=GRAZ2015_mean, BA2015_mean)) +
  geom_hex(bins = 100 ) +
  scale_fill_viridis_c( ) +
  labs(title = NULL, x = "Grazing probability 2015 (%)", y = "Burned area 2015 (%)") + 
  my_theme
p2


# grazing probability ca. 2015 fire frequency 2000-2015


p3 <- ggplot(data = spatial_3km, aes(x=GRAZ2015_mean, FireFreq_mean)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c( ) +
  labs(title = NULL, x = "Grazing probability 2015 (%)", y = "Fire frequency 2000-2015 (%)") + 
  my_theme
p3


library(ggpubr)

ggarrange(p1, p2, p3, 
          nrow=2, ncol=2, 
          labels = c("a)", "b)", "c)" ),
          hjust = -0.3
          )

ggsave("Figures/04-Bivariate-scatterplots.pdf", plot=last_plot(), device="pdf",
       width=18, height=18, units="cm", dpi=1200)

ggsave("Figures/04-Bivariate-scatterplots.png", plot=last_plot(), device="png",
       width=18, height=18, units="cm", dpi=1200)


