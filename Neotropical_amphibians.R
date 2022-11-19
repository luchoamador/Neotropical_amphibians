##############################################################
##############################################################
#Determinant of nucleotide diversity in Neotropical Amphibans#
#setwd("C:/Users/Luis Amador/Dropbox/Amphibia_Neotropical")

matrix <- read.csv(file = "amphibia_Neotropical.csv", sep = ";")
colnames(matrix)
matrix <- matrix
matrix[is.na(matrix) | matrix == "Inf"] <- NA # Replace NaN & Inf with NA

library(ggplot2)
library(ggpubr)

library(visdat)
library(tidyverse)
library(lattice)
library(DHARMa)
library(performance)
library(MuMIn)
library(piecewiseSEM)
library(MASS)
library(ggExtra)
library(Rmisc)
library(emmeans)
library(sjPlot)
library(bbmle)
library(glmmTMB)
library(ordinal)
library(car)
library(ecolottery)
library(naniar)
library(vcd)
library(gvlma)
library(generalhoslem)
library(randomForest)

#Density of Nucleotide Diversity
d <- density(matrix$pi)
plot(d, main="Kernel Density of Nucleotide Diversity")
polygon(d, col = "red", border = "blue")
rug(matrix$pi, col = "brown")
dl <- density(log(matrix$pi))
plot(dl, main="Kernel Density of Nucleotide Diversity")
polygon(dl, col = "red", border = "blue")
rug(matrix$pi, col = "brown")

############Linear Models####################################
###Amphibian body size
# Fit regression line
size <- lm(pi~Body_size_mm, data = matrix)
summary(size)
anova(size)

pi_size <- ggplot(matrix, aes(x=(Body_size_mm), y=pi,  shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Body size (mm)", 
                                                             y="Genetic diversity (pi)")+
  theme_minimal()

coeff=coefficients(size)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
eq
# Plot
pi_size + geom_abline(intercept = 2.793722e-02, slope = -1.953164e-05, linetype="dashed", size=0.5)+
  ggtitle(eq)

###Anura
anura <- read.csv(file = "Anura.csv", sep=";")
an1 <- lm(pi~Body_size_mm, data = anura)
anova(an1)
summary(an1)

anu_pi_size <- ggplot(anura, aes(x=(Body_size_mm), y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Body size (mm)", 
  y="Genetic diversity (pi)")+
  theme_minimal()

coeff=coefficients(an1)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu_pi_size + geom_abline(intercept = 2.922940e-02, slope = -3.657869e-05, 
                      linetype="dashed", size=0.5)+
  ggtitle(eq)

###############################
##Caudata
caudata <- read.csv(file = "Caudata.csv", sep = ";")
ca1 <- lm(pi~Body_size_mm, data = caudata)
anova(ca1)
summary(ca1)

cau_pi_size <- ggplot(caudata, aes(x=(Body_size_mm), y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Body size (mm)", 
  y="Genetic diversity (pi)")+
  theme_minimal()

coeff=coefficients(ca1)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau_pi_size + geom_abline(intercept = 0.0078384058, slope = 0.0001417394, 
                          linetype="dashed", size=0.5)+
  ggtitle(eq)


all1 <- pi_size + geom_abline(intercept = 2.793722e-02, slope = -1.953164e-05, 
                          linetype="dashed", size=0.5)+ggtitle(eq)
anu1 <- anu_pi_size + geom_abline(intercept = 2.922940e-02, slope = -3.657869e-05, 
                          linetype="dashed", size=0.5)+ggtitle(eq)
cau1 <- cau_pi_size + geom_abline(intercept = 0.0078384058, slope = 0.0001417394, 
                          linetype="dashed", size=0.5)+ggtitle(eq)

ggarrange(all1, anu1, cau1, labels = c("A", "B", "C"), ncol = 1, nrow = 3)
####################################################################
####ANURA Development mode
par(mfrow=c(1,1))

ADm <- lm(pi~Development_Mode, data = matrix)
anova(ADm)
summary(ADm)

all2 <- ggplot(matrix, aes(x=(Development_Mode), y=pi,shape=Order, fill=Order))+
  geom_boxplot(notch = TRUE)+labs(x="Development Mode", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#anura
aDm <- lm(pi~Development_Mode, data = anura)
anova(aDm)
summary(aDm)

anu2 <- ggplot(anura, aes(x=(Development_Mode), y=pi,shape=Order, fill=Order))+
  geom_boxplot(notch = TRUE)+labs(x="Development Mode", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#caudata

cau2 <- ggplot(caudata, aes(x=(Development_Mode), y=pi,shape=Order, fill=Order))+
  geom_boxplot(notch = TRUE)+labs(x="Development Mode", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

ggarrange(all2, anu2, cau2, labels = c("A", "B", "C"), ncol = 3, nrow = 1)

##############################################################
#######Amphibia habitat#######################################
am_h <- lm(pi ~ Microhabitat, data=matrix)
anova(am_h)
summary(am_h)

all3 <- ggplot(matrix, aes(x=Microhabitat, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Habitat", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#anura
an_h <- lm(pi ~ Microhabitat, data=anura)
anova(an_h)
summary(an_h)

anu3 <- ggplot(anura, aes(x=Microhabitat, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Habitat", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#Caudata
ca_h <- lm(pi ~ Microhabitat, data=caudata)
anova(ca_h)
summary(ca_h)

cau3 <- ggplot(caudata, aes(x=Microhabitat, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Habitat", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

ggarrange(all3, anu3, cau3, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

##############################################################
#######Amphibia activity#######################################
am_a <- lm(pi ~ Activity, data=matrix)
anova(am_a)
summary(am_a)

all4 <- ggplot(matrix, aes(x=Activity, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Activity", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#anura
an_a <- lm(pi ~ Activity, data=anura)
anova(an_a)
summary(an_a)

anu4 <- ggplot(anura, aes(x=Activity, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Activity", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

cau4 <- ggplot(caudata, aes(x=Activity, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Activity", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

#Caudata all species are nocturnal

ggarrange(all4, anu4, cau4, labels = c("A", "B", "C"), ncol = 3, nrow = 1)

##############################################################
#######Amphibia number of sequences#######################################

an_seqs <- lm(pi ~ Seqs, data=matrix)
anova(an_seqs)
summary(an_seqs)

all5 <- ggplot(matrix, aes(x=Seqs, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Number of sequences", 
             y="Genetic diversity (pi)")+ theme_minimal()

coeff=coefficients(an_seqs)
coeff
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all5 <- all5 + geom_abline(intercept = 0.0238057780, slope = 0.0001605617, 
    linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura
an_seqs <- lm(pi ~ Seqs, data=anura)
anova(an_seqs)
summary(an_seqs)

anu5 <- ggplot(anura, aes(x=Seqs, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Number of sequences", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_seqs)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu5 <- anu5 + geom_abline(intercept = 0.0242478420, slope = 0.0001623265, 
                    linetype="dashed", size=0.5)+ ggtitle(eq)

#Caudata
c_seqs <- lm(pi ~ Seqs, data=caudata)
anova(c_seqs)
summary(c_seqs)

c_sp <- ggplot(caudata, aes(x=Seqs, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Number of sequences", 
       y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(c_seqs)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau5 <- c_sp + geom_abline(intercept = 0.0233871917, slope = -0.0002070629, 
                   linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all5, anu5, cau5, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

##########################################################################
#######Amphibia Elevation mean#######################################
am_el <- lm(pi ~ Elevation_mean, data=matrix)
anova(am_el)
summary(am_el)

amele <- ggplot(matrix, aes(x=Elevation_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Elevation", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_el)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all6 <- amele + geom_abline(intercept = 3.090135e-02, slope = -4.449711e-06, 
                   linetype="dashed", size=0.5)+ ggtitle(eq)

#anura
an_el <- lm(pi ~ Elevation_mean, data=anura)
anova(an_el)
summary(an_el)

anele <- ggplot(anura, aes(x=Elevation_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Elevation", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_el)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu6 <- anele + geom_abline(intercept = 3.060783e-02, slope = -3.790475e-06, 
                            linetype="dashed", size=0.5)+ ggtitle(eq)

#caudata
ca_el <- lm(pi ~ Elevation_mean, data=caudata)
anova(ca_el)
summary(ca_el)

caele <- ggplot(caudata, aes(x=Elevation_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Elevation", 
             y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_el)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau6 <- caele + geom_abline(intercept = 3.370568e-02, slope = -7.218771e-06, 
                            linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all6, anu6, cau6, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

###Minimum_Elevation

am_minel <- lm(pi ~ Elevation_min, data=matrix)
anova(am_minel)
summary(am_minel)

amminele <- ggplot(matrix, aes(x=Elevation_min, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Elevation",                           y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_minel)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all7 <- amminele + geom_abline(intercept = 3.192748e-02, slope = -1.028587e-05, 
                    linetype="dashed", size=0.5)+ ggtitle(eq)

#anura
an_minel <- lm(pi ~ Elevation_min, data=anura)
anova(an_minel)
summary(an_minel)

anminele <- ggplot(anura, aes(x=Elevation_min, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Elevation",                           y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_minel)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu7 <- anminele + geom_abline(intercept = 3.184713e-02, slope = -1.052934e-05, 
                               linetype="dashed", size=0.5)+ ggtitle(eq)

#caudata
ca_minel <- lm(pi ~ Elevation_min, data=caudata)
anova(ca_minel)
summary(ca_minel)

caminele <- ggplot(caudata, aes(x=Elevation_min, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Elevation",                           y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_minel)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau7 <- caminele + geom_abline(intercept = 3.438494e-02, slope = -1.102808e-05, 
                              linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all7, anu7, cau7, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

####Maximum elevation

am_maxel <- lm(pi ~ Elevation_max, data=matrix)
anova(am_maxel)
summary(am_maxel)

ammaxele <- ggplot(matrix, aes(x=Elevation_max, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Elevation",                           y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_maxel)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all8 <- ammaxele + geom_abline(intercept = 2.804654e-02, slope = -2.671310e-08, 
                       linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura

an_mx_el <- lm(pi ~ Elevation_max, data=anura)
anova(an_mx_el)
summary(an_mx_el)

anmaxele <- ggplot(anura, aes(x=Elevation_max, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Elevation", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_mx_el)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu8 <- anmaxele + geom_abline(intercept = 2.778631e-02, slope = 8.242393e-07, 
                               linetype="dashed", size=0.5)+ ggtitle(eq)

#Caudata

ca_mx_el <- lm(pi ~ Elevation_max, data=caudata)
anova(ca_mx_el)
summary(ca_mx_el)
camaxele <- ggplot(caudata, aes(x=Elevation_max, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Elevation", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_mx_el)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau8 <- camaxele + geom_abline(intercept = 1.642412e-02, slope = 1.869884e-06, 
                               linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all8, anu8, cau8, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

##########################################################################
#######Amphibia LAtitude mean#######################################
am_lat <- lm(pi ~ Lat_mean, data=matrix)
anova(am_lat)
summary(am_lat)
amlat <- ggplot(matrix, aes(x=Lat_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Latitude", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_lat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all9 <- amlat + geom_abline(intercept = 0.0272700313, slope = 0.0001095448, 
                               linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura
an_lat <- lm(pi ~ Lat_mean, data=anura)
anova(an_lat)
summary(an_lat)
anlat <- ggplot(anura, aes(x=Lat_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Latitude", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_lat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu9 <- anlat + geom_abline(intercept = 0.0282673850, slope = 0.0001418137, 
                            linetype="dashed", size=0.5)+ ggtitle(eq)

#Caudata
ca_lat <- lm(pi ~ Lat_mean, data = caudata)
anova(ca_lat)
summary(ca_lat)

calat <- ggplot(caudata, aes(x=Lat_mean, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Mean Latitude", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_lat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau9 <- calat + geom_abline(intercept = 0.003758430, slope = 0.001263873, 
                            linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all9, anu9, cau9, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

#Minimum Latitude"

am_ilat <- lm(pi ~ Lat_realmin, data=matrix)
anova(am_ilat)
summary(am_ilat)

amminlat <- ggplot(matrix, aes(x=Lat_realmin, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Latitude", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_ilat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all10 <- amminlat + geom_abline(intercept = 2.721016e-02, slope = 7.198868e-05, 
                            linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura

an_ilat <- lm(pi ~ Lat_realmin, data=anura)
anova(an_ilat)
summary(an_ilat)

anminlat <- ggplot(anura, aes(x=Lat_realmin, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Latitude", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_ilat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu10 <- anminlat + geom_abline(intercept = 0.028173144, slope = 0.000105755, 
                                linetype="dashed", size=0.5)+ ggtitle(eq)
#Caudata

ca_ilat <- lm(pi ~ Lat_realmin, data=caudata)
anova(ca_ilat)
summary(ca_ilat)

caminlat <- ggplot(caudata, aes(x=Lat_realmin, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Min Latitude",          y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_ilat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau10 <- caminlat + geom_abline(intercept = 0.003217001, slope = 0.001268071, 
                                linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all10, anu10, cau10, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

###Maximum Latitude

am_xlat <- lm(pi ~ Lat_realmax, data=matrix)
anova(am_xlat)
summary(am_xlat)

amxlat <- ggplot(matrix, aes(x=Lat_realmax, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Latitude", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(am_xlat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all11 <- amxlat + geom_abline(intercept = 0.0272507497, slope = 0.0001281249, 
                                linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura

an_xlat <- lm(pi ~ Lat_realmax, data=anura)
anova(an_xlat)
summary(an_xlat)

anxlat <- ggplot(anura, aes(x=Lat_realmax, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Latitude", y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(an_xlat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu11 <- anxlat + geom_abline(intercept = 0.0282249798, slope = 0.0001541099, 
                              linetype="dashed", size=0.5)+ ggtitle(eq)

#Caudata

ca_xlat <- lm(pi ~ Lat_realmax, data=caudata)
anova(ca_xlat)
summary(ca_xlat)
caxlat <- ggplot(caudata, aes(x=Lat_realmax, y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Max Latitude",          y="Genetic diversity (pi)")+  theme_minimal()
coeff=coefficients(ca_xlat)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau11 <- caxlat + geom_abline(intercept = 0.005228133, slope = 0.001190412, 
                              linetype="dashed", size=0.5)+ ggtitle(eq)


ggarrange(all11, anu11, cau11, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

##########################################################################
#######Amphibia Species range#######################################
am_sr <- lm(pi ~ Area_km2, data=matrix)
anova(am_sr)
summary(am_sr)

sr1 <- ggplot(matrix, aes(x=log(Area_km2), y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Range size (log)", 
                     y="Genetic diversity (pi)")+  theme_minimal()

coeff=coefficients(am_sr)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
all12 <- sr1 + geom_abline(intercept = 2.510185e-02, slope = 3.408207e-09, 
                  linetype="dashed", size=0.5)+ ggtitle(eq)

#Anura
an_sr <- lm(pi ~ Area_km2, data=anura)
anova(an_sr)
summary(an_sr)
#Residual standard error: 0.03183 on 228 degrees of freedom
#Multiple R-squared:  0.1082,	Adjusted R-squared:  0.1043 
#F-statistic: 27.67 on 1 and 228 DF,  p-value: 3.315e-07 ***
sr2 <- ggplot(anura, aes(x=log(Area_km2), y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Range size (log)", 
                          y="Genetic diversity (pi)")+ theme_minimal()

coeff=coefficients(an_sr)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
anu12 <- sr2 + geom_abline(intercept = 2.561003e-02, slope = 3.305602e-09, 
                  linetype="dashed", size=0.5)+ ggtitle(eq)


#Caudata
ca_sr <- lm(pi ~ Area_km2, data=caudata)
anova(ca_sr)
summary(ca_sr)

sr3 <- ggplot(caudata, aes(x=log(Area_km2), y=pi, shape=Order, fill=Order))+
  geom_point(size=5, color="black",shape=21, alpha=0.7)+labs(x="Range size (log)", 
  y="Genetic diversity (pi)")+ theme_minimal()
coeff=coefficients(ca_sr)
coeff
# Equation of the line : 
eq = paste0("y = ", coeff[2],1, "*x + ", coeff[1],1)
# Plot
# Change line type, color and size
cau12 <- sr3 + geom_abline(intercept = 1.784711e-02, slope = 1.546709e-07, 
                  linetype="dashed", size=0.5)+ ggtitle(eq)

ggarrange(all12, anu12, cau12, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

#########################################################################
############ Orders##################################################
an_o <- lm(pi ~ Order, data=matrix)
anova(an_o)
summary(an_o)

ggplot(matrix, aes(x=Order, y=pi,shape=Order, fill=Order))+
  geom_boxplot(notch = FALSE)+labs(x="Order", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", alpha=0.5)+
  theme_minimal()

###Family

#############Amphibian families

fam <- lm(pi ~ Family, data=matrix)
anova(fam)
fams <- aov(matrix$pi ~ matrix$Family)
anova(fams)
summary(fams)
coefficients(fams)
TukeyHSD(fams, conf.level = 0.95)

matrix$Family <- factor(matrix$Family, levels=c("Ranidae", "Leptodactylidae","Bufonidae", "Odontophrynidae","Hylidae","Alsodidae","Hylodidae", "Batrachylidae","Cycloramphidae",
 "Telmatobiidae", "Rhinodermatidae","Strabomantidae", "Craugastoridae", "Eleutherodactylidae",  "Brachycephalidae",  "Aromobatidae",  "Dendrobatidae","Plethodontidae"))

ggplot(matrix, aes(x=Family, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Family", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", alpha=0.5)+
  scale_fill_manual(values = c("#F08080", "#20B2AA"), name="Order")+
  theme_minimal()+theme(axis.text.x = element_text(size = 6.5)) 

#Anuran families
an.fam <- lm(pi ~ Family, data=anura)
anova(an.fam)
summary(an.fam)

####Conservation status
uicn <- lm(pi ~ IUCN, data=matrix)
anova(uicn)
summary(uicn)

all13 <- ggplot(matrix, aes(x=IUCN, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="IUCN Conservation Status", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

uicn_a <- lm(pi ~ IUCN, data=anura)
anova(uicn_a)
summary(uicn_a)

anu13 <- ggplot(anura, aes(x=IUCN, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="IUCN Conservation Status", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

uicn_c <- lm(pi ~ IUCN, data=caudata)
anova(uicn_c)
summary(uicn_c)

cau13 <- ggplot(caudata, aes(x=IUCN, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="IUCN Conservation Status", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

ggarrange(all13, anu13, cau13, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

#####Biogeographic_Unit#####
bu1 <- lm(pi ~ Biogeographic_Unit, data = matrix)
anova(bu1)
summary(bu1)

all14 <- ggplot(matrix, aes(x=Biogeographic_Unit, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Biogeographic Unit", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

bu2 <- lm(pi ~ Biogeographic_Unit, data = anura)
anova(bu2)
summary(bu2)

anu14 <- ggplot(anura, aes(x=Biogeographic_Unit, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Biogeographic Unit", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

bu3 <- lm(pi ~ Biogeographic_Unit, data = caudata)
anova(bu3)
summary(bu3)

cau14 <- ggplot(caudata, aes(x=Biogeographic_Unit, y=pi,shape=Order, fill=Order))+
  geom_boxplot()+labs(x="Biogeographic Unit", y="Genetic diversity (pi)")+
  scale_shape_manual(values = c(21,21))+
  geom_point(position = "jitter", fill="white", size=2)+
  theme_minimal()

ggarrange(all14, anu14, cau14, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

##########################################
###combined plots#####
# escala de colores perceptualmente uniforme 
ggplot(matrix,aes(x=Elevation_min,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(matrix,aes(x=Lat_mean,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(matrix,aes(log(x=Area_km2),y=pi,fill=Lat_mean))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(matrix,aes(log(x=Area_km2),y=pi,fill=Elevation_mean))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(matrix,aes(log(x=Area_km2),y=pi,fill=Elevation_min))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()



#Anura
ggplot(anura,aes(x=Lat_mean,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+ 
  scale_fill_viridis_c()+theme_bw()

ggplot(anura,aes(log(x=Area_km2),y=pi,fill=Lat_mean))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(anura,aes(x=Elevation_min,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(anura,aes(log(x=Area_km2),y=pi,fill=Elevation_min))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()


#Caudata
ggplot(caudata,aes(x=Lat_mean,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(caudata,aes(log(x=Area_km2),y=pi,fill=Elevation_min))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(caudata,aes(x=Elevation_min,y=pi,fill=Area_km2))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

ggplot(caudata,aes(log(x=Area_km2),y=Lat_mean,fill=pi))+
  geom_point(color="black",size=5,shape=21,alpha=0.7)+
  scale_fill_viridis_c()+theme_bw()

##Multiple linear regression
mlr1 <- lm(pi ~ Lat_mean + Area_km2, data = matrix)
summary(mlr1)
#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 2.539e-02  2.126e-03  11.943   <2e-16 ***
#  Lat_mean    1.151e-04  1.054e-04   1.092   0.2760    
#Area_km2    3.433e-09  1.207e-09   2.845   0.0048 ** 
#Residual standard error: 0.03202 on 253 degrees of freedom
#Multiple R-squared:  0.03499,	Adjusted R-squared:  0.02736 
#F-statistic: 4.587 on 2 and 253 DF,  p-value: 0.01104

mlr2 <- lm(pi ~ Area_km2 + Lat_mean, data = anura)
summary(mlr2)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 2.624e-02  2.377e-03  11.040  < 2e-16 ***
#  Area_km2    3.285e-09  1.257e-09   2.614  0.00955 ** 
#  Lat_mean    1.378e-04  1.146e-04   1.202  0.23060    
#Residual standard error: 0.03317 on 227 degrees of freedom
#Multiple R-squared:  0.03553,	Adjusted R-squared:  0.02703 
#F-statistic: 4.181 on 2 and 227 DF,  p-value: 0.01648

mlr3 <- lm(pi ~ Area_km2 + Lat_mean, data = caudata)
summary(mlr3)
#Coefficients:
# Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 3.269e-03  8.889e-03   0.368   0.7164  
#Area_km2    1.395e-07  6.214e-08   2.244   0.0347 *
#  Lat_mean    1.078e-03  6.036e-04   1.786   0.0873 .
#Residual standard error: 0.01659 on 23 degrees of freedom
#Multiple R-squared:  0.2924,	Adjusted R-squared:  0.2309 
#F-statistic: 4.752 on 2 and 23 DF,  p-value: 0.01873

#Interaction
mlr4 <- lm(pi ~ Area_km2 + Lat_mean + Area_km2:Lat_mean, data = anura)
summary(mlr4)

mlr5 <- lm(pi ~ Area_km2 + Lat_mean + Area_km2:Lat_mean, data = caudata)
summary(mlr5)




##################Nucleotide diversity#################


###IMPORT SEQUENCE DATA
#setwd("~/Downloads/phylogatr-results/Amphibia/")
setwd("~/Documents/NuevoTropico//")
myseqs <- list.files(path = "Neotr_amphibians_CytBSeqs/", pattern = ".fasta", full.names = TRUE, recursive = TRUE)

#Calculate genetic diversity
for (f in myseqs) {
  
  #get species name and gene name
  sp <- basename(f)
  sp <- file_path_sans_ext(sp)
  
  #read in fasta file
  seq <- fasta2DNAbin(f, quiet = TRUE)
  
  #get number of sequence
  n <- nrow(seq)
  
  #calculate nucleotide diversity
  pi <- nuc.div(seq)
  
  #write data to file
  write.table(data.frame(sp, pi, n), file = "neotropical_amphibians_pi.csv", quote = FALSE, row.names = FALSE, 
              col.names = !file.exists("neotropical_amphibians_pi.csv"), append = TRUE, sep = ",")
}

warnings()
##Outlier analysis
#First import the output from previous step
data <- read.csv("neotropical_amphibians_pi.csv")

lower_bound <- quantile(data$pi, 0.025)
lower_bound

upper_boud <- quantile(data$pi, 0.975)
upper_boud

outliers <- which(data$pi < lower_bound | data$pi > upper_boud)
data[outliers,]

# nucleotide diversity Is defined as the average number of #nucleotide differences per site in pairwise comparisons among DNA sequences



##################Random Forests###########################
###########################################################
#All amphibians
names(matrix) <- c("Species_name","pi","Microhabitat", 
                   "Leaves","Flowers","Seeds","Fruits","Arthro","Vert", 
                   "Activity","Wet_warm","Wet_cold","Dry_warm",              
                   "Dry_cold","Body_mass_g","Age_at_maturity_min_y",
                   "Age_at_maturity_max_y","Body_size_mm",
                   "Size_at_maturity_min_mm","Size_at_maturity_max_mm",
                   "Longevity_max_y","Litter_size_min_n","Litter_size_max_n",
                   "Reproductive_output_y","Offspring_size_min_mm",          
                   "Offspring_size_max_mm", "Development_Mode","Seqs",
                   "Max_Species_SVL_mm","Max_MaleSVL_mm","Max_FemaleSVL_mm",
                   "Vertebrate_Recorded_In_Diet.","Invertebrate_Recorded_In_Diet.", 
                   "PlantMaterial_Recorded_In_Diet.","Elevation_mean",     
                   "Elevation_min","Elevation_max","Lat_mean","Lat_realmax",
                   "Lat_realmin","Area_km2","Order","Family","Genus","Species",
                   "Scientific_name","Biogeographic_Unit","OBS","IUCN")
colnames(matrix)
df <- matrix[c(-1,-4,-5,-6,-7,-8,-9,-11,-12,-13,-14,-15,-16,-17,-19,-20,-21,-22,
               -23,-24,-25,-26,-29,-30,-31,-32,-33,-34,-35,-37,-39,-40,
               -44,-45,-46,-48)]
matrix$pi <- as.factor(matrix$pi)

set.seed(3483)
fit_1 <- randomForest(pi~.,data = df, proximity=TRUE, importance=TRUE, 
                      ntree=1000, nPerm=10)
fit_1

varImpPlot(fit_1)
## Look at variable importance:
round(importance(fit_1), 2)


fit_1$importance
fit_1$importanceSD
fit_1$y


#####Anura random forests######
names(anura) <- c("Species_name","pi","Microhabitat", 
                  "Leaves","Flowers","Seeds","Fruits","Arthro","Vert", 
                  "Activity","Wet_warm","Wet_cold","Dry_warm",              
                  "Dry_cold","Body_mass_g","Age_at_maturity_min_y",
                  "Age_at_maturity_max_y","Body_size_mm",
                  "Size_at_maturity_min_mm","Size_at_maturity_max_mm",
                  "Longevity_max_y","Litter_size_min_n","Litter_size_max_n",
                  "Reproductive_output_y","Offspring_size_min_mm",          
                  "Offspring_size_max_mm", "Development_Mode","Seqs",
                  "Max_Species_SVL_mm","Max_MaleSVL_mm","Max_FemaleSVL_mm",
                  "Vertebrate_Recorded_In_Diet.","Invertebrate_Recorded_In_Diet.", 
                  "PlantMaterial_Recorded_In_Diet.","Elevation_mean",     
                  "Elevation_min","Elevation_max","Lat_mean","Lat_realmax",
                  "Lat_realmin","Area_km2", "Order","Family","Genus","Species",
                  "Scientific_name","Biogeographic_Unit","OBS","IUCN")
df_A <- anura[c(-1,-4,-5,-6,-7,-8,-9,-11,-12,-13,-14,-15,-16,-17,-19,-20,-21,-22,
                -23,-24,-25,-26,-29,-30,-31,-32,-33,-34,-35,-37,-39,-40,-42,
                -44,-45,-46,-48)]
anura$pi <- as.factor(anura$pi)
set.seed(3634)

rf_A <- randomForest(pi~., data=df_A, proximity=TRUE, importance=TRUE, ntree=1000, nPerm=10) 
rf_A
print(rf_A)
anu.avg <- 9.23 + 9.74 + 9.77
anu.avg/3 #[1] 9.58

importance(rf_A)

importance(rf_A, type = 2)
round(importance(rf_A))
round(importance(rf_A, 1))
varImpPlot(rf_A)


#####Caudata random forests######
names(caudata) <- c("Species_name","pi","Microhabitat", 
                    "Leaves","Flowers","Seeds","Fruits","Arthro","Vert", 
                    "Activity","Wet_warm","Wet_cold","Dry_warm",              
                    "Dry_cold","Body_mass_g","Age_at_maturity_min_y",
                    "Age_at_maturity_max_y","Body_size_mm",
                    "Size_at_maturity_min_mm","Size_at_maturity_max_mm",
                    "Longevity_max_y","Litter_size_min_n","Litter_size_max_n",
                    "Reproductive_output_y","Offspring_size_min_mm",          
                    "Offspring_size_max_mm", "Development_Mode","Seqs",
                    "Max_Species_SVL_mm","Max_MaleSVL_mm","Max_FemaleSVL_mm",
                    "Vertebrate_Recorded_In_Diet.","Invertebrate_Recorded_In_Diet.", 
                    "PlantMaterial_Recorded_In_Diet.","Elevation_mean",     
                    "Elevation_min","Elevation_max","Lat_mean","Lat_realmax",
                    "Lat_realmin","Area_km2","Order","Family","Genus","Species",
                    "Scientific_name","Biogeographic_Unit","OBS","IUCN")
df_C <- caudata[c(-1,-4,-5,-6,-7,-8,-9,-11,-12,-13,-14,-15,-16,-17,-19,-20,-21,-22, -23,-24,-25,-26,-27, -29,-30,-31,-32,-33,-34,-35,-37,-39,-40,-42,-43, -44,-45,-46,-48)]
caudata$pi <- as.factor(caudata$pi)
set.seed(9078)

rf_C <- randomForest(pi~., data=df_C, proximity=TRUE, importance=TRUE, ntree=1000, nPerm=10) 
print(rf_C)

cau.avg <- 25.15 + 25.33 + 25.93
cau.avg/3 #[1] 25.47

importance(rf_C)

importance(rf_C, type = 2)
round(importance(rf_C))
round(importance(rf_C, 1))
varImpPlot(rf_C)


###Phylogenetic approach

library(ape)
library(phylobase)
library(ggplot2)
library(ggtree)
library(geiger)
library(phytools)

amph.tree <- read.tree(file = "amph_shl_new_Consensus_7238.tre")
class(amph.tree) #phylo
plot(amph.tree)

amph.tree$tip.label

tips <- c("Alsodes_barrioi", 
            "Alsodes_coppingeri", 
            "Alsodes_gargola", 
            "Alsodes_hugoi", 
            "Alsodes_nodosus", 
            "Alsodes_pehuenche", 
            "Alsodes_tumultuosus", 
            "Eupsophus_altor", 
            "Eupsophus_calcaratus", 
            "Eupsophus_contulmoensis", 
            "Eupsophus_insularis", 
            "Eupsophus_migueli", 
            "Eupsophus_nahuelbutensis", 
            "Eupsophus_roseus", 
            "Eupsophus_septentrionalis", 
            "Allobates_brunneus", 
            "Allobates_caeruleodactylus", 
            "Allobates_femoralis", 
            "Allobates_gasconi", 
            "Allobates_hodli", 
            "Allobates_talamancae", 
            "Allobates_tinae", 
            "Allobates_trilineatus", 
            "Anomaloglossus_baeobatrachus", 
            "Anomaloglossus_degranvillei", 
            "Anomaloglossus_roraima", 
            "Mannophryne_olmonae", 
            "Mannophryne_trinitatis", 
            "Atelognathus_nitoi", 
            "Batrachyla_leptopus", 
            "Batrachyla_taeniata", 
            "Brachycephalus_brunneus", 
            "Brachycephalus_crispus", 
            "Brachycephalus_didactylus", 
            "Brachycephalus_ephippium", 
            "Brachycephalus_garbeana", 
            "Brachycephalus_guarani", 
            "Brachycephalus_hermogenesi", 
            "Brachycephalus_izecksohni", 
            "Brachycephalus_margaritatus", 
            "Brachycephalus_nodoterga", 
            "Brachycephalus_pitanga", 
            "Brachycephalus_sulfuratus", 
            "Ischnocnema_guentheri", 
            "Ischnocnema_henselii", 
            "Ischnocnema_parva", 
            "Atelopus_barbotini", 
            "Atelopus_flavescens", 
            "Atelopus_franciscus", 
            "Atelopus_glyphus", 
            "Atelopus_hoogmoedi", 
            "Atelopus_limosus", 
            "Atelopus_senex", 
            "Atelopus_varius", 
            "Atelopus_zeteki", 
            "Incilius_campbelli", 
            "Incilius_coccifer", 
            "Incilius_macrocristatus", 
            "Incilius_nebulifer", 
            "Incilius_occidentalis", 
            "Incilius_valliceps", 
            "Melanophryniscus_alipioi", 
            "Melanophryniscus_biancae", 
            "Melanophryniscus_milanoi", 
            "Melanophryniscus_pachyrhynus", 
            "Melanophryniscus_rubriventris", 
            "Melanophryniscus_stelzneri", 
            "Melanophryniscus_xanthostomus", 
            "Peltophryne_empusa", 
            "Peltophryne_guentheri", 
            "Peltophryne_peltocephala", 
            "Rhinella_arenarum", 
            "Rhinella_diptycha", 
            "Rhinella_dorbignyi", 
            "Rhinella_granulosa", 
            "Rhinella_horribilis", 
            "Rhinella_icterica", 
            "Rhinella_major", 
            "Rhinella_marina", 
            "Craugastor_augusti", 
            "Craugastor_crassidigitus", 
            "Craugastor_fitzingeri", 
            "Craugastor_occidentalis", 
            "Craugastor_raniformis", 
            "Craugastor_talamancae", 
            "Cycloramphus_acangatan", 
            "Cycloramphus_bolitoglossus", 
            "Cycloramphus_boraceiensis", 
            "Cycloramphus_dubius", 
            "Cycloramphus_eleutherodactylus", 
            "Adelphobates_quinquevittatus", 
            "Ameerega_bassleri", 
            "Ameerega_hahneli", 
            "Ameerega_pepperi", 
            "Ameerega_petersi", 
            "Ameerega_pongoensis", 
            "Ameerega_trivittata", 
            "Ameerega_yoshina", 
            "Andinobates_cassidyhornae", 
            "Andinobates_claudiae", 
            "Andinobates_fulguritus", 
            "Andinobates_opisthomelas", 
            "Andinobates_victimatus", 
            "Andinobates_virolinensis", 
            "Colostethus_brachistriatus", 
            "Colostethus_fraterdanieli", 
            "Dendrobates_auratus", 
            "Dendrobates_leucomelas", 
            "Dendrobates_tinctorius", 
            "Epipedobates_anthonyi", 
            "Epipedobates_boulengeri", 
            "Epipedobates_machalilla", 
            "Epipedobates_tricolor", 
            "Hyloxalus_bocagei", 
            "Hyloxalus_chlorocraspedus", 
            "Hyloxalus_elachyhistus", 
            "Hyloxalus_idiomelus", 
            "Hyloxalus_nexipus", 
            "Hyloxalus_subpunctatus", 
            "Leucostethus_jota", 
            "Oophaga_granulifera", 
            "Oophaga_histrionica", 
            "Oophaga_pumilio", 
            "Phyllobates_aurotaenia", 
            "Phyllobates_bicolor", 
            "Phyllobates_lugubris", 
            "Phyllobates_terribilis", 
            "Phyllobates_vittatus", 
            "Ranitomeya_amazonica", 
            "Ranitomeya_fantastica", 
            "Ranitomeya_imitator", 
            "Ranitomeya_reticulata", 
            "Ranitomeya_sirensis", 
            "Ranitomeya_uakarii", 
            "Ranitomeya_variabilis", 
            "Ranitomeya_ventrimaculata", 
            "Silverstoneia_nubicola", 
            "Eleutherodactylus_atkinsi", 
            "Eleutherodactylus_auriculatus", 
            "Eleutherodactylus_coqui", 
            "Eleutherodactylus_cubanus", 
            "Eleutherodactylus_cuneatus", 
            "Eleutherodactylus_dimidiatus", 
            "Eleutherodactylus_eileenae", 
            "Eleutherodactylus_etheridgei", 
            "Eleutherodactylus_feichtingeri", 
            "Eleutherodactylus_glamyrus", 
            "Eleutherodactylus_guanahacabibes", 
            "Eleutherodactylus_iberia", 
            "Eleutherodactylus_jaumei", 
            "Eleutherodactylus_limbatus", 
            "Eleutherodactylus_orientalis", 
            "Eleutherodactylus_planirostris", 
            "Eleutherodactylus_portoricensis", 
            "Eleutherodactylus_rogersi", 
            "Eleutherodactylus_simulans", 
            "Eleutherodactylus_tonyi", 
            "Eleutherodactylus_varleyi", 
            "Bokermannohyla_saxicola", 
            "Dendropsophus_branneri", 
            "Dendropsophus_decipiens", 
            "Dendropsophus_labialis", 
            "Dendropsophus_minutus", 
            "Hyla_eximia", 
            "Hypsiboas_aguilari", 
            "Hypsiboas_albomarginatus", 
            "Hypsiboas_balzani", 
            "Hypsiboas_callipleura", 
            "Hypsiboas_gladiator", 
            "Hypsiboas_marianitae", 
            "Hypsiboas_pulchellus", 
            "Hypsiboas_riojanus", 
            "Lysapsus_boliviana", 
            "Lysapsus_limellum", 
            "Osteocephalus_buckleyi", 
            "Osteocephalus_taurinus", 
            "Osteopilus_septentrionalis", 
            "Phytotriades_auratus", 
            "Plectrohyla_bistincta", 
            "Pseudis_fusca", 
            "Pseudis_paradoxa", 
            "Pseudis_tocantins", 
            "Scinax_cruentommus", 
            "Scinax_granulatus", 
            "Scinax_nebulosus", 
            "Scinax_ruber", 
            "Scinax_squalirostris", 
            "Sphaenorhynchus_caramaschii", 
            "Hylodes_asper", 
            "Adenomera_ajurauna", 
            "Adenomera_andreae", 
            "Adenomera_chicomendesi", 
            "Adenomera_diptyx", 
            "Adenomera_heyeri", 
            "Adenomera_hylaedactyla", 
            "Adenomera_lutzi", 
            "Adenomera_marmorata", 
            "Adenomera_saci", 
            "Adenomera_simonstuarti", 
            "Adenomera_thomei", 
            "Crossodactylodes_bokermanni", 
            "Crossodactylodes_izecksohni", 
            "Engystomops_pustulosus", 
            "Leptodactylus_albilabris", 
            "Pleurodema_cinereum", 
            "Pleurodema_nebulosum", 
            "Proceratophrys_boiei", 
            "Proceratophrys_melanopogon", 
            "Phyllomedusa_ayeaye", 
            "Phyllomedusa_bicolor", 
            "Phyllomedusa_megacephala", 
            "Phyllomedusa_oreades", 
            "Phyllomedusa_rohdei", 
            "Rana_magnaocularis", 
            "Rana_yavapaiensis", 
            "Insuetophrynus_acarpicus", 
            "Pristimantis_achatinus", 
            "Pristimantis_altamazonicus", 
            "Pristimantis_altamnis", 
            "Pristimantis_ardilae", 
            "Pristimantis_croceoinguinis", 
            "Pristimantis_kichwarum", 
            "Pristimantis_luscombei", 
            "Pristimantis_w-nigrum", 
            "Telmatobius_bolivianus", 
            "Telmatobius_chusmisensis", 
            "Telmatobius_fronteriensis", 
            "Telmatobius_halli", 
            "Telmatobius_marmoratus", 
            "Telmatobius_pefauri", 
            "Bolitoglossa_awajun", 
            "Bolitoglossa_celaque", 
            "Bolitoglossa_copinhorum", 
            "Bolitoglossa_gomezi", 
            "Bolitoglossa_lincolni", 
            "Bolitoglossa_mexicana", 
            "Bolitoglossa_morio", 
            "Bolitoglossa_nympha", 
            "Bolitoglossa_occidentalis", 
            "Bolitoglossa_pesrubra", 
            "Bolitoglossa_rufescens", 
            "Bolitoglossa_subpalmata", 
            "Bolitoglossa_yariguiensis", 
            "Bradytriton_silus", 
            "Cryptotriton_veraepacis", 
            "Nototriton_abscondens", 
            "Nototriton_limnospectator", 
            "Oedipina_poelzi", 
            "Oedipina_uniformis", 
            "Pseudoeurycea_bellii", 
            "Pseudoeurycea_leprosa", 
            "Pseudoeurycea_lineola", 
            "Pseudoeurycea_robertsi", 
            "Pseudoeurycea_werleri", 
            "Thorius_maxillabrochus", 
            "Thorius_narisovalis")

my.tree <- phylo4(amph.tree)
class(my.tree)
plot(subset(my.tree, tips.include=tips))
sub.tree <- subset(my.tree, tips.include=tips)
class(sub.tree)
new.tree <- as(sub.tree, "phylo")
class(new.tree)
old.tree <- as.character(new.tree$tip.label)
class(old.tree)
write.tree(new.tree, file = "subset_neotropical.tree")
x <- fastBM(neot.tree)
class(x)

my.pi <- read.csv(file = "ordenar_data.txt", header=FALSE, sep = "", row.names = 1)

neot.tree <- read.newick(file = "neotropical_tree.newick")
class(neot.tree)





## ContMap (phytools)

## extract character of interest
Nucleotide.diversity <- setNames(my.pi$V2,
                          rownames(my.pi))
## create "contMap" object
AmphNeot.contMap<-contMap(neot.tree,
                       Nucleotide.diversity,plot=FALSE,res=200)
## change color scheme
AmphNeot.contMap<-setMap(AmphNeot.contMap,
                       c("#0000FF","#1E90FF", "#87CEEB",	"#F0FFFF","#FFA07A","#FF8C00","#FF4500","#FF0000"))
AmphNeot.contMap<-setMap(AmphNeot.contMap,
                         c("#00BFFF","#FF4500"))
plot(AmphNeot.contMap,fsize=c(0.22,0.8),
     leg.txt="Nucleotide diversity",lwd=2.2, outline=FALSE, sig=2)
par(mar=c(4.1,4.1,4.1,2.1)) ## reset margins to default


phylosig(neot.tree, Nucleotide.diversity, nsim = 999, method = 'lambda')
#Phylogenetic signal K : 0.109084 
#Phylogenetic signal lambda : 0.614976 
#logL(lambda) : 490.991



