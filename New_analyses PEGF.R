



# ---------------------------------------------
# New Analysis: 24 de Junio 2022
# 24 Jun 2022
# Pablo E. Gutiérrez-Fonseca
# pabloe.gutierrezfonseca@gmail.com
# ---------------------------------------------
#  


# Loading libraries -------------------------------------------------------

library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmulti)
library(MASS)
library(fitdistrplus)
library(corrplot)
library(sjPlot)
library(arm)
library(bbmle)

# Fisher's alpha of moth data ---------------------------------------------

G_matrix <- read.csv('Geo_site_matrix.csv')
G_fisher <- fisher.alpha(G_matrix) 
G_fisher <- as.data.frame(G_fisher)

A_matrix <- read.csv('Arc_site_matrix.csv')
A_fisher <- fisher.alpha(A_matrix) 
A_fisher <- as.data.frame(A_fisher)


# NMDS of plant data ------------------------------------------------------

# Extracting the NMDS axis coordinates for the plant species matrix,
# to be used as proxy for plant species composition

plants <- read.csv("Plant_matrix.csv")
habitat <- read.csv("Habitat.csv")

set.seed(15)
PlantsOrd <- metaMDS(plants, distance = "bray", k = 2, trymax=100)
summary(PlantsOrd)
PlantsOrd

par(mfrow=c(1,1))
plot(PlantsOrd)

#Extract NMDS axis scores
nms_axis <- as.data.frame(scores(PlantsOrd, 'sites'))  # Pablo 24 de Junio

# Creating data frame with all variables ----------------------------------
covar <- read.csv("covar.csv")
data_all <- cbind(covar, nms_axis, G_fisher, A_fisher)


# Response variables ------------------------------------------------------

shapiro.test(data_all$G_fisher) # not normal
shapiro.test(data_all$A_fisher) # normal

par(mfrow=c(2,2))
descdist(data_all$G_fisher, discrete=FALSE, boot=500) # GAMMA
descdist(data_all$A_fisher, discrete=FALSE, boot=500) # UNIFORM

plot(data_all$G_fisher, as.factor(data_all$Habitat))
plot(data_all$A_fisher, as.factor(data_all$Habitat))



# Before GLM, scaled variables

# non-correlated floristic variables: Vegetation Diversity, NMDS 1 and NMDS 2
# non-correlated structural variables: Understory complexity, vertical complexity and canopy cover

# First, let's rescale the predictor variables, so they are all in the same scale
data_all$VegDiversity = rescale(data_all$VegDiversity)
data_all$NMDS1 = rescale(data_all$NMDS1)
data_all$NMDS2 = rescale(data_all$NMDS2)
data_all$UnderComplex = rescale(data_all$UnderComplex)
data_all$CanopyCover = rescale(data_all$CanopyCover)
data_all$VerticalComplex = rescale(data_all$VerticalComplex)


# GEOMETRIDAE Floristic models --------------------------------------------

gf.null <- glmer(G_fisher ~ 1 + (1|Habitat), 
                 data = data_all, family=Gamma)

gf1 <- glmer(G_fisher ~ VegDiversity + NMDS1 + NMDS2 + (1|Habitat),
             data = data_all, family=Gamma)

gf2 <- glmer(G_fisher ~ VegDiversity + NMDS1 + (1|Habitat),
             data = data_all, family=Gamma)  

gf3 <- glmer(G_fisher ~ VegDiversity + NMDS2 + (1|Habitat),
             data = data_all, family=Gamma)

gf4 <- glmer(G_fisher ~ NMDS1 + NMDS2 +(1|Habitat),
             data = data_all, family=Gamma) # is Singular
ranef(gf4)
isSingular(gf4, tol = 1e-4)

gf5 <- glmer(G_fisher ~ VegDiversity + (1|Habitat),
             data = data_all, family=Gamma)
isSingular(gf5, tol = 1e-4)

gf6 <- glmer(G_fisher ~ NMDS1 + (1|Habitat),
             data = data_all, family=Gamma) # is Singular
ranef(gf6)

gf7 <- glmer(G_fisher ~ NMDS2 + (1|Habitat),
             data = data_all, family=Gamma)

anova(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7,test="F")
AIC(gf1,gf2,gf5)
bbmle::AICctab(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7, base = T,weights = T)
