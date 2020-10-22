# -------------------------------------------------------------------------------
# Using NMDS of moths as response variable for linear models 
# 22 Oct 2020
# AAR
# -------------------------------------------------------------------------------
# 


# Loading libraries -------------------------------------------------------

library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmulti)
library(MASS)
library(sjPlot)
library(arm)
library(AICcmodavg)
library(glmmTMB)


# Geometridae NMDS --------------------------------------------------------

# Extracting the NMDS axis coordinates for the Geometridae species matrix,
# to be used as response variable

gmatrix <- read.csv("Geo_full_matrix.csv")
code <- read.csv("Geo_NMDScode.csv")

set.seed(15)
GeoOrd <- metaMDS(gmatrix,distance = "bray", k = 2,trymax=100)
summary(GeoOrd)
GeoOrd
#Extract NMDS axis scores
nms_axis <- scores(PlantsOrd, choices=c(1,2))
nms_axis <- as.data.frame(nms_axis)

