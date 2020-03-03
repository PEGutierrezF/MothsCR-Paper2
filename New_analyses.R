# Data analysis for paper 2 of La Gamba moth project
# Evaluate effects of vegetation parameters on 
# moth diversity in 4 land use types in La Gamba CR

# AAR

# Fisher's alpha of moth data ---------------------------------------------

library(vegan)

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

#install.packages("MASS")
library(MASS)

PlantsOrd <- metaMDS(plants,distance = "bray", k = 3,trymax=100)
#Extract NMDS axis scores
nms_axis <- scores(PlantsOrd, choices=c(1,2))
nms_axis <- as.data.frame(nms_axis)


# Creating data frame with all variables ----------------------------------

covar <- read.csv("covar.csv")
data_all <- cbind(covar, nms_axis, G_fisher, A_fisher)


# Checking distribution of data -------------------------------------------

library(fitdistrplus)
descdist(data_all$SpRichness, discrete=FALSE, boot=500) # UNIFORM
descdist(data_all$VegDensity, discrete=FALSE, boot=500) # GAMMA
descdist(data_all$VegDiversity, discrete=FALSE, boot=500) # BETA
descdist(data_all$UnderDensity, discrete=FALSE, boot=500) # UNIFORM
descdist(data_all$UnderComplex, discrete=FALSE, boot=500) # BETA, close to NORMAL or UNIFORM
descdist(data_all$UnderCover, discrete=FALSE, boot=500) # BETA, close to UNIFORM
descdist(data_all$VerticalComplex, discrete=FALSE, boot=500) # BETA
descdist(data_all$CanopyCover, discrete=FALSE, boot=500) # BETA
descdist(data_all$BasalArea, discrete=FALSE, boot=500) # BETA
descdist(data_all$Elevation, discrete=FALSE, boot=500) # GAMMA
descdist(data_all$Slope, discrete=FALSE, boot=500) # UNIFORM
descdist(data_all$Moonlight, discrete=FALSE, boot=500) # BETA
descdist(data_all$NMDS1, discrete=FALSE, boot=500) # BETA
descdist(data_all$NMDS2, discrete=FALSE, boot=500) # BETA,close to NORMAL or UNIFORM

descdist(data_all$G_fisher, discrete=FALSE, boot=500) # GAMMA
descdist(data_all$A_fisher, discrete=FALSE, boot=500) # UNIFORM

shapiro.test(data_all$SpRichness) # not normal
shapiro.test(data_all$VegDensity) # not normal
shapiro.test(data_all$VegDiversity) # not normal
shapiro.test(data_all$UnderDensity) # normal
shapiro.test(data_all$UnderComplex) # normal
shapiro.test(data_all$UnderCover) # not normal
shapiro.test(data_all$VerticalComplex) # not normal
shapiro.test(data_all$CanopyCover) # not normal
shapiro.test(data_all$BasalArea) # not normal
shapiro.test(data_all$Elevation) # not normal
shapiro.test(data_all$Slope) # not normal
shapiro.test(data_all$Moonlight) # normal
shapiro.test(data_all$NMDS1) # not normal
shapiro.test(data_all$NMDS2) # normal

shapiro.test(data_all$G_fisher) # not normal
shapiro.test(data_all$A_fisher) # normal

# Checking for collinearity of variables ----------------------------------

library(corrplot)
var_cor=subset(data_all, select = c("SpRichness","VegDensity","VegDiversity",
                             "UnderDensity","UnderComplex","UnderCover",
                             "VerticalComplex","CanopyCover",
                             "BasalArea","NMDS1","NMDS2"))
var_cor <- as.matrix(var_cor)
c <- corrplot(cor(var_cor,use="pairwise.complete.obs", method = "spearman"), 
         is.corr = TRUE, method = "number", tl.cex=1, 
         type = "upper", order = "hclust")

# res1 <- cor.mtest(var_cor, conf.level = 0.95)
# corrplot(cor(var_cor,use="pairwise.complete.obs", method = "spearman"),
#         method = "square", type = "upper", order = "hclust",
#         p.mat = res1$p, sig.level = 0.05)

# More than 0.6 correlation coefficient is bad

str(c)
head(c)
correlations <- as.data.frame(c)

# testing correlations between floristic and structural variables separately

# floristic

var_cor_flor=subset(data_all, select = c("SpRichness","VegDensity","VegDiversity",
                                    "UnderDensity","NMDS1","NMDS2"))
var_cor_flor <- as.matrix(var_cor_flor)
c_flor <- corrplot(cor(var_cor_flor,use="pairwise.complete.obs", method = "spearman"), 
              is.corr = TRUE, method = "number", tl.cex=1, 
              type = "upper", order = "hclust")

# structural

var_cor_struc=subset(data_all, select = c("UnderComplex","UnderCover",
                                    "VerticalComplex","CanopyCover",
                                    "BasalArea"))
var_cor_struc <- as.matrix(var_cor_struc)
corrplot(cor(var_cor_struc,use="pairwise.complete.obs", method = "spearman"), 
              is.corr = TRUE, method = "number", tl.cex=1, 
              type = "upper", order = "hclust")



# NEXT STEPS -----------------------------

# do model with all variables, just to try it even though its not the way to go
# do floristic model and stuctural model separately with all variables, just to try it as well
# do floristic model and structural model separately with pre-selected variables
# perform model selection glmulti for all of these
# determine best final model for arctiines and geometrids

# Geometridae - gamma, use glmer
# Arctiinae - normal, use lmer

# QUESTION: ADD SITE AS RANDOM EFFECT???

# Model with all variables that are not correlated ------------------------

# According to the correlation plot that includes all variables, there were 5 
# variables that were selected as not being correlated:
# Vegetation Diversity, Understory Complexity, Canopy Cover, NMDS 1 and NMDS 2

#install.packages("lme4")
library(lme4)
#install.packages("lmerTest")
library(lmerTest)
#install.packages("glmulti")
library(glmulti)


# Functions
lmer.glmulti=function(formula, data, random = "",...) { ### glmulti
  lmer(paste(deparse(formula), random),data=data,...)
}

glmer.glmulti=function(formula, data, random = "",...) {
  glmer(paste(deparse(formula), random),data=data, REML=F,...)
}

# GEOMETRIDAE

Geom_all <- glmulti(G_fisher ~ VegDiversity+UnderComplex+CanopyCover+NMDS1+NMDS2,
                  level=1, fitfunc=glmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                  data=data_all, method ="h")
plot(Geom_all, type="s")

Geom_final <- glmer(G_fisher ~1+NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(Geom_final)

# ARCTIINAE

Arc_all <- glmulti(A_fisher ~ VegDiversity+UnderComplex+CanopyCover+NMDS1+NMDS2,
                    level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                   data=data_all, method ="h")
plot(Arc_all, type="s")

Arc_final <- lmer(A_fisher ~1+UnderComplex+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
                    data = data_all); summary(Arc_final)


# Floristic model ---------------------------------------------------------


