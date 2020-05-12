# Data analysis for paper 2 of La Gamba moth project
# Evaluate effects of vegetation parameters on 
# moth diversity in 4 land use types in La Gamba CR

# AAR


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
PlantsOrd <- metaMDS(plants,distance = "bray", k = 3,trymax=100)
summary(PlantsOrd)
#Extract NMDS axis scores
nms_axis <- scores(PlantsOrd, choices=c(1,2))
nms_axis <- as.data.frame(nms_axis)
# write.csv(nms_axis,"C:\\Users\\aalonsor\\OneDrive - University of Vermont\\Documents\\UVM\\Classes\\Spring 2020\\Computational Biology\\AlonsoRodzBio381\\nms_axis.csv", row.names = TRUE)


# Creating data frame with all variables ----------------------------------

covar <- read.csv("covar.csv")
data_all <- cbind(covar, nms_axis, G_fisher, A_fisher)


# Checking distribution of data -------------------------------------------


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

# transform G_fisher
data_all$logGfisher <- log(data_all$G_fisher)
shapiro.test(data_all$logGfisher)
descdist(data_all$logGfisher, discrete=FALSE, boot=500) 

# Checking for collinearity of variables ----------------------------------

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

# QUESTIONS: 
# ADD SITE AS RANDOM EFFECT???
# HOW TO CONDUCT A SPATIAL AUTOCORRELATION TEST WITH LAT LONG OF SITES?
# OKAY TO REMOVE BASAL AREA FROM ANALYSES?
# HOW TO CHECK FOR MODEL VALIDATION



# Model with all variables that are not correlated ------------------------

# According to the correlation plot that includes all variables, there were 5 
# variables that were selected as not being correlated:
# Vegetation Diversity, Understory Complexity, Canopy Cover, NMDS 1 and NMDS 2

# Functions
lmer.glmulti=function(formula, data, random = "",...) { ### glmulti
  lmer(paste(deparse(formula), random),data=data,...)
}

glmer.glmulti=function(formula, data, random = "",...) {
  glmer(paste(deparse(formula), random),data=data, REML=F,...)
}

# GEOMETRIDAE

Geom_all <- glmulti(G_fisher ~ VegDiversity+UnderComplex+CanopyCover+NMDS1+NMDS2,
                    level=1, fitfunc=glmer.glmulti, random=c("+(1|Moonlight)"), 
                    data=data_all, method ="h", crit = "aicc")

Geom_all <- glmulti(G_fisher ~ VegDiversity+UnderComplex+CanopyCover+NMDS1+NMDS2,
                  level=1, fitfunc=glmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                  data=data_all, method ="h", crit = "aicc")
print(Geom_all)
AIC <- weightable(Geom_all)
AIC[1:7,]
plot(Geom_all, type="s")

Geom_final <- glmer(G_fisher ~1+NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(Geom_final)

# ARCTIINAE

Arc_all <- glmulti(A_fisher ~ VegDiversity+UnderComplex+CanopyCover+NMDS1+NMDS2,
                    level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                   data=data_all, method = "h", crit = "aicc")
print(Arc_all)
AIC <- weightable(Arc_all)
AIC[1:2,]
plot(Arc_all, type="s")


Arc_final <- lmer(A_fisher ~1+UnderComplex+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
                    data = data_all); summary(Arc_final)



# Floristic model ---------------------------------------------------------

# According to the correlation plot that includes floristic variables, there were 3
# variables that were selected as not being correlated:
# Vegetation Diversity, NMDS 1 and NMDS 2

# GEOMETRIDAE

Geom_flor <- glmulti(G_fisher ~ VegDiversity+NMDS1+NMDS2,
                    level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                    data=data_all, method ="h", crit = "aicc")
# ERROR: function above does not work with glmer.glmulti. WHY?

print(Geom_flor)
AIC <- weightable(Geom_flor)
AIC[1:7,]
plot(Geom_flor, type="s")

Geom_flor_final <- glmer(G_fisher ~1+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
                    data = data_all, family=Gamma(link = log)); summary(Geom_flor_final)

# ARCTIINAE

Arc_flor <- glmulti(A_fisher ~ VegDiversity+NMDS1+NMDS2,
                   level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                   data=data_all, method = "h", crit = "aicc")
print(Arc_flor)
AIC <- weightable(Arc_flor)
AIC[1:2,]
plot(Arc_flor, type="s")


Arc_flor_final <- lmer(A_fisher ~1+NMDS1+(1|Moonlight)+(1|Habitat),
                  data = data_all); summary(Arc_flor_final)


# Structural Model --------------------------------------------------------

# According to the correlation plot that includes structural variables, there were 4
# variables that were selected as not being correlated:
# Understory complexity, vertical complexity, canopy cover and basal area
# Decided to remove basal area from model, because it is severely high in oil palm

# GEOMETRIDAE

Geom_str <- glmulti(G_fisher ~ UnderComplex + CanopyCover + VerticalComplex,
                     level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                     data=data_all, method ="h", crit = "aicc")
# ERROR: function above does not work with glmer.glmulti. WHY?

print(Geom_str)
AIC <- weightable(Geom_str)
AIC[1:7,]
plot(Geom_str, type="s")

Geom_str_final <- glmer(G_fisher ~1+UnderComplex + VerticalComplex+(1|Moonlight)+(1|Habitat),
                         data = data_all, family=Gamma(link = log)); summary(Geom_str_final)

# ARCTIINAE

Arc_str <- glmulti(A_fisher ~ UnderComplex + CanopyCover + VerticalComplex,
                    level=1, fitfunc=lmer.glmulti, random=c("+(1|Moonlight)","+(1|Habitat)"), 
                    data=data_all, method = "h", crit = "aicc")
print(Arc_str)
AIC <- weightable(Arc_str)
AIC[1:3,]
plot(Arc_str, type="s")


Arc_str_final <- lmer(A_fisher ~1+UnderComplex + VerticalComplex+(1|Moonlight)+(1|Habitat),
                       data = data_all); summary(Arc_str_final)

#------------------------------------------------------------------------------------------
#NOTES

# tab_model de paquete sjplots para poner todos los modelos individuales
# en una misma tabla con sus aicc

# nested vs crossed random effects -- habitat/site o site/habitat... como ponerlos
# ben bolker - buscar su pagina


#------------------------------------------------------------------------------------------
# Running models individually, without glmulti --------------------------------------------

plot(G_fisher~Moonlight, data = data_all)
scatter.smooth(x=data_all$Moonlight, y=data_all$G_fisher)

plot(G_fisher~Habitat, data = data_all)

plot(G_fisher~Code, data = data_all)


# Floristic model ------------------------------------------------------------------------

# GEOMETRIDAE

gf.null <- glmer(G_fisher ~ 1+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf.null)

gf1 <- glmer(G_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf1)

gf2 <- glmer(G_fisher ~ VegDiversity+NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf2)

gf3 <- glmer(G_fisher ~ VegDiversity+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf3)

gf4 <- glmer(G_fisher ~ NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf4) # is Singular

gf5 <- glmer(G_fisher ~ VegDiversity+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf5)

gf6 <- glmer(G_fisher ~ NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf6) # is Singular

gf7 <- glmer(G_fisher ~ NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf7)

gf8 <- glmer(G_fisher ~ (VegDiversity*NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf8)

gf9 <- glmer(G_fisher ~ (VegDiversity*NMDS1)+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf9)

gf10 <- glmer(G_fisher ~ (VegDiversity*NMDS2)+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf10)

gf11 <- glmer(G_fisher ~ (NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
            data = data_all, family=Gamma(link = log)); summary(gf11)

ggplot(data_all,aes(x=Habitat,y=G_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=G_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7,gf8,gf9,gf10,gf11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)


# ARCTIINAE


af.null <- lmer(A_fisher ~ 1+(1|Moonlight)+(1|Habitat),
                data = data_all); summary(af.null)  # is Singular

af1 <- lmer(A_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af1)  # is Singular

af2 <- lmer(A_fisher ~ VegDiversity+NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af2)  # is Singular

af3 <- lmer(A_fisher ~ VegDiversity+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af3)  # is Singular

af4 <- lmer(A_fisher ~ NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af4)  # is Singular

af5 <- lmer(A_fisher ~ VegDiversity+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af5)

af6 <- lmer(A_fisher ~ NMDS1+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af6)  # is Singular

af7 <- lmer(A_fisher ~ NMDS2+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(af7)  # failed to converge

af8 <- lmer(A_fisher ~ (VegDiversity*NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
           data = data_all); summary(af8)  # is Singular

af9 <- lmer(A_fisher ~ (VegDiversity*NMDS1)+(1|Moonlight)+(1|Habitat),
           data = data_all); summary(af9)  

af10 <- lmer(A_fisher ~ (VegDiversity*NMDS2)+(1|Moonlight)+(1|Habitat),
           data = data_all); summary(af10)  # is Singular

af11 <- lmer(A_fisher ~ (NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
           data = data_all); summary(af11)  # is Singular

tab_model(af.null,af1,af2,af3,af4,af5,af6,af7,af8,af9,af10,af11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(a.null,a1,a2,a3,a4,a5,a6,a7)


# Structural model -------------------------------------------------------------------------

# GEOMETRIDAE

gs.null <- glmer(G_fisher ~ 1+(1|Moonlight)+(1|Habitat),
                 data = data_all, family=Gamma(link = log)); summary(gs.null)

gs1 <- glmer(G_fisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs1)

gs2 <- glmer(G_fisher ~ UnderComplex + CanopyCover +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs2)

gs3 <- glmer(G_fisher ~ UnderComplex+VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs3)

gs4 <- glmer(G_fisher ~ CanopyCover + VerticalComplex +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs4)

gs5 <- glmer(G_fisher ~ UnderComplex+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs5)

gs6 <- glmer(G_fisher ~ CanopyCover+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs6) 

gs7 <- glmer(G_fisher ~ VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs7)

gs8 <- glmer(G_fisher ~ (UnderComplex * CanopyCover * VerticalComplex) +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs8)

gs9 <- glmer(G_fisher ~ (UnderComplex * CanopyCover) +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs9)

gs10 <- glmer(G_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs10)

gs11 <- glmer(G_fisher ~ (CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs11)



ggplot(data_all,aes(x=Habitat,y=G_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=G_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(gs.null,gs1,gs2,gs3,gs4,gs5,gs6,gs7,gs8,gs9,gs10,gs11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)


# ARCTIINAE

as.null <- lmer(A_fisher ~ 1+(1|Moonlight)+(1|Habitat),
                data = data_all); summary(as.null) 

as1 <- lmer(A_fisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as1)  

as2 <- lmer(A_fisher ~ UnderComplex + CanopyCover+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as2)  

as3 <- lmer(A_fisher ~ UnderComplex+VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as3) 

as4 <- lmer(A_fisher ~ CanopyCover + VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as4) 

as5 <- lmer(A_fisher ~ UnderComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as5)

as6 <- lmer(A_fisher ~ CanopyCover+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as6)  

as7 <- lmer(A_fisher ~ VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as7)  

as8 <- lmer(A_fisher ~ (UnderComplex*CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as8)  

as9 <- lmer(A_fisher ~ (UnderComplex*CanopyCover)+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as9)  

as10 <- lmer(A_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as10) 

as11 <- lmer(A_fisher ~ (CanopyCover*VerticalComplex)+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as11) 


tab_model(as.null,as1,as2,as3,as4,as5,as6,as7,as8,as9,as10,as11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(a.null,a1,a2,a3,a4,a5,a6,a7)





