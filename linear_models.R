# -------------------------------------------------------------------------------
# Linear models for CR moth paper 2 
# 13 Aug 2020
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


# loading data ------------------------------------------------------------

data_all <- read.csv('data_all.csv')
data_all <- as.data.frame(data_all)


# Notes -------------------------------------------------------------------

# Geometridae - gamma distribution, use glmer
# or I could transform the Geometridae data to be able to use lmer instead of glmer
# Arctiinae - normal distribution, use lmer

# non-correlated floristic variables: Vegetation Diversity, NMDS 1 and NMDS 2
# non-correlated structural variables: Understory complexity, vertical complexity and canopy cover


# -------------------------------------------------------------------------

# First, let's rescale the predictor variables, so they are all in the same scale
data_all$VegDiversity = rescale(data_all$VegDiversity)
data_all$NMDS1 = rescale(data_all$NMDS1)
data_all$NMDS2 = rescale(data_all$NMDS2)
data_all$UnderComplex = rescale(data_all$UnderComplex)
data_all$CanopyCover = rescale(data_all$CanopyCover)
data_all$VerticalComplex = rescale(data_all$VerticalComplex)


# Individual models (without glmulti) -------------------------------------------------------

# Floristic models ------------------------------------------------------------------------

# GEOMETRIDAE floristic models --------------------------------------------------------------

gf.null <- glmer(G_fisher ~ 1+(1|Moonlight)+(1|Habitat),
                 data = data_all, family=Gamma(link = inverse)); summary(gf.null)

gf1 <- glmer(G_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gf1)  # failed to converge with inverse

gf2 <- glmer(G_fisher ~ VegDiversity+NMDS1+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf2)  

gf3 <- glmer(G_fisher ~ VegDiversity+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf3)

gf4 <- glmer(G_fisher ~ NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf4) 

gf5 <- glmer(G_fisher ~ VegDiversity+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf5)

gf6 <- glmer(G_fisher ~ NMDS1+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf6) 

gf7 <- glmer(G_fisher ~ NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf7)

# gf8 <- glmer(G_fisher ~ (VegDiversity*NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
#             data = data_all, family=Gamma(link = log)); summary(gf8) # failed to converge

# gf9 <- glmer(G_fisher ~ (VegDiversity*NMDS1)+(1|Moonlight)+(1|Habitat),
#             data = data_all, family=Gamma(link = log)); summary(gf9)

# gf10 <- glmer(G_fisher ~ (VegDiversity*NMDS2)+(1|Moonlight)+(1|Habitat),
#              data = data_all, family=Gamma(link = log)); summary(gf10)

# gf11 <- glmer(G_fisher ~ (NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
#              data = data_all, family=Gamma(link = log)); summary(gf11)

ggplot(data_all,aes(x=Habitat,y=G_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=G_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7)

# According to this, the model that includes (1+VegDiversity+NMDS2) has the lowest AIC , 
# but there is another model that is close (1+VegDiversity)... 
# Need to compare with another method gf3 vs gf5
# Using glmulti, the best model was (1 + NMDS1 + NMDS2)


# here is another method of comparing models: https://www.scribbr.com/statistics/akaike-information-criterion/

models <- list(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7)
model.names <- c("gf.null","gf1","gf2","gf3","gf4","gf5","gf6","gf7")
aictab(cand.set = models, modnames = model.names)




# ARCTIINAE floristic models --------------------------------------------------------------
# Had to remove Moonlight as a random effect, because the variance was too close to cero (isSingular error)

af.null <- lmer(A_fisher ~ 1+(1|Habitat),
                data = data_all); summary(af.null)  

af1 <- lmer(A_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Habitat),
            data = data_all); summary(af1)  

af2 <- lmer(A_fisher ~ VegDiversity+NMDS1+(1|Habitat),
            data = data_all); summary(af2) 

af3 <- lmer(A_fisher ~ VegDiversity+NMDS2+(1|Habitat),
            data = data_all); summary(af3)  

af4 <- lmer(A_fisher ~ NMDS1+NMDS2+(1|Habitat),
            data = data_all); summary(af4)  

af5 <- lmer(A_fisher ~ VegDiversity+(1|Habitat),
            data = data_all); summary(af5)

af6 <- lmer(A_fisher ~ NMDS1+(1|Habitat),
            data = data_all); summary(af6) 

af7 <- lmer(A_fisher ~ NMDS2+(1|Habitat),
            data = data_all); summary(af7) 

# af8 <- lmer(A_fisher ~ (VegDiversity*NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
#            data = data_all); summary(af8)  # is Singular

# af9 <- lmer(A_fisher ~ (VegDiversity*NMDS1)+(1|Moonlight)+(1|Habitat),
#            data = data_all); summary(af9)  

# af10 <- lmer(A_fisher ~ (VegDiversity*NMDS2)+(1|Moonlight)+(1|Habitat),
#             data = data_all); summary(af10)  # is Singular

# af11 <- lmer(A_fisher ~ (NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
#             data = data_all); summary(af11)  # is Singular

ggplot(data_all,aes(x=Habitat,y=A_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=A_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(af.null,af1,af2,af3,af4,af5,af6,af7, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(af.null,af1,af2,af3,af4,af5,af6,af7)

# According to this, the model that includes (1 + VegDiversity+NMDS1+NMDS2) has the lowest AIC
# But, there is another model that is very close (1 + VegDiversity+NMDS1). Need to check which is best.
# Using glmulti, the best model was also (1 + NMDS1)

# si la varianza de Moonlight es cerca de 0 en todos los modelos de Arctiinae, tal vez puedo
# remover Moonlight de todos los modelos y justificarlo bien. 
# Tal vez la Luna afecta las Arctiinae mucho menos que a las Geometridae... Revisar esto

# Structural model -------------------------------------------------------------------------

# GEOMETRIDAE structural models --------------------------------------------------------------

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
             data = data_all, family=Gamma(link = log)); summary(gs5) # does not converge

gs6 <- glmer(G_fisher ~ CanopyCover+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs6) # does not converge

gs7 <- glmer(G_fisher ~ VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs7)

gs8 <- glmer(G_fisher ~ (UnderComplex * CanopyCover * VerticalComplex) +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs8) # does not converge

gs9 <- glmer(G_fisher ~ (UnderComplex * CanopyCover) +(1|Moonlight)+(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs9) # does not converge

gs10 <- glmer(G_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
              data = data_all, family=Gamma(link = log)); summary(gs10)

gs11 <- glmer(G_fisher ~ (CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
              data = data_all, family=Gamma(link = log)); summary(gs11) # does not converge



ggplot(data_all,aes(x=Habitat,y=G_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=G_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(gs.null,gs1,gs2,gs3,gs4,gs5,gs6,gs7,gs8,gs9,gs10,gs11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

# According to this, the null model has the lowest AIC
# With glmulti, the best model included (1 + UnderComplex + VerticalComplex)

# ARCTIINAE

as.null <- lmer(A_fisher ~ 1+(1|Moonlight)+(1|Habitat),
                data = data_all); summary(as.null) # is Singular

as1 <- lmer(A_fisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as1)  

as2 <- lmer(A_fisher ~ UnderComplex + CanopyCover+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as2)  

as3 <- lmer(A_fisher ~ UnderComplex+VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as3) 

as4 <- lmer(A_fisher ~ CanopyCover + VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as4) 

as5 <- lmer(A_fisher ~ UnderComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as5)  # is Singular

as6 <- lmer(A_fisher ~ CanopyCover+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as6)  

as7 <- lmer(A_fisher ~ VerticalComplex+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as7)  

as8 <- lmer(A_fisher ~ (UnderComplex*CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as8)  # in Singular 

as9 <- lmer(A_fisher ~ (UnderComplex*CanopyCover)+(1|Moonlight)+(1|Habitat),
            data = data_all); summary(as9)  

as10 <- lmer(A_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(as10) 

as11 <- lmer(A_fisher ~ (CanopyCover*VerticalComplex)+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(as11) 


tab_model(as.null,as1,as2,as3,as4,as5,as6,as7,as8,as9,as10,as11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(a.null,a1,a2,a3,a4,a5,a6,a7)

# According to this, the model that includes (1+UnderComplex*VerticalComplex) has the lowest AIC
# With glmulti, the best model was (1+UnderComplex+VerticalComplex), but I have
# yet to conduct interactions (*) using glmulti. 



# -------------------------------------------------------------------------

# Geometridae models using log transformed G_fisher. ----------------------
# (lmer instead of glmer with gamma dist)

# Floristic

gft.null <- lmer(logGfisher ~ 1+(1|Moonlight)+(1|Habitat),
                 data = data_all); summary(gft.null)  # is Singular

gft1 <- lmer(logGfisher ~ VegDiversity+NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft1)  # is Singular

gft2 <- lmer(logGfisher ~ VegDiversity+NMDS1+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft2) 

gft3 <- lmer(logGfisher ~ VegDiversity+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft3) 

gft4 <- lmer(logGfisher ~ NMDS1+NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft4)  # is Singular

gft5 <- lmer(logGfisher ~ VegDiversity+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft5)

gft6 <- lmer(logGfisher ~ NMDS1+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft6)  # is Singular

gft7 <- lmer(logGfisher ~ NMDS2+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft7)  

gft8 <- lmer(logGfisher ~ (VegDiversity*NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft8)  

gft9 <- lmer(logGfisher ~ (VegDiversity*NMDS1)+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gft9)  

gft10 <- lmer(logGfisher ~ (VegDiversity*NMDS2)+(1|Moonlight)+(1|Habitat),
              data = data_all); summary(gft10)  

gft11 <- lmer(logGfisher ~ (NMDS1*NMDS2)+(1|Moonlight)+(1|Habitat),
              data = data_all); summary(gft11) 

tab_model(gft.null,gft1,gft2,gft3,gft4,gft5,gft6,gft7,gft8,gft9,gft10,gft11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

# According to this, the null model has the lowest AIC
# Using glmulti, the best model was (1 + NMDS1)

# Structural

gst.null <- lmer(logGfisher ~ 1+(1|Moonlight)+(1|Habitat),
                 data = data_all); summary(gst.null) 

gst1 <- lmer(logGfisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst1)  # is Singular

gst2 <- lmer(logGfisher ~ UnderComplex + CanopyCover+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst2)  # is Singular

gst3 <- lmer(logGfisher ~ UnderComplex+VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst3)  # is Singular

gst4 <- lmer(logGfisher ~ CanopyCover + VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst4) 

gst5 <- lmer(logGfisher ~ UnderComplex+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst5)  # is Singular

gst6 <- lmer(logGfisher ~ CanopyCover+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst6)  

gst7 <- lmer(logGfisher ~ VerticalComplex+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst7)  

gst8 <- lmer(logGfisher ~ (UnderComplex*CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst8)  # in Singular 

gst9 <- lmer(logGfisher ~ (UnderComplex*CanopyCover)+(1|Moonlight)+(1|Habitat),
             data = data_all); summary(gst9)  # does not converge

gst10 <- lmer(logGfisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
              data = data_all); summary(gst10) 

gst11 <- lmer(logGfisher ~ (CanopyCover*VerticalComplex)+(1|Moonlight)+(1|Habitat),
              data = data_all); summary(gst11)  # is Singular


tab_model(gst.null,gst1,gst2,gst3,gst4,gst5,gst6,gst7,gst8,gst9,gst10,gst11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

# According to this, the null model has the lowest AIC
# Using glmulti, the best model was (1 + VerticalComplex)


