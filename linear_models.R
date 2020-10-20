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

install.packages("glmmTMB")
library(glmmTMB)

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


# Checking random effects -------------------------------------------------

# Moonlight
plot(data_all$G_fisher, data_all$Moonlight)
plot(data_all$A_fisher, data_all$Moonlight)

cor.test(data_all$G_fisher, data_all$Moonlight, method=c("spearman"))
cor.test(data_all$A_fisher, data_all$Moonlight, method=c("pearson"))

# Moonlight is not significantly affecting moth diversity, so it will no longer be included
# as a random effect in any of the models.

# Habitat
plot(data_all$G_fisher, as.factor(data_all$Habitat))
plot(data_all$A_fisher, as.factor(data_all$Habitat))



# Individual models (without glmulti) -------------------------------------------------------

# Floristic models ------------------------------------------------------------------------

# GEOMETRIDAE floristic models --------------------------------------------------------------

gf.null <- glmer(G_fisher ~ 1+(1|Habitat),
                 data = data_all, family=Gamma(link = inverse)); summary(gf.null)

gf1 <- glmer(G_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf1)  

gf2 <- glmer(G_fisher ~ VegDiversity+NMDS1+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf2)  

gf3 <- glmer(G_fisher ~ VegDiversity+NMDS2+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf3)

gf4 <- glmer(G_fisher ~ NMDS1+NMDS2+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf4) 

gf5 <- glmer(G_fisher ~ VegDiversity+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf5)

gf6 <- glmer(G_fisher ~ NMDS1+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gf6) 

gf7 <- glmer(G_fisher ~ NMDS2+(1|Habitat),
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
        

ggplot(data_all,aes(x=VegDiversity,y=G_fisher)) + geom_jitter() + 
        geom_point(alpha=0.2) + geom_smooth(method = "lm")
plot(gf5)    # para ver los residuales?

tab_model(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7)

# here is another method of comparing models: https://www.scribbr.com/statistics/akaike-information-criterion/

models <- list(gf.null,gf1,gf2,gf3,gf4,gf5,gf6,gf7)
model.names <- c("gf.null","gf1","gf2","gf3","gf4","gf5","gf6","gf7")
aictab(cand.set = models, modnames = model.names)

# According to tab_model and aictab, the model that includes (1+VegDiversity) has the 
# lowest AIC, but the null model is very close with a deltaAIC<2... 
# Need to compare with another method null vs gf5 (e.g. model averaging), or just mention both
# Using glmulti, the best model was (1 + NMDS1 + NMDS2)


# This method below is to compare between the best candidate models
# this is called model averaging
# https://danstich.github.io/stich/classes/BIOL678/06_modelSelection.html

# Cand.mod <- list()
# Cand.mod[[1]]<-gf3
# Cand.mod[[2]]<-gf5
# 
# 
# Modnames <- c("VegDiversity+NMDS2","VegDiversity") 
# aictab(cand.set = Cand.mod, modnames = Modnames)
# evidence(aictab(cand.set = Cand.mod, modnames = Modnames))
# 
# confset(cand.set = Cand.mod, modnames = Modnames, second.ord = TRUE,
#         method = "raw")
# 
# 
# modavg(Cand.mod, "VegDiversity", modnames=Modnames, c.hat = 1, gamdisp = NULL,
#        conf.level = 0.95, second.ord = TRUE, nobs = NULL,
#        exclude = list("VegDiversity+NMDS2"), warn = TRUE, uncond.se = "revised",
#        parm.type = NULL)
# 
# modavg(Cand.mod, "NMDS2", modnames=Modnames, c.hat = 1, gamdisp = NULL,
#        conf.level = 0.95, second.ord = TRUE, nobs = NULL,
#        exclude = list(""), warn = TRUE, uncond.se = "revised",
#        parm.type = NULL)

# From these results, it seems that NMDS2 is not significant (CI passes by 0)
# This might mean that NMDS2 should not be included at all, in which case the 
# gf5 model might be the most parsimonious, which only includes VegDiversity.

# ** another option to try could be using 'model.avg'


# ARCTIINAE floristic models --------------------------------------------------------------
# Had to remove Moonlight as a random effect, because the variance was too close to cero (isSingular error)

af.null <- lmer(A_fisher ~ 1+(1|Habitat),
                data = data_all, REML = FALSE); summary(af.null)  

af1 <- lmer(A_fisher ~ VegDiversity+NMDS1+NMDS2+(1|Habitat),
            data = data_all, REML = FALSE); summary(af1)     # is singular con el REML=FALSE

ggplot(data_all, aes(x=NMDS1, y=A_fisher)) + geom_point()
plot(af1)     # para ver los residuales?
# El problema es que la palma me esta dividiendo mucho los datos

af2 <- lmer(A_fisher ~ VegDiversity+NMDS1+(1|Habitat),
            data = data_all, REML = FALSE); summary(af2) 

af3 <- lmer(A_fisher ~ VegDiversity+NMDS2+(1|Habitat),
            data = data_all, REML = FALSE); summary(af3)  

af4 <- lmer(A_fisher ~ NMDS1+NMDS2+(1|Habitat),
            data = data_all, REML = FALSE); summary(af4)    # is singular con REML = FALSE

af5 <- lmer(A_fisher ~ VegDiversity+(1|Habitat),
            data = data_all, REML = FALSE); summary(af5)

af6 <- lmer(A_fisher ~ NMDS1+(1|Habitat),
            data = data_all, REML = FALSE); summary(af6)    # is singular

af7 <- lmer(A_fisher ~ NMDS2+(1|Habitat),
            data = data_all, REML = FALSE); summary(af7) 

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

models <- list(af.null,af1,af2,af3,af4,af5,af6,af7)
model.names <- c("af.null","af1","af2","af3","af4","af5","af6","af7")
aictab(cand.set = models, modnames = model.names)

# According to this, the model that includes (1 + VegDiversity+NMDS1+NMDS2) has the lowest AIC
# But, there is another model that is very close (1 + VegDiversity+NMDS1). Need to check which is best.
# Using glmulti, the best model was also (1 + NMDS1)

# This method below is to compare between the best candidate models
# this is called model averaging
# https://danstich.github.io/stich/classes/BIOL678/06_modelSelection.html

Cand.mod <- list()
Cand.mod[[1]]<-af1
Cand.mod[[2]]<-af2


Modnames <- c("VegDiversity+NMDS1+NMDS2","VegDiversity+NMDS1") 
aictab(cand.set = Cand.mod, modnames = Modnames)
evidence(aictab(cand.set = Cand.mod, modnames = Modnames))

confset(cand.set = Cand.mod, modnames = Modnames, second.ord = TRUE,
        method = "raw")


modavg(Cand.mod, "VegDiversity", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

modavg(Cand.mod, "NMDS1", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

modavg(Cand.mod, "NMDS2", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

# si la varianza de Moonlight es cerca de 0 en todos los modelos de Arctiinae, tal vez puedo
# remover Moonlight de todos los modelos y justificarlo bien. 
# Tal vez la Luna afecta las Arctiinae mucho menos que a las Geometridae... Revisar esto con estadistica!




# Structural model -------------------------------------------------------------------------

# GEOMETRIDAE structural models --------------------------------------------------------------

gs.null <- glmer(G_fisher ~ 1+(1|Habitat),
                 data = data_all, family=Gamma(link = inverse)); summary(gs.null)

gs1 <- glmer(G_fisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Habitat),
             data = data_all, family=Gamma(link = log)); summary(gs1)  # does not converge with inverse

gs2 <- glmer(G_fisher ~ UnderComplex + CanopyCover +(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs2)

gs3 <- glmer(G_fisher ~ UnderComplex+VerticalComplex+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs3)

gs4 <- glmer(G_fisher ~ CanopyCover + VerticalComplex +(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs4)

gs5 <- glmer(G_fisher ~ UnderComplex+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs5) # does not converge

gs6 <- glmer(G_fisher ~ CanopyCover+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs6) # does not converge

gs7 <- glmer(G_fisher ~ VerticalComplex+(1|Habitat),
             data = data_all, family=Gamma(link = inverse)); summary(gs7)

# gs8 <- glmer(G_fisher ~ (UnderComplex * CanopyCover * VerticalComplex) +(1|Moonlight)+(1|Habitat),
#             data = data_all, family=Gamma(link = log)); summary(gs8) # does not converge

# gs9 <- glmer(G_fisher ~ (UnderComplex * CanopyCover) +(1|Moonlight)+(1|Habitat),
#             data = data_all, family=Gamma(link = log)); summary(gs9) # does not converge

# gs10 <- glmer(G_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
#              data = data_all, family=Gamma(link = log)); summary(gs10)

# gs11 <- glmer(G_fisher ~ (CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
#              data = data_all, family=Gamma(link = log)); summary(gs11) # does not converge



ggplot(data_all,aes(x=Habitat,y=G_fisher)) + geom_jitter() + geom_boxplot(alpha=0.2) 
ggplot(data_all,aes(x=Moonlight,y=G_fisher)) + geom_jitter() + geom_point(alpha=0.2) 


tab_model(gs.null,gs1,gs2,gs3,gs4,gs5,gs6,gs7, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

models <- list(gs.null,gs1,gs2,gs3,gs4,gs5,gs6,gs7)
model.names <- c("gs.null","gs1","gs2","gs3","gs4","gs5","gs6","gs7")
aictab(cand.set = models, modnames = model.names)

# According to this, the null model has the lowest AIC
# There is no competing model. Does this mean that the null model is the best?
# With glmulti, the best model included (1 + UnderComplex + VerticalComplex)


# ARCTIINAE structural models --------------------------------------------------------------
# removed Moonlight as random effect

as.null <- lmer(A_fisher ~ 1+(1|Habitat),
                data = data_all, REML = FALSE); summary(as.null) 

as1 <- lmer(A_fisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Habitat),
            data = data_all, REML = FALSE); summary(as1)    # is Singular
plot(as1)     # para ver los residuales?

as2 <- lmer(A_fisher ~ UnderComplex + CanopyCover+(1|Habitat),
            data = data_all, REML = FALSE); summary(as2)  

as3 <- lmer(A_fisher ~ UnderComplex+VerticalComplex+(1|Habitat),
            data = data_all, REML = FALSE); summary(as3) 

as4 <- lmer(A_fisher ~ CanopyCover + VerticalComplex+(1|Habitat),
            data = data_all, REML = FALSE); summary(as4) 

as5 <- lmer(A_fisher ~ UnderComplex+(1|Habitat),
            data = data_all, REML = FALSE); summary(as5)  

as6 <- lmer(A_fisher ~ CanopyCover+(1|Habitat),
            data = data_all, REML = FALSE); summary(as6)  

as7 <- lmer(A_fisher ~ VerticalComplex+(1|Habitat),
            data = data_all, REML = FALSE); summary(as7)  
plot(as7)      # para ver los residuales?

# as8 <- lmer(A_fisher ~ (UnderComplex*CanopyCover*VerticalComplex) +(1|Moonlight)+(1|Habitat),
#            data = data_all); summary(as8)  # in Singular 

# as9 <- lmer(A_fisher ~ (UnderComplex*CanopyCover)+(1|Moonlight)+(1|Habitat),
#            data = data_all); summary(as9)  

# as10 <- lmer(A_fisher ~ (UnderComplex*VerticalComplex)+(1|Moonlight)+(1|Habitat),
#             data = data_all); summary(as10) 

# as11 <- lmer(A_fisher ~ (CanopyCover*VerticalComplex)+(1|Moonlight)+(1|Habitat),
#             data = data_all); summary(as11) 


tab_model(as.null,as1,as2,as3,as4,as5,as6,as7, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

anova(a.null,a1,a2,a3,a4,a5,a6,a7)

models <- list(as.null,as1,as2,as3,as4,as5,as6,as7)
model.names <- c("as.null","as1","as2","as3","as4","as5","as6","as7")
aictab(cand.set = models, modnames = model.names)

# According to this, the model that includes (1+UnderComplex + CanopyCover + VerticalComplex) 
# has the lowest AIC, but there are several competing models. I need to use model averaging.
# With glmulti, the best model was (1+UnderComplex+VerticalComplex).


# This method below is to compare between the best candidate models
# this is called model averaging
# https://danstich.github.io/stich/classes/BIOL678/06_modelSelection.html

Cand.mod <- list()
Cand.mod[[1]]<-as1
Cand.mod[[2]]<-as4
Cand.mod[[3]]<-as2
Cand.mod[[4]]<-as3


Modnames <- c("UnderComplex + CanopyCover + VerticalComplex",
              "CanopyCover + VerticalComplex", 
              "UnderComplex + CanopyCover",
              "UnderComplex + VerticalComplex") 
aictab(cand.set = Cand.mod, modnames = Modnames)
evidence(aictab(cand.set = Cand.mod, modnames = Modnames))

confset(cand.set = Cand.mod, modnames = Modnames, second.ord = TRUE,
        method = "raw")


modavg(Cand.mod, "UnderComplex", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

modavg(Cand.mod, "CanopyCover", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

modavg(Cand.mod, "VerticalComplex", modnames=Modnames, c.hat = 1, gamdisp = NULL,
       conf.level = 0.95, second.ord = TRUE, nobs = NULL,
       exclude = list(""), warn = TRUE, uncond.se = "revised",
       parm.type = NULL)

# all the variables pass by cero... None of them are significant?


# -------------------------------------------------------------------------

# Geometridae models using log transformed G_fisher. ----------------------
# (lmer instead of glmer with gamma dist)

# Floristic

gft.null <- lmer(logGfisher ~ 1+(1|Habitat),
                 data = data_all, REML = FALSE); summary(gft.null)  

gft1 <- lmer(logGfisher ~ VegDiversity+NMDS1+NMDS2+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft1)  # is Singular

gft2 <- lmer(logGfisher ~ VegDiversity+NMDS1+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft2) 

gft3 <- lmer(logGfisher ~ VegDiversity+NMDS2+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft3) 

gft4 <- lmer(logGfisher ~ NMDS1+NMDS2+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft4)  # is Singular

gft5 <- lmer(logGfisher ~ VegDiversity+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft5)

gft6 <- lmer(logGfisher ~ NMDS1+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft6)  # is Singular

gft7 <- lmer(logGfisher ~ NMDS2+(1|Habitat),
             data = data_all, REML = FALSE); summary(gft7)  


tab_model(gft.null,gft1,gft2,gft3,gft4,gft5,gft6,gft7,gft8,gft9,gft10,gft11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

# According to this, the null model has the lowest AIC
# Using glmulti, the best model was (1 + NMDS1)

# Structural

gst.null <- lmer(logGfisher ~ 1+(1|Habitat),
                 data = data_all, REML = FALSE); summary(gst.null) 

gst1 <- lmer(logGfisher ~ UnderComplex + CanopyCover + VerticalComplex +(1|Habitat),
             data = data_all, REML = FALSE); summary(gst1)  # is Singular

gst2 <- lmer(logGfisher ~ UnderComplex + CanopyCover+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst2)  

gst3 <- lmer(logGfisher ~ UnderComplex+VerticalComplex+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst3)  

gst4 <- lmer(logGfisher ~ CanopyCover + VerticalComplex+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst4) # is Singular

gst5 <- lmer(logGfisher ~ UnderComplex+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst5)  

gst6 <- lmer(logGfisher ~ CanopyCover+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst6)  

gst7 <- lmer(logGfisher ~ VerticalComplex+(1|Habitat),
             data = data_all, REML = FALSE); summary(gst7)  


tab_model(gst.null,gst1,gst2,gst3,gst4,gst5,gst6,gst7,gst8,gst9,gst10,gst11, show.aic = TRUE, show.aicc = TRUE, show.fstat = TRUE)

# According to this, the null model has the lowest AIC
# Using glmulti, the best model was (1 + VerticalComplex)


