################################################# LMER models ###############################################
library(lme4)
library(glmmTMB)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(ggplot2)
library(scales)


# read in data
water1 <- read.csv('data/watersheds_characterized.csv',header=T)
water2 <- read.csv('data/watersheds_characterized1.csv',header=T)
water <- merge(x = water1, y= water2[, c("sample_id", "bc")], by = "sample_id", all.x = TRUE)

# drop duplicate values
water_nodup <- water[which(water$duplicate == 'no'),]


# first check what the data and histograms
water_nodup$logcond_us_cm = log(water_nodup$cond_us_cm) 
water_nodup$log_suva = log(water_nodup$SUVA) 
water_nodup$log_e2e3 = log(water_nodup$E2_E3) 

hist(water_nodup$TOC.mgL)
hist(water_nodup$logcond_us_cm)
hist((water_nodup$Comp.7)^(1/4))


# test normality 
shapiro.test(water_nodup$hix)

# relevel the factors
water_nodup$treatment <- relevel(as.factor(water_nodup$treatment), ref='control')
water_nodup$trip <- relevel(as.factor(water_nodup$trip), ref='june')


# fit linear model to test how DOC concentrations vary
# test interaction that estimates different trip value for each treatment
m1 <- lmer(TOC.mgL ~ treatment * trip + (1|watershed), data=water_nodup) # add ,REML = FALSE for aic weights
m2 <- lmer(TN.mgL ~ treatment * trip + (1|watershed), data=water_nodup)
m3 <- lmer(ph ~ treatment * trip  + (1|watershed), data=water_nodup)
m4 <- lmer(orp_mv ~ treatment * trip  + (1|watershed), data=water_nodup)
m5 <- lmer(logcond_us_cm ~ treatment * trip + (1|watershed), data=water_nodup)
m6 <- lmer(temp ~ treatment * trip  + (1|watershed), data=water_nodup)
m7 <- lmer(do2_mg_l ~ treatment * trip  + (1|watershed), data=water_nodup)

m8 <- lmer(h.c.wtavg ~ treatment * trip + (1|watershed), data=water_nodup) # H:C
m9 <- lmer(o.c.wtavg ~ treatment * trip + (1|watershed), data=water_nodup) # O:C
m10 <- lmer(GFE ~ treatment * trip + (1|watershed), data=water_nodup) # GFE
m11 <- lmer(AI_Mod ~ treatment * trip + (1|watershed), data=water_nodup) # AI mod
m12 <- lmer(bc ~ treatment * trip + (1|watershed), data=water_nodup) 
m13 <- lmer(NOSC ~ treatment * trip + (1|watershed), data=water_nodup) # NOSC

water_nodup$biot <- water_nodup$transformation/water_nodup$X
m14a <- lmer(biot ~ treatment * trip + (1|watershed), data=water_nodup) 

m15 <- lmer(bix ~ treatment * trip  + (1|watershed), data=water_nodup) 
m16 <- lmer(hix ~ treatment * trip  + (1|watershed), data=water_nodup) 
m17 <- lmer(log_suva ~ treatment * trip + (1|watershed), data=water_nodup)
m18 <- lmer(SR ~ treatment * trip + (1|watershed), data=water_nodup)
m19 <- lmer(fi ~ treatment * trip + (1|watershed), data=water_nodup)
m20 <- lmer(log_e2e3~ treatment * trip + (1|watershed), data=water_nodup)

# Component Analysis Models
m21 <- lmer(Comp.1 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m22 <- lmer(Comp.2 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m23 <- lmer(Comp.3 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m24 <- lmer(Comp.4 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m25 <- lmer(Comp.5 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m26 <- lmer(Comp.6 ~ haifls_mod * trip + (1|watershed), data=water_nodup)
m27 <- lmer(I(Comp.7^(1/4)) ~ haifls_mod * trip + (1|watershed), data=water_nodup)


# Model Summaries
summary(m1)
aic1 <- AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27) # put in AIC table for results write up
pairs(emmeans(m1,'treatment'))
pairs(emmeans(m14,~treatment|trip, type = "response")) #for glmer model
pairs(emmeans(m1,~treatment|trip)) # use for the differences when writing up results
emmeans(m1, ~treatment|trip) # use for means in results write up


# great result - check that the model meets assumptions
# fitted vs residual values should look like the night sky = no directional trend and even spread of points across x-axis
plot(m14b)
hist(resid(m21))
# check what the R2 value of the model - how well does the model fit the data
# we have R2m = marginal R2 (fixed effects only) and R2c = conditional R2 (fixed + random effects)
# the R2m = 20% of variation in TOC explained by trip and treatment, with 69% (R2c) explained by "watershed" effect 
# overall we explain 89% of the variation in TOC with our model
r.squaredGLMM(m4)


# look at alternatives to treatment
## create new variable for burn_mean
water_nodup$burn_mean_mod <- water_nodup$burn_mean
water_nodup$ifls_mod <- water_nodup$hw_burn_mean_ifls
hist(water_nodup$ifls_mod)
water_nodup$haifls_mod <- water_nodup$hw_burn_mean_haifls
hist(water_nodup$haifls_mod)
# set everything unburned at zero
water_nodup$burn_mean_mod[which(water_nodup$burn_mean_mod < 200)] <- 0
water_nodup$ifls_mod[which(water_nodup$ifls_mod < 200)] <- 0
water_nodup$haifls_mod[which(water_nodup$haifls_mod < 176)] <- 0

mca <- lmer(TOC.mgL ~ burn_mean * trip + (1|watershed), data=water_nodup)
mcb <- lmer(TOC.mgL ~ burn_mean_mod * trip + (1|watershed), data=water_nodup)
summary(water_nodup$burn_mean)





################################################### Box plots ######################################################
# Response variable individual plots
toc <- ggplot(water_nodup) +
 aes(x = treatment, y = TOC.mgL, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "DOC (mg/L)", tag = "a.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 65, label = "* (Aug only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

tn <- ggplot(water_nodup) +
 aes(x = treatment, y = TN.mgL, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "DN (mg/L)", tag = "f.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")  +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

ph <- ggplot(water_nodup) +
 aes(x = treatment, y = ph, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "pH", tag = "b.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 8, label = "* (Aug only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

temp <- ggplot(water_nodup) +
 aes(x = treatment, y = temp, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "Temperature (°C)", tag = "d.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 20, label = "* (All Months)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

orp <- ggplot(water_nodup) +
 aes(x = treatment, y = orp_mv, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "ORP (mV)", tag = "e.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 425, label = "* (All Months)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

do2 <- ggplot(water_nodup) +
 aes(x = treatment, y = do2_mg_l, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "DO2 (mg/L)") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

logcond <- ggplot(water_nodup) +
 aes(x = treatment, y = logcond_us_cm, fill = treatment) +
 geom_boxplot() +
 scale_fill_brewer(palette = "Pastel2", 
 direction = 1) +
 labs(y = "Conductivity (µs/cm)", tag = "c.") +
 theme_classic() +
 theme(plot.title = element_text(hjust = 0.5),
       axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 6, label = "* (Aug only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn")) +
 scale_y_continuous(trans = scales::log_trans())



########################################### OPtical #######################################################
suva <- ggplot(water_nodup) +
  aes(x = treatment, y = log_suva, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "SUVA", tag = "e.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))+
  scale_y_continuous(trans = scales::log_trans())

sr <- ggplot(water_nodup) +
  aes(x = treatment, y = SR, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "SR", tag = "f.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

fi <- ggplot(water_nodup) +
  aes(x = treatment, y = fi, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "FI", tag = "b.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 1.4, label = "* (July only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

hix <- ggplot(water_nodup) +
  aes(x = treatment, y = hix, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "HIX", tag = "d.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

bix <- ggplot(water_nodup) +
  aes(x = treatment, y = bix, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "BIX", tag = "a.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 0.525, label = "* (July only)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

e2e3 <- ggplot(water_nodup) +
  aes(x = treatment, y = log_e2e3, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "E2:E3", tag = "c.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))+
  annotate("text", x = "burn", y = 2, label = "* (All Months)") +
  scale_y_continuous(trans = scales::log_trans())



################################################### Components ###########################################################
comp1 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.1, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.1", tag = "c.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp2 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.2, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.2", tag = "a.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 0.3, label = "* (Aug only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp3 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.3, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.3", tag = "d.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp4 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.4, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.4", tag = "b.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 0.5, label = "* (Aug only)") +
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp5 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.5, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.5", tag = "e.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp6 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.6, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.6", tag = "f.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

comp7 <- ggplot(water_nodup) +
  aes(x = treatment, y = Comp.7, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Comp.7", tag = "g.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))



############################################## FTMS Graphing for models ############################################################################
# Response variable individual plots
h.c <- ggplot(water_nodup) +
  aes(x = treatment, y = h.c.wtavg, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "H:C", tag = "a.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12),)+
  annotate("text", x = "burn", y = 1.15, label = "* (All Months)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

o.c <- ggplot(water_nodup) +
  aes(x = treatment, y = o.c.wtavg, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "O:C") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14), legend.position = "none")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

gfe <- ggplot(water_nodup) +
  aes(x = treatment, y = GFE, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "GFE (J)", tag = "f.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 66, label = "* (July only)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

ai_mod <- ggplot(water_nodup) +
  aes(x = treatment, y = AI_Mod, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "AI Mod", tag = "b.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 0.375, label = "* (All Months)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

nosc <- ggplot(water_nodup) +
  aes(x = treatment, y = NOSC, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "NOSC", tag = "e.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12),  legend.position = "none")+
  annotate("text", x = "burn", y = 0.1, label = "* (July only)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

biot <- ggplot(water_nodup) +
  aes(x = treatment, y = biot, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Biotransformations", tag = "d.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none")+
  annotate("text", x = "burn", y = 0.08, label = "* (All Months)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

bc <- ggplot(water_nodup) +
  aes(x = treatment, y = bc*100, fill = treatment) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel2", 
                    direction = 1) +
  labs(y = "Abundance of black carbon (%)", tag = "c.") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), legend.position = "none") +
  annotate("text", x = "burn", y = 3.25, label = "* (All Months)")+
  scale_x_discrete(labels = c("control" = "Unburned","burn" = "Burn"))

#################################################### combining into a panel with grid.arrange #################################################################
library(gridExtra)

#legend <- gtable::gtable_filter(ggplotGrob(for_legend), "guide-box") #if needed
box_chem <- grid.arrange(toc, ph, logcond, temp, orp, tn, nrow=2) 
box_dom <- grid.arrange(h.c, ai_mod, bc, biot, nosc, gfe, nrow=2) 
box_opt <- grid.arrange(bix, fi, e2e3, hix, suva, sr, nrow=2)
box_com <- grid.arrange(comp2, comp4, comp1, comp3, comp5, comp6, comp7, nrow=3)






################################### Table ####################################################
library(sjPlot)

table1 <- tab_model(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14b,m15,m16,m17,m18,m19,m20, 
          show.df = TRUE, show.se = TRUE, show.ci = FALSE, show.est = FALSE, show.icc = FALSE, show.ngroups = FALSE, show.re.var = TRUE, show.stat = TRUE,
          pred.labels = c("Intercept",
                              "Treatment: Burn",
                              "Trip: August",
                              "Trip: July",
                              "Treatment: Burn x Trip: August",
                              "Treatment: Burn x Trip: July"),
          dv.labels = c("DOC","DN","pH","ORP","Conductivity","Temp","DO2","H:C","O:C","GFE","AI Mod","Black Carbon-like","NOSC","Putative Transformations","BIX","HIX","SUVA","SR","FI","E2 E3"))
