### Initial FTMS R Processing  ###
################################# Non ftmsRanalysis package pre processing ##########################################################
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Reading in data 
ftms3 <- read_csv("data/ftms3.csv")
water <- read.csv('data/watersheds_characterized.csv',header=T)

water$duplicate <- as.factor(water$duplicate) 
water$treatment <- as.factor(water$treatment)
water$trip <- as.factor(water$trip)

# Filtering on treatment 
burn_names <- filter({water}, treatment=="burn", duplicate=="no") 
burn_names2 <- burn_names$ftms_name

control_names <- filter({water}, treatment=="control", duplicate=="no") 
control_names2 <- control_names$ftms_name

# Weighted Averages H:C and O:C
rel_absb <- apply(ftms3[,c(burn_names2)],2,function(x){x/sum(x,na.rm=T)})
sample_h.c.wtavg_b <- apply(rel_absb,2,function(x){sum(x*ftms3$H.C,na.rm=T)})
sample_o.c.wtavg_b <- apply(rel_absb,2,function(x){sum(x*ftms3$O.C,na.rm=T)})


rel_absc <- apply(ftms3[,c(control_names2)],2,function(x){x/sum(x,na.rm=T)})
sample_h.c.wtavg_c <- apply(rel_absc,2,function(x){sum(x*ftms3$H.C,na.rm=T)})
sample_o.c.wtavg_c <- apply(rel_absc,2,function(x){sum(x*ftms3$O.C,na.rm=T)})





######################################### non ftmsRanalysis VanKrevlan plot ############################################################
# all samples plot 
with(ftms3,plot(O.C,H.C))

# burn VK plot 
VK_sub1 <- apply(ftms3[,c(burn_names2)],1,function(x){as.numeric(sum(x,na.rm=T)>0)}) #subseting based on names
VK_OC_sub1 <- ftms3$O.C*VK_sub1
VK_HC_sub1 <- ftms3$H.C*VK_sub1
VK_OC_sub1 <- VK_OC_sub1[which(VK_HC_sub1!=0)]
VK_HC_sub1 <- VK_HC_sub1[which(VK_HC_sub1!=0)]

plot(VK_OC_sub1,VK_HC_sub1,pch=19,col=rgb(0,0,0,.1))

# control VK plot
VK_sub_c <- apply(ftms3[,c(control_names2)],1,function(x){as.numeric(sum(x,na.rm=T)>0)}) #subseting based on names
VK_OC_sub_c <- ftms3$O.C*VK_sub_c
VK_HC_sub_c <- ftms3$H.C*VK_sub_c
VK_OC_sub_c <- VK_OC_sub_c[which(VK_HC_sub_c!=0)]
VK_HC_sub_c <- VK_HC_sub_c[which(VK_HC_sub_c!=0)]

par(mfrow=c(1,2))
plot(VK_OC_sub1,VK_HC_sub1,pch=19,col=rgb(0,0,0,.1))
plot(VK_OC_sub_c,VK_HC_sub_c,pch=19,col=rgb(0,0,0,.1)) 

# Notch histogram for H:C and O:C
par(mfrow=c(1,2))
boxplot(c(sample_h.c.wtavg_b,sample_h.c.wtavg_c)~c(rep('b',length(sample_h.c.wtavg_b)),rep('c',length(sample_h.c.wtavg_c))),notch=T)
boxplot(c(sample_o.c.wtavg_b,sample_o.c.wtavg_c)~c(rep('b',length(sample_o.c.wtavg_b)),rep('c',length(sample_o.c.wtavg_c))),notch=T)






########################################## Guided package for processing ftms data ####################################################
library(ftmsRanalysis)

# Import data
hotw_edata <- read.csv("data/hotw_edata.csv")
hotw_fdata <- read.csv("data/hotw_fdata.csv")
hotw_emeta <- read.csv("data/hotw_emeta.csv")


# Combine into one object
peakObj <- as.peakData(hotw_edata, hotw_fdata, hotw_emeta,             
                       edata_cname="mz", fdata_cname="ftms_name", 
                       mass_cname="mz", c_cname="C", h_cname="H", 
                       o_cname="O", n_cname="N", s_cname="S", 
                       p_cname="P", isotopic_cname = "C13", 
                       isotopic_notation = "3") # needs to be 3 to keep in all masses from ftms3 file
plot(peakObj) 


# Pre processing
peakObj <- edata_transform(peakObj, data_scale="log2") 

# For presence/absence transformation:
edata_transform(peakObj, data_scale="pres")

# Calculating meta data 
peakObj <- compound_calcs(peakObj)
head(peakObj$e_meta)

peakObj <- assign_elemental_composition(peakObj)
table(peakObj$e_meta[,getElCompColName(peakObj)])

peakObj <- assign_class(peakObj, boundary_set = "bs1") # Will need this for stacked compound class graphing 
table(peakObj$e_meta[, getBS1ColName(peakObj)])



# VanKrevlan  of one sample 
one_sample <- subset(peakObj, samples="Sample14_erin_mass_listhotw_0010_5ppm_8m_000002.csv") # use burn_names and control_names to subset burned and unburned
summary(one_sample)
vanKrevelenPlot(one_sample, title="rb01 example")

# VanKrevlan of each treatment 
control_plot <- subset(peakObj, samples=c(control_names2)) # use burn_names and control_names2 to subset burned and unburned
vanKrevelenPlot(control_plot, title="control")

burn_plot <- subset(peakObj, samples=c(burn_names2)) # use burn_names and control_names2 to subset burned and unburned
vanKrevelenPlot(burn_plot, title="burn")



##################################################### NOSC Graph #####################################################################
#devtools::install_github("delta-rho/datadr") #if needed
library(datadr)

# NOSC Distribution/ histogram 
peakObj <- group_designation(peakObj, main_effects=c("treatment"))
getGroupDF(peakObj)

group_summary <- summarizeGroups(peakObj, summary_functions = 
                                   c("n_present", "prop_present"))
head(group_summary$e_data)

densityPlot(peakObj, samples=FALSE, groups=c("burn","control"), variable="NOSC", #my nosc plot to share with Andrew
            title="Comparison of NOSC Between Treatments") 

# Unique peaks + G-test stats test
byGroup <- divideByGroupComparisons(peakObj, 
                                    comparisons = "all")[[1]]$value

treat_unique <- summarizeGroupComparisons(byGroup, 
                                         summary_functions="uniqueness_gtest", 
                                         summary_function_params=list(
                                           uniqueness_gtest=list(pres_fn="nsamps", 
                                                                 pres_thresh=2, pvalue_thresh=0.05)))

head(treat_unique$e_data)

p <- vanKrevelenPlot(treat_unique, colorCName = "uniqueness_gtest", showVKBounds = FALSE) # showVKBounds is for the compound class boxes

# in the plot treat_unique, set visible = FALSE for the third trace
# which happens to be the points for 'observed in both'
plotly::style(p, visible = FALSE, traces = c(3))





####################################################### Unique MF Plots by month ###############################################################
# Making lists for filtering later
june_names <- filter({water}, duplicate=="no", trip=="june")
june_names2 <- june_names$ftms_name
july_names <- filter({water}, duplicate=="no", trip=="july")
july_names2 <- july_names$ftms_name
august_names <- filter({water}, duplicate=="no", trip=="august")
august_names2 <- august_names$ftms_name

june <- subset(peakObj, samples=c(june_names2))
july <- subset(peakObj, samples=c(july_names2))
august <- subset(peakObj, samples=c(august_names2))


# June
byGroup_june <- divideByGroupComparisons(june, 
                                    comparisons = "all")[[1]]$value

june_unique <- summarizeGroupComparisons(byGroup_june, 
                                          summary_functions="uniqueness_gtest", 
                                          summary_function_params=list(
                                            uniqueness_gtest=list(pres_fn="nsamps", 
                                                                  pres_thresh=2, pvalue_thresh=0.05)))

june_unique$e_data$uniqueness_gtest <- relevel(as.factor(june_unique$e_data$uniqueness_gtest), ref='Unique to control')
cpal2 <- scales::col_factor("Set2", levels = levels(june_unique$e_data$uniqueness_gtest))
june_vk <- vanKrevelenPlot(june_unique, colorPal = cpal2, colorCName = "uniqueness_gtest", showVKBounds = FALSE)
june_fin <- plotly::style(june_vk, visible = FALSE, traces = c(3))  %>% 
  plotly::style(showlegend = F)

# July
byGroup_july <- divideByGroupComparisons(july, 
                                         comparisons = "all")[[1]]$value

july_unique <- summarizeGroupComparisons(byGroup_july, 
                                         summary_functions="uniqueness_gtest", 
                                         summary_function_params=list(
                                           uniqueness_gtest=list(pres_fn="nsamps", 
                                                                 pres_thresh=2, pvalue_thresh=0.05)))

july_unique$e_data$uniqueness_gtest <- relevel(as.factor(july_unique$e_data$uniqueness_gtest), ref='Unique to control')
cpal3 <- scales::col_factor("Set2", levels = levels(july_unique$e_data$uniqueness_gtest))
july_vk <- vanKrevelenPlot(july_unique, colorPal = cpal3,colorCName = "uniqueness_gtest", showVKBounds = FALSE)
july_fin <- plotly::style(july_vk, visible = FALSE, traces = c(3)) %>% 
  plotly::style(showlegend = F)

# August
byGroup_august <- divideByGroupComparisons(august, 
                                         comparisons = "all")[[1]]$value

august_unique <- summarizeGroupComparisons(byGroup_august, 
                                         summary_functions="uniqueness_gtest", 
                                         summary_function_params=list(
                                           uniqueness_gtest=list(pres_fn="nsamps", 
                                                                 pres_thresh=2, pvalue_thresh=0.05)))

august_unique$e_data$uniqueness_gtest <- relevel(as.factor(august_unique$e_data$uniqueness_gtest), ref='Unique to control')
cpal <- scales::col_factor("Set2", levels = levels(august_unique$e_data$uniqueness_gtest))
august_vk <- vanKrevelenPlot(august_unique, colorPal = cpal, colorCName = "uniqueness_gtest", showVKBounds = FALSE)
aug_fin <- plotly::style(august_vk, visible = FALSE, traces = c(3))

library(plotly)
subplot(june_fin, july_fin, aug_fin, nrows = 1)


########### GG plot Unique MFs ################
# use if you want to add more different things to the graph like black carbon moleclues a different shape

## make a data set in long form
monthly_vk <- peakObj$e_meta %>% left_join(june_unique$e_data %>% rename("june" = "uniqueness_gtest"))%>% 
  left_join(july_unique$e_data %>% rename("july" = "uniqueness_gtest"))%>% 
  left_join(august_unique$e_data %>% rename("august" = "uniqueness_gtest")) %>% 
  pivot_longer(cols = c("june","july","august"), names_to = "month", values_to = "unique_group") %>% 
  filter(!is.na(unique_group)) %>% 
  mutate(bc = ifelse(AI_Mod <= 0.66, "no","yes")) %>% 
  filter(!unique_group == "Observed in Both")

monthly_vk$month <- factor(monthly_vk$month, levels = c("june", "july", "august"))
month_names <- as_labeller(c("june" = "a. June", "july" = "b. July", "august" = "c. August")) # maybe take off a.b.c. based on comments

ggplot(monthly_vk, aes(x = O.C, y = H.C, shape = bc, color = bc)) +
  geom_point(size = 3, alpha = 0.4, aes(fill = unique_group)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_color_manual(values = c("transparent", "black")) +
  scale_fill_manual(values = c("#66c2a5","#fc8d62")) +
  facet_grid(.~month, labeller = month_names) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21)), 
         shape = "none", 
         color = "none") +
  labs(fill = "Uniqueness G Test")+
  theme(legend.position = "none", axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),strip.text.x = element_text(size = 12),
         panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.spacing = unit(1, "lines"))
# 1000X375 pixels to make all three graphs square

uni_count <- filter({monthly_vk}, month=="august", unique_group=="Unique to burn") # for counting in each month by type of unique molecule


# t tests to compare unique mf's in each month
june2 <- filter(monthly_vk, month=="june")
july2 <- filter(monthly_vk, month=="july")
aug2 <- filter(monthly_vk, month=="august")
hist(june2$OtoC_ratio)

#June
t.test(HtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=june2)
t.test(OtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=june2)
#July
t.test(HtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=july2)
t.test(OtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=july2)
#August
t.test(HtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=aug2)
t.test(OtoC_ratio~unique_group, paired=FALSE, var.equal=TRUE, data=aug2)





########################################### Getting other useful data from peakObj #############################################################
# GFE                                                                        
sample_gfe_b <- apply(rel_absb,2,function(x){sum(x*peakObj[["e_meta"]]$GFE,na.rm=T)}) 
sample_gfe_c <- apply(rel_absc,2,function(x){sum(x*peakObj[["e_meta"]]$GFE,na.rm=T)}) 

# NOSC
nosc_b <- apply(rel_absb,2,function(x){sum(x*peakObj[["e_meta"]]$NOSC,na.rm=T)}) 
nosc_c <- apply(rel_absc,2,function(x){sum(x*peakObj[["e_meta"]]$NOSC,na.rm=T)}) 

# AI Mod
ai_mod_b <- apply(rel_absb,2,function(x){sum(x*peakObj[["e_meta"]]$AI_Mod,na.rm=T)}) 
ai_mod_c <- apply(rel_absc,2,function(x){sum(x*peakObj[["e_meta"]]$AI_Mod,na.rm=T)})

# DBE_1
dbe1_b <- apply(rel_absb,2,function(x){sum(x*peakObj[["e_meta"]]$DBE_1,na.rm=T)}) 
dbe1_c <- apply(rel_absc,2,function(x){sum(x*peakObj[["e_meta"]]$DBE_1,na.rm=T)})

write.csv(bs1_c, file = "bs1c.csv") 





############################################ Tidy Verse to figure out BS1 Classes #######################################################
# Needs:
# peakObj from ftmsRanalysis
# subset sample names ie burn_names2 and control_names2
library(tidyverse)

# Combining subclasses into the most prominent class
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Carbohydrate;Amino Sugar'] <- 'Carbohydrate'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Lignin;Amino Sugar'] <- 'Lignin'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Lignin;Tannin'] <- 'Lignin'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Lipid;Protein'] <- 'Lipid'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Protein;Amino Sugar'] <- 'Protein'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Protein;Lignin'] <- 'Protein'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Protein;Lignin;Amino Sugar'] <- 'Protein'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Tannin;Cond Hydrocarbon'] <- 'Tannin'
peakObj$e_meta$bs1_class[peakObj$e_meta$bs1_class == 'Unsat Hydrocarbon;Cond Hydrocarbon'] <- 'Unsat Hydrocarbon'


relative_abundance <- apply(peakObj[["e_data"]][,c(burn_names2, control_names2)],2,function(x){x/sum(x,na.rm=T)}) %>%  # relative abundance
  as.data.frame() %>%  # turn into df to use tidyverse functions
  rownames_to_column("id") %>% # giving id column to match to e_meta file
  pivot_longer(-"id", names_to = "sample", values_to = "value") %>%  # pivoting to long style to have all relative abundances in one column
  left_join(peakObj[["e_meta"]] %>% rownames_to_column("id")) %>%  # joining to e meta
  left_join(peakObj[["f_data"]], by = c("sample" = "ftms_name")) %>%  # joining to f data by sample name
  filter(!is.na(value))# removing na's to clean up


stack <- relative_abundance %>% 
  group_by(trip, treatment) %>%
  # find what the total amount of counts is for the group
  mutate(overall_abund = sum(value)) %>%
           #Then we are grouping by the type of compound, and getting the relative amount of that compound compared to the totals
           group_by(treatment, trip, bs1_class) %>%
           summarize(value = sum(value)/unique(overall_abund)) %>% 
  filter(!bs1_class == "Other")

stack$treatment <- factor(stack$treatment, levels = c("control","burn"))
stack$trip <- factor(stack$trip, levels = c("june", "july", "august"))

ggplot(stack, aes(x = interaction(treatment, trip), value*100, fill = bs1_class)) + # to get june.c july.c aug.c... swap treatment and trip
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette="YlOrRd") +
  theme_classic() +
  labs(x = " June                                       July                                     August", 
       y = "% relative abundance", 
       fill = "",) +
  scale_x_discrete(labels = c("control.june" = "Unburned",
                              "control.july" = "Unburned",
                              "control.august" = "Unburned",
                              "burn.june" = "Burn",
                              "burn.july" = "Burn",
                              "burn.august" = "Burn"))+
  theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12))
