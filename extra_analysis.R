### Extra Analysis ###
########################### Black Carbon ########################################
# Using FTMS3 from ICBMOcean
# Needs: 
# ftms3
library(readr)
library(dplyr)

ftms3 <- read_csv("data/ftms3.csv")

ftms3_w_bc <- ftms3 %>% 
  mutate(bc = if_else(condition = ftms3$AI.mod <= 0.66,
                         true = "no",
                         false = "yes")) # Characterizing BC molecules based on aimod


water <- read.csv('data/watersheds_characterized.csv',header=T)

water$duplicate <- as.factor(water$duplicate) 
water$treatment <- as.factor(water$treatment)
water$trip <- as.factor(water$trip)


# Filtering on treatment
burn_names <- filter({water}, treatment=="burn", duplicate=="no") 
burn_names2 <- burn_names$ftms_name

control_names <- filter({water}, treatment=="control", duplicate=="no") 
control_names2 <- control_names$ftms_name


# Weighted Averages of BC
rel_absb <- apply(ftms3_w_bc[,c(burn_names2)],2,function(x){x/sum(x,na.rm=T)})
bc_b <- apply(rel_absb,2,function(x){sum(x*ftms3_w_bc$bc,na.rm=T)})



rel_absc <- apply(ftms3_w_bc[,c(control_names2)],2,function(x){x/sum(x,na.rm=T)})
bc_c <- apply(rel_absc,2,function(x){sum(x*ftms3_w_bc$bc,na.rm=T)})





#################################################### Other ###################################################################
# Count rows with all NA's except the first column
## for reporting peaks in results section
count <- sum(rowSums(is.na(peakObj[["e_data"]][, -1])) == (ncol(peakObj[["e_data"]]) - 1))
print(count)



# of X ± Y per sample in burned and Z ± ZZ in unburned sites
library(Rmisc)
df <-  subset(ftms3, select = c(control_names2)) #or burn_names2
count1 <- colSums(!is.na(df)) %>% 
  as.data.frame()
print(count1)
summarySE(count1, measurevar = ".")

df2 <-  subset(ftms3, select = c(burn_names2)) #or burn_names2
count2 <- colSums(!is.na(df2)) %>% 
  as.data.frame()
print(count2)
summarySE(count2, measurevar = ".")

write.csv(count1, file = "df_control.csv")
write.csv(count2, file = "df_burn.csv")



# checking number of mfs
hist(mfs$count)

t.test(count~type, paired=FALSE , var.equal=TRUE, data=mfs) # checking which r thinks is the first group
t.test(count~type, alternative="greater", paired=FALSE , var.equal=TRUE, data=mfs)



# counting unique MF's
view(treat_unique[["e_data"]])
table(treat_unique[["e_data"]]$uniqueness_gtest)