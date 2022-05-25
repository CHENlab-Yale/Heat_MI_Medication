##########################################################################################
#### Triggering of myocardial infarction by heat exposure modified by medication intake
#### R codes for main model analysis
#### Kai Chen, May 2022
#### Helmholtz Zentrum Muenchen; Yale School of Public Health
##########################################################################################
####Load R packages for this analysis
library(dlnm); library(splines); library(gnm)

##Function to calculate the p-value for difference between 2 RRs with 95%CI 
pvalue_dif_RRs <- function(RR1, RR1.low, RR1.high, RR2, RR2.low, RR2.high){
  ##First calculate the est and se for each period
  #RR1
  Est1 <- log(RR1)
  SE1 <- (log(RR1.high) - log(RR1.low))/2/1.96
  #RR2
  Est2 <- log(RR2)
  SE2 <- (log(RR2.high) - log(RR2.low))/2/1.96
  ##Then calculate the difference between log relative risks
  dif_est <- Est1 - Est2
  dif_se <- sqrt(SE1^2 + SE2^2)
  z_score <- dif_est/dif_se
  pvalue <- 2*pnorm(-abs(z_score))
  return(pvalue)
  ##End
}

#######################################################################
###1 . Read the Augsburg MI dataset for warmseason during 2001-2014 ###
#######################################################################
##Given the MI registry data is not public available, we don't share it here.
# data

#######################
## 2.Main analysis ####
#######################
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun = "ns"##natural cubic regression
# SPECIFICATION OF THE LAG FUNCTION
# DEFINE THE CROSSBASIS using ns with 4df for temperature and p-spline for lag structure
df.temp <- 3

###Constract the crossbasis function
cb <- onebasis(data$Tmean,fun=varfun, df=df.temp)
##use nonlinear term for RH
RH.ob <- onebasis(data$RH, fun = "ns", df=3)

##Obtain MT.NonfatalMI from the overall data (without missing MI)
MMT.NonfatalMI <- 7.5 ## degree C

######################################
##Effect modification by Medication ##
######################################
#Note that based on Table1 summary, only drug intakse with >10% in Nonfatal MI will be included in this analysis
##Drug list include Antiplatelet drugs (ANTPLB); ACE_inhibitors (ACEB); Beta-blocker (BBLOKB);
##Calcium_channel_blocker (CABLOB); 
##Diuretics(DIURB); Statins(STATIB)
Subgroup.names <- c("ANTPLB","NoANTPLB","ACEB","NoACEB", "BBLOKB", "NoBBLOKB",  "CABLOB","NoCABLOB",
                    "DIURB","NoDIURB","STATIB","NoSTATIB", 
                    "Anticoagulants", "NoAnticoagulants","Antihypertensives", "NoAntihypertensives",
                    "Lipidlowering", "NoLipidlowering", "Antidiabetics", "NoAntidiabetics")
##corresponding Drug group
Drug.groups <- c("Antiplatelet",  "ACE inhibitors", "Beta blocker",  "Calcium channel blocker",
                 "Diuretics", "Statin", "Anticoagulants","Antihypertensives","Lipidlowering","Antidiabetics")


##dataframe to store the results
TempRR <- matrix(rep(NA,(10*2)*8),ncol=8) # 10 for 10 droug.group; 2 for users and non-users of medication
colnames(TempRR) <- c("Drug.group","Medication","RR","RRlow","RRhigh","p.dif")
TempRR <- as.data.frame(TempRR)

######Loof for each subgroup
for (i in 1:length(Drug.groups)){
  
  #########################
  ##Medications Subgroup ##
  #########################
  ##Medication Subgroup 1 (Yes)
  model1.Subgroup <- gnm(as.formula(paste0(Subgroup.names[2*i-1],"~ cb + RH")),data=data, offset=log_pop, family=quasipoisson, eliminate=factor(stratum))
  MMT1.Subgroup <- MMT.NonfatalMI ## 7.5°C
  
  ##prediction
  pred1.Subgroup <- crosspred(cb,model1.Subgroup, at=c(min(data$Tmean, na.rm=T):max(data$Tmean, na.rm=T),
                                                       Tmean.1st,Tmean.2.5th, Tmean.5th, Tmean.10th,Tmean.90th, Tmean.95th,
                                                       Tmean.97.5th, Tmean.99th),cen= MMT1.Subgroup, bylag = 1, cumul = T)
  
  ##Medication Subgroup 2 (No)
  model2.Subgroup <- gnm(as.formula(paste0(Subgroup.names[2*i],"~ cb + RH")),data=data, offset=log_pop, family=quasipoisson, eliminate=factor(stratum))
  MMT2.Subgroup <- MMT.NonfatalMI
  ##prediction
  pred2.Subgroup <- crosspred(cb,model2.Subgroup, at=c(min(data$Tmean, na.rm=T):max(data$Tmean, na.rm=T),
                                                       Tmean.1st,Tmean.2.5th, Tmean.5th, Tmean.10th,Tmean.90th, Tmean.95th,
                                                       Tmean.97.5th, Tmean.99th),cen= MMT2.Subgroup, bylag = 1, cumul = T)
  ################################
  ###Calculate the heat effects ##
  ################################
  ####Drug.group
  TempRR$Drug.group[(2*i-1):(2*i)] <- Drug.groups[i]
  
  ####Medication
  TempRR$Medication[2*i-1] <- Subgroup.names[2*i-1]
  #95th vs MMT
  TempRR$RR[2*i-1] <- pred1.Subgroup$allRRfit[as.character(Tmean.95th)]
  TempRR$RRlow[2*i-1] <- pred1.Subgroup$allRRlow[as.character(Tmean.95th)]
  TempRR$RRhigh[2*i-1]<- pred1.Subgroup$allRRhigh[as.character(Tmean.95th)]
  
  ####Medication
  TempRR$Medication[2*i] <- Subgroup.names[2*i]
  #95th vs MMT
  TempRR$RR[2*i] <- pred2.Subgroup$allRRfit[as.character(Tmean.95th)]
  TempRR$RRlow[2*i] <- pred2.Subgroup$allRRlow[as.character(Tmean.95th)]
  TempRR$RRhigh[2*i]<- pred2.Subgroup$allRRhigh[as.character(Tmean.95th)]

  ##Heat (99th vs MMT)
  ###Check the difference between users and nonusers of a medication 
  #95th vs MMT
  TempRR$p.dif[2*i-1] <- pvalue_dif_RRs(RR1=TempRR$RR[2*i-1], RR1.low =TempRR$RRlow[2*i-1], RR1.high =TempRR$RRhigh[2*i-1], 
                                             RR2 =TempRR$RR[2*i], RR2.low =TempRR$RRlow[2*i], RR2.high =TempRR$RRhigh[2*i])
  
  ###End analysis for the subgroup
  cat(paste("End calculation for drug.group:", Drug.groups[i],"!\n"))
}

##Get the results for heat effects stratified by medication intakes
TempRR
