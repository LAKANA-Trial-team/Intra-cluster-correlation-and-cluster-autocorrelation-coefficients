


############# Longitudinal CAC
#### Exchangeable, Nested-exchageable, Exponential decay

library(lme4)
library(glmmTMB)
library(tidyverse)
library(boot)


get_AR1 = function(data = growth_data,outcome = "WAZ", clusterID = "anon_villid",visit = "visit", covariates = c("intv")) {
  
  data[,c(visit, clusterID, covariates)] <- 
    lapply(data[,c(visit, clusterID, covariates)], factor)
  
  fixed_effects = paste(c(visit,covariates),collapse = " + ")
  
  formula_ar1 = as.formula(paste(outcome," ~ ",fixed_effects, "+ ar1(",visit,"+ 0 | ",clusterID,")"))
  
  fit_AR <- glmmTMB(formula = formula_ar1, REML=T, 
                    data = data, family = gaussian)
  V3 <- VarCorr(fit_AR)
  
  ar1 <- round(attr(V3$cond$anon_villid,"correlation")[1,2],4)
  
  return(ar1)
  
}


### Design effect
# Parallel
de_parallel <- function(m, icc) {
  1 + (m - 1) * icc
}


corrStuct = function(data,outcome,clusterID,visit,covariates) {
  
  data[,c(visit, clusterID, covariates)] <- 
    lapply(data[,c(visit, clusterID, covariates)], factor)
  
  data = data[which(!is.na(data[,outcome])),]
  
  fixed_effects = paste(c(visit,covariates),collapse = " + ")
  
  formula_exch = as.formula(paste(outcome," ~ ",fixed_effects, "+ (1 |",clusterID,")"))
  formula_nest_exch = as.formula(paste(outcome," ~ ",fixed_effects, "+ (1 |",clusterID,") + (1|",visit,":",clusterID,")"))
  formula_ar1 = as.formula(paste(outcome," ~ ",fixed_effects, "+ ar1(",visit,"+ 0 | ",clusterID,")"))
  
  # Fit three models
  # Exchangeable model
  fit_HH <- lmer(formula = formula_exch, data = data, REML=T)
  # Nested exchangeable model
  fit_HG <- lmer(formula = formula_nest_exch, data = data, REML=T)
  # Exponential decay model
  fit_AR <- glmmTMB(formula = formula_ar1, REML=T, 
                    data = data, family = gaussian)
  
  
  # Extract variance component estimates
  V1 <- as.data.frame(VarCorr(fit_HH))
  V2 <- as.data.frame(VarCorr(fit_HG))
  V3 <- VarCorr(fit_AR)
  
  # Estimate within-period ICC for exchangeable model
  wpicc_hh <- round(V1[,"vcov"][which(V1[,"grp"] == clusterID)]/sum(V1[,"vcov"]),4)
  
  # Estimate within-, between-period ICC and cac for nested exchangeable model
  wpicc_hg <- round((V2[,"vcov"][which(V2[,"grp"] == clusterID)] + V2[,"vcov"][which(V2[,"grp"] == paste0(visit,":",clusterID))])/sum(V2[,"vcov"]),6)
  bpicc_hg <- round((V2[,"vcov"][which(V2[,"grp"] == clusterID)])/sum(V2[,"vcov"]),4)
  cac_hg <- round(bpicc_hg/wpicc_hg,4)
  
  # Estimate within-period ICC and decay rate for exponential decay model
  wpicc_ar <- round(V3$cond$anon_villid[1,1]/sum(V3$cond$anon_villid[1,1],(attr(V3$cond,"sc"))^2),6)
  ar1 <- round(attr(V3$cond$anon_villid,"correlation")[1,2],4)
  
  
  res = data.frame(
    "Model" = c("Exch","Nested-Exch","Exp. Decay"),
    "AIC" = c(AIC(fit_HH),AIC(fit_HG),AIC(fit_AR)),
    "BIC" = c(BIC(fit_HH),BIC(fit_HG),BIC(fit_AR)),
    "Bp-ICC" = c(wpicc_hh,bpicc_hg,"-"),
    "Wp-ICC" = c(wpicc_hh,wpicc_hg,wpicc_ar),
    "CAC" = c(1,cac_hg,ar1),
    
    "Best.model.AIC" = ifelse(c(AIC(fit_HH),AIC(fit_HG),AIC(fit_AR))==min(c(AIC(fit_HH),AIC(fit_HG),AIC(fit_AR))),1,0),
    "Best.model.BIC" = ifelse(c(BIC(fit_HH),BIC(fit_HG),BIC(fit_AR))==min(c(BIC(fit_HH),BIC(fit_HG),BIC(fit_AR))),1,0)
    
  )
  
  
  res$min_AIC = res[which(res$AIC== min(res$AIC)),"AIC"]
  res$delta_AIC = round(res$AIC - res$min_AIC,1)
  
  return(list("Outcome" = outcome,"Results" = res) )
  
}




data = haven::read_dta(file = "longitudinal main set.dta")
data$sev_underweight = ifelse(data$waz_who_ < -3,1,0)


mec = "LAKANA EED mechanistic data, with EE score.csv"
mec_data = readr::read_csv(mec) %>% as.data.frame()

## Mec data into long format
mec_data$EE_1 = mec_data$V_1EE
mec_data$EE_2 = mec_data$V_2EE


mec_data$lowEE_1 = ifelse(mec_data$EE_1 < median(mec_data$EE_1,na.rm=T),1,0) %>% as.factor()
mec_data$lowEE_2 = ifelse(mec_data$EE_2 < median(mec_data$EE_2,na.rm=T),1,0) %>% as.factor()

varyingVec = c("AAT.results_1","EE_1","lowEE_1","MPO.results_1", "NEO.results_1" ,
               "AAT.results_2","EE_2","lowEE_2", "MPO.results_2","NEO.results_2") 

varyingN = c("AAT.results","EE","lowEE", "MPO.results","NEO.results")

mec_data_long = reshape(mec_data[,c("anon_id","VillageID","AAT.results_1",
                                    "AAT.results_2", "NEO.results_1",
                                    "NEO.results_2", "MPO.results_1", 
                                    "MPO.results_2", "EE_1","EE_2",
                                    "lowEE_1","lowEE_2",
                                    "interventionbinary"
)], direction="long", 
timevar="visit", varying = varyingVec, v.names=varyingN,
idvar="anon_id",ids = mec_data$anon_id)

mec_data_long$logAAT = log(as.numeric(mec_data_long$AAT.results) )
mec_data_long$logMPO = log(as.numeric(mec_data_long$MPO.results))
mec_data_long$logNEO = log(as.numeric(mec_data_long$NEO.results))



growth = "LAKANA Growth outcomes data_de-ident_HIPAA_20251003.csv" 
growth_data = readr::read_csv(growth) %>% as.data.frame()



amr_data1 = readxl::read_excel("LAKANA-AMR data_de-ident_HIPAA_2026-06-03.xlsx",
                               sheet = "E. coli",
                               guess_max = 203054) %>% as.data.frame()

amr_data1$AzithromycinEcoliR16 = as.numeric(amr_data1$AzithromycinEcoliR16)
amr_data1$AzithromycinEcoliR32 = as.numeric(amr_data1$AzithromycinEcoliR32)



varlist1 = c()
for(i in grep("RC",names(amr_data1),value=T) ) {
  
  amr_data1[,i] = as.numeric( amr_data1[,i] )
  
  varlist1 = append(varlist1,i)
  
}



amr_data1$anon_villid = amr_data1$anon_village


ecoli_azi32_res = corrStuct(data = amr_data1 ,outcome = "AzithromycinEcoliR32", 
                            clusterID = "anon_villid",visit = "Visit", covariates = c("Intervention_group"))$Results
ecoli_azi32_res$prevalence = sum(amr_data1[,"AzithromycinEcoliR32"],na.rm=T) / length(which(!is.na(amr_data1[,"AzithromycinEcoliR32"])))
ecoli_azi32_res$outcome = "Ecoli azi32"
ecoli_azi32_res$period_length = "12 months"
ecoli_azi32_res$n_timepoints = "4"



#import excel "/LAKANA-AMR data_de-ident_HIPAA_2026-06-03.xlsx", sheet("S. pneumoniae") firstrow
amr_data2 = readxl::read_excel("LAKANA-AMR data_de-ident_HIPAA_2026-06-03.xlsx",
                               sheet = "S. pneumoniae",
                               guess_max = 203054) %>% as.data.frame()


amr_data2$erythtomycin_outcome = as.numeric( ifelse(amr_data2$Erythromycin15SP<19 & 
                                                      !is.na(amr_data2$Erythromycin15SP),1,0) )


varlist2 = c()
for(i in grep("SPRC",names(amr_data2),value=T) ) {
  
  amr_data2[,i] = as.numeric( amr_data2[,i] )
  
  varlist2 = append(varlist2,i)
  
}

amr_data2$anon_villid = amr_data2$anon_village



erythromycin_res = corrStuct(data = amr_data2 ,outcome = "erythtomycin_outcome", 
                             clusterID = "anon_villid",visit = "Visit", covariates = c("Intervention_group"))$Results
erythromycin_res$prevalence = sum(amr_data2[,varlist2[3]],na.rm=T) / length(which(!is.na(amr_data2[,varlist2[3]])))
erythromycin_res$outcome = "Spneumo Erythromycin"
erythromycin_res$period_length = "12 months"
erythromycin_res$n_timepoints = "4"



# # # Results
## Growth
#WAZ
waz_res = corrStuct(data = growth_data,outcome = "WAZ", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
waz_res$prevalence = NA
waz_res$outcome = "WAZ"
waz_res$period_length = "3 months"
waz_res$n_timepoints = "4"

#LAZ
laz_res = corrStuct(data = growth_data,outcome = "LAZ", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
laz_res$prevalence = NA
laz_res$outcome = "LAZ"
laz_res$period_length = "3 months"
laz_res$n_timepoints = "4"
#WLZ
wlz_res = corrStuct(data = growth_data,outcome = "WLZ", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
wlz_res$prevalence = NA
wlz_res$outcome = "WLZ"
wlz_res$period_length = "3 months"
wlz_res$n_timepoints = "4"
#MUACZ
muacz_res = corrStuct(data = growth_data,outcome = "MUACZ", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
muacz_res$prevalence = NA
muacz_res$outcome = "MUACZ"
muacz_res$period_length = "3 months"
muacz_res$n_timepoints = "4"

growth_data$underweight = ifelse(growth_data$WAZ< -2,1,0 )
growth_data$stunted = ifelse(growth_data$LAZ< -2,1,0 )
growth_data$wasting = ifelse(growth_data$WLZ< -2,1,0 )
growth_data$lowmuac = ifelse(growth_data$MUACZ< -2,1,0 )

# Underweight
underweight_res = corrStuct(data = growth_data,outcome = "underweight", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
underweight_res$prevalence = sum(growth_data[,"underweight"],na.rm=T) / length(which(!is.na(growth_data[,"underweight"])))
underweight_res$outcome = "WAZ < -2"
underweight_res$period_length = "3 months"
underweight_res$n_timepoints = "4"

# Stunted
stunted_res = corrStuct(data = growth_data,outcome = "stunted", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
stunted_res$prevalence = sum(growth_data[,"stunted"],na.rm=T) / length(which(!is.na(growth_data[,"stunted"])))
stunted_res$outcome = "LAZ < -2"
stunted_res$period_length = "3 months"
stunted_res$n_timepoints = "4"

# Wasting
wasting_res = corrStuct(data = growth_data,outcome = "wasting", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
wasting_res$prevalence = sum(growth_data[,"wasting"],na.rm=T) / length(which(!is.na(growth_data[,"wasting"])))
wasting_res$outcome = "WLZ < -2"
wasting_res$period_length = "3 months"
wasting_res$n_timepoints = "4"

# low MUAC
lowmuac_res = corrStuct(data = growth_data,outcome = "lowmuac", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
lowmuac_res$prevalence = sum(growth_data[,"lowmuac"],na.rm=T) / length(which(!is.na(growth_data[,"lowmuac"])))
lowmuac_res$outcome = "MUACZ < -2"
lowmuac_res$period_length = "3 months"
lowmuac_res$n_timepoints = "4"




growth_data$sev_underweight = ifelse(growth_data$WAZ< -3,1,0 )
growth_data$sev_stunted = ifelse(growth_data$LAZ< -3,1,0 )
growth_data$sev_wasting = ifelse(growth_data$WLZ< -3,1,0 )
growth_data$sev_lowmuac = ifelse(growth_data$MUACZ< -3,1,0 )

# Underweight
sev_underweight_res = corrStuct(data = growth_data,outcome = "sev_underweight", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
sev_underweight_res$prevalence = sum(growth_data[,"sev_underweight"],na.rm=T) / length(which(!is.na(growth_data[,"sev_underweight"])))
sev_underweight_res$outcome = "WAZ < -3"
sev_underweight_res$period_length = "3 months"
sev_underweight_res$n_timepoints = "4"

# Stunted
sev_stunted_res = corrStuct(data = growth_data,outcome = "sev_stunted", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
sev_stunted_res$prevalence = sum(growth_data[,"sev_stunted"],na.rm=T) / length(which(!is.na(growth_data[,"sev_stunted"])))
sev_stunted_res$outcome = "LAZ < -3"
sev_stunted_res$period_length = "3 months"
sev_stunted_res$n_timepoints = "4"

# Wasting
sev_wasting_res = corrStuct(data = growth_data,outcome = "sev_wasting", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
sev_wasting_res$prevalence = sum(growth_data[,"sev_wasting"],na.rm=T) / length(which(!is.na(growth_data[,"sev_wasting"])))
sev_wasting_res$outcome = "WLZ < -3"
sev_wasting_res$period_length = "3 months"
sev_wasting_res$n_timepoints = "4"

# low MUAC
sev_lowmuac_res = corrStuct(data = growth_data,outcome = "sev_lowmuac", clusterID = "anon_villid",visit = "visit", covariates = c("intv"))$Results
sev_lowmuac_res$prevalence = sum(growth_data[,"sev_lowmuac"],na.rm=T) / length(which(!is.na(growth_data[,"sev_lowmuac"])))
sev_lowmuac_res$outcome = "MUACZ < -3"
sev_lowmuac_res$period_length = "3 months"
sev_lowmuac_res$n_timepoints = "4"




## Main Trial
# Mortality
data$anon_villid = data$villageid
mort_res = corrStuct(data = data,outcome = "death_", clusterID = "anon_villid",visit = "visit", covariates = c("intv","calendarmonth_"))$Results
mort_res$prevalence = sum(data[,"death_"],na.rm=T) / length(which(!is.na(data[,"death_"])))
mort_res$outcome = "mortality"
mort_res$period_length = "3 months"
mort_res$n_timepoints = "8"
# WAZ 
main_waz_res = corrStuct(data = data,outcome = "waz_who_", clusterID = "anon_villid",visit = "visit", covariates = c("intv","calendarmonth_"))$Results
main_waz_res$prevalence = NA
main_waz_res$outcome = "main trial WAZ"
main_waz_res$period_length = "3 months"
main_waz_res$n_timepoints = "8"



# Underweight 
main_underweight_res = corrStuct(data = data,outcome = "underweight_", clusterID = "anon_villid",visit = "visit", covariates = c("intv","calendarmonth_"))$Results
main_underweight_res$prevalence = sum(data[,"underweight_"],na.rm=T) / length(which(!is.na(data[,"underweight_"])))
main_underweight_res$outcome = "main trial underweight"
main_underweight_res$period_length = "3 months"
main_underweight_res$n_timepoints = "8"

main_sev_underweight_res = corrStuct(data = data,outcome = "sev_underweight", clusterID = "anon_villid",visit = "visit", covariates = c("intv","calendarmonth_"))$Results
main_sev_underweight_res$prevalence = sum(data[,"sev_underweight"],na.rm=T) / length(which(!is.na(data[,"sev_underweight"])))
main_sev_underweight_res$outcome = "main trial severe underweight"
main_sev_underweight_res$period_length = "3 months"
main_sev_underweight_res$n_timepoints = "8"


mec_data_long$anon_villid = mec_data_long$VillageID
mec_data_long$logMPO = ifelse(is.infinite(mec_data_long$logMPO),NA,mec_data_long$logMPO)
mec_data_long$logNEO = ifelse(is.infinite(mec_data_long$logNEO),NA,mec_data_long$logNEO)
# AAT
aat_res = corrStuct(data = mec_data_long,outcome = "logAAT", clusterID = "anon_villid",visit = "visit", covariates = c("interventionbinary"))$Results
aat_res$prevalence = NA
aat_res$outcome = "AAT"
aat_res$period_length = "2 weeks"
aat_res$n_timepoints = "2"
# MPO 
mpo_res = corrStuct(data = mec_data_long,outcome = "logMPO", clusterID = "anon_villid",visit = "visit", covariates = c("interventionbinary"))$Results
mpo_res$prevalence = NA
mpo_res$outcome = "MPO"
mpo_res$period_length = "2 weeks"
mpo_res$n_timepoints = "2"
# NEO 
neo_res = corrStuct(data = mec_data_long,outcome = "logNEO", clusterID = "anon_villid",visit = "visit", covariates = c("interventionbinary"))$Results
neo_res$prevalence = NA
neo_res$outcome = "NEO"
neo_res$period_length = "2 weeks"
neo_res$n_timepoints = "2"
# EE
mec_data_long$logEE = log(mec_data_long$EE)
mec_data_long$logEE = ifelse(is.infinite(mec_data_long$logEE),NA,mec_data_long$logEE)
ee_res = corrStuct(data = mec_data_long,outcome = "logEE", clusterID = "anon_villid",visit = "visit", covariates = c("interventionbinary"))$Results
ee_res$prevalence = NA
ee_res$outcome = "EE"
ee_res$period_length = "2 weeks"
ee_res$n_timepoints = "2"
# Low EE
mec_data_long$lowEE = as.numeric(mec_data_long$lowEE)-1
ee_bin_res = corrStuct(data = mec_data_long,outcome = "lowEE", clusterID = "anon_villid",visit = "visit", covariates = c("interventionbinary"))$Results
ee_bin_res$prevalence = sum(mec_data_long[,"lowEE"],na.rm=T) / length(which(!is.na(mec_data_long[,"lowEE"])))
ee_bin_res$outcome = "Low EE"
ee_bin_res$period_length = "2 weeks"
ee_bin_res$n_timepoints = "2"




pause = matrix(ncol=ncol(ee_res),nrow=1) %>% as.data.frame()
names(pause) = names(ee_res)
pause = ifelse(is.na(pause)," ",pause)


c(which(mort_res$Best.model.AIC==1),
  which(main_waz_res$Best.model.AIC==1),
  which(main_underweight_res$Best.model.AIC==1),
  which(main_underweight_res$Best.model.AIC==1),
  which(waz_res$Best.model.AIC==1),
  which(laz_res$Best.model.AIC==1),
  which(wlz_res$Best.model.AIC==1),
  which(muacz_res$Best.model.AIC==1),
  which(underweight_res$Best.model.AIC==1),
  which(stunted_res$Best.model.AIC==1),
  which(wasting_res$Best.model.AIC==1),
  which(lowmuac_res$Best.model.AIC==1),
  which(sev_underweight_res$Best.model.AIC==1),
  which(sev_stunted_res$Best.model.AIC==1),
  which(sev_wasting_res$Best.model.AIC==1),
  which(sev_lowmuac_res$Best.model.AIC==1),
  which(aat_res$Best.model.AIC==1),
  which(mpo_res$Best.model.AIC==1),
  which(neo_res$Best.model.AIC==1),
  which(ee_res$Best.model.AIC==1),
  which(ee_bin_res$Best.model.AIC==1)) %>% table()




full = rbind.data.frame(mort_res,pause,main_waz_res,pause,main_underweight_res,pause,main_sev_underweight_res,pause,
                        laz_res,pause,wlz_res,pause,muacz_res,pause,
                        stunted_res,pause,wasting_res,pause,lowmuac_res,pause,
                        sev_stunted_res,pause,sev_wasting_res,pause,sev_lowmuac_res,pause,
                        aat_res,pause,mpo_res,pause,neo_res,pause,ee_res,pause,ee_bin_res, pause,
                        ecoli_azi32_res,pause, erythromycin_res,pause
)



# Treatment Arms
L = 2 #c(2:8)

# Times
t = 4 #c(2:20)

# Number of individuals
m = 10


A <- matrix(c(
  1,1,1,1,
  0,0,0,0
), nrow = 2, byrow = TRUE)

B <- sum(A)
C <- sum(rowSums(A)^2)
D <- sum(colSums(A)^2)

rho = as.numeric(full$Wp.ICC) #c(wpicc_hg,wpicc_ar)
cac = as.numeric(full$CAC) #c(cac_hg,ar1)

totalCorr = (m*rho*cac) / (1 + (m-1)*rho)
deff_r_r = ((L^2)*(1-totalCorr)*(1 + t*totalCorr)) / 
  4*(L*B - D + ((B^2) + L*t*B - t*D - L*C)*totalCorr)



A <- rbind(c(0, 1, 1, 1),
           c(0, 0, 1, 1),
           c(0, 0, 0, 1))

# B: sum of all entries
B <- sum(A)

C <- sum(rowSums(A)^2)

D <- sum(colSums(A)^2)


L=nrow(A);
Tee=L


r = (m*rho*cac) / (1+(m-1)*rho)
deff_r = ( 3*L*(1-r)*(1+L*r))/((L-1)*(2+L*r)) 






deff_c = 1 + ( m-1 ) * rho

full$Deff_c = deff_c
full$Deff_r = deff_r
full$"Deff_r x Deff_c" = deff_r * deff_c



par(mfrow=c(2,1))

x <- log(c(.1/100,1/100, 10/100, 40/100,50/100))
y_mod = exp(-2.80+x*0.91)


y_mod_lb = exp(-5.20+x*0.91)
y_mod_ub = exp(0.20+x*0.91)


# Make an empty log-log plot
plot(exp(x)*100, y_mod,main = "HTA context:\n ICC model estimates ",
     log = "xy",
     type = "n",     # empty plot
     xlab = "Prevalence (%)",
     ylab = "Wp-ICC",
     xlim = c(0.01, 100),
     ylim = c(1.0E-06, 1),
     cex.lab=1.7,
     cex.axis = 1.7,
     cex.main = 1.7)

# Add diagonal line y = x
lines(exp(x)*100, y_mod, col = "black", lwd = 2,lty=2)

lines(exp(x)*100, y_mod_lb, col = "grey", lwd = 2,lty=2)
lines(exp(x)*100, y_mod_ub, col = "grey", lwd = 2,lty=2)


legend("topleft",
       legend = c("1-11-month-old mortality","WAZ < -2","WAZ < -3","LAZ < -2","LAZ < -3",
                  "WLZ < -2","WLZ < -3","MUAC-Z < -2","S. Pneumoaie","E. Coli (MIC ≤ 32)","Below median EE score"),
       pch = c(1,2,3,4,5,6,7,8,10,12,13),
       col=1)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="mortality" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="mortality" ),]$Wp.ICC),cex = 1.8, pch = 1, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial underweight" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial underweight" ),]$Wp.ICC),cex = 1.8, pch = 2, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial severe underweight" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial severe underweight" ),]$Wp.ICC),cex = 1.8, pch = 3, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 4, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -3" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 5, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 6, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -3" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 7, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 8, col = 1)

# points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -3" ),]$prevalence)*100,
#        y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 9, col = 1)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Spneumo Erythromycin" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Spneumo Erythromycin" ),]$Wp.ICC),cex = 1.8, pch = 10, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Ecoli azi32" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Ecoli azi32" ),]$Wp.ICC),cex = 1.8, pch = 12, col = 1)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Low EE" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Low EE" ),]$Wp.ICC),cex = 1.8, pch = 13, col = 1)






x <- log(c(.001/100,.01/100,.1/100,1/100, 10/100, 40/100, 50/100))
y_mod = exp(-2.10+x*0.61)


y_mod_lb = exp(-5.82+x*0.61)
y_mod_ub = exp(0.46+x*0.61)

# Make an empty log-log plot
plot(exp(x)*100, y_mod, main = "GPRD context:\n ICC model estimates",
     log = "xy",
     type = "n",     # empty plot
     xlab = "Prevalence (%)",
     ylab = "Wp-ICC",
     xlim = c(0.0001, 100),
     ylim = c(1.0E-06, 1),
     cex.lab = 1.7,
     cex.axis = 1.7,
     cex.main = 1.7)

# Add diagonal line y = x
lines(exp(x)*100, y_mod, col = "black", lwd = 2,lty=2)

lines(exp(x)*100, y_mod_lb, col = "grey", lwd = 2,lty=2)
lines(exp(x)*100, y_mod_ub, col = "grey", lwd = 2,lty=2)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="mortality" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="mortality" ),]$Wp.ICC),cex = 1.8, pch = 1, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial underweight" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial underweight" ),]$Wp.ICC),cex = 1.8, pch = 2, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial severe underweight" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="main trial severe underweight" ),]$Wp.ICC),cex = 1.8, pch = 3, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 4, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -3" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="LAZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 5, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 6, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -3" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="WLZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 7, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -2" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -2" ),]$Wp.ICC),cex = 1.8, pch = 8, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -3" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="MUACZ < -3" ),]$Wp.ICC),cex = 1.8, pch = 9, col = 1)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Spneumo Erythromycin" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Spneumo Erythromycin" ),]$Wp.ICC),cex = 1.8, pch = 10, col = 1)

points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Ecoli azi32" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Ecoli azi32" ),]$Wp.ICC),cex = 1.8, pch = 12, col = 1)


points(x=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Low EE" ),]$prevalence)*100,
       y=as.numeric(full[which(full$Model == "Exch" & full$outcome=="Low EE" ),]$Wp.ICC),cex = 1.8, pch = 13, col = 1)





openxlsx::write.xlsx(full,"Results.xlsx",overwrite = T)


