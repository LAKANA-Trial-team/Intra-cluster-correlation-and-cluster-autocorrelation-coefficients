

## Simulation

### Simulation result inspection and tabulation

library(tidyverse)
library(plotrix)
library(Hmisc)


#Replace with your own paths
icc_no_intv_path = "/ICC no intervention results"
icc_intv_path = "/ICC with intervention results"




low_icc_noIntv = haven::read_dta(paste0(icc_no_intv_path,"/low_onegroup ICC results_no intv.dta"))
hi_icc_noIntv = haven::read_dta(paste0(icc_no_intv_path,"/onegroup ICC results_no intv.dta"))



low_icc_Intv = haven::read_dta(paste0(icc_intv_path,"/low_ICC results_intv.dta"))
hi_icc_Intv = haven::read_dta(paste0(icc_intv_path,"/ICC results_intv.dta"))




agg = aggregate(low_icc_noIntv,by=list(low_icc_noIntv$scenario),FUN = function(x) c("mean" = mean(x,na.rm=T),
                                                                      "sd" = sqrt(var(x,na.rm=T)),
                                                                      "n" = length(na.omit(x)),
                                                                      "2.5%" = quantile(na.omit(x), c(.025)),
                                                                      "97.5%" = quantile(na.omit(x), c(.975))) )

par(mar = c(5, 7, 4, 2)) 
plot(agg$Group.1,agg$identityICC[,1],ylim = c(-.0,.42), pch = 18,xlim=c(0.5,11),
     xlab = substitute(paste(bold("Simulation scenario"))), 
     ylab = substitute(paste(bold("Mean ICC estimate (Monte Carlo 95% CI )"))),
     main = "Simulation results, 1000 iteration rounds / scenario", 
     xaxt = "n",cex.lab=2,cex.main=2,cex.axis = 2)

abline(h = 0.01,lty=2)
segments(agg$Group.1,agg$identityICC[,4],agg$Group.1,agg$identityICC[,5],col=1)




#######################
#######################
#######################
#######################


res_low_icc_noIntv = low_icc_noIntv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(identityICC,na.rm=T), 
                   n_na = sum(is.na(identityICC)),
                   "2.5%" = quantile(na.omit(identityICC), c(.025)),
                   "97.5%" = quantile(na.omit(identityICC), c(.975)),
                   "Q1" = quantile(na.omit(identityICC), c(.25)),
                   "Q3" = quantile(na.omit(identityICC), c(.75))) %>% 
  as.data.frame() %>% na.omit()



res_low_icc_noIntv$target_ICC = 0.01
res_low_icc_noIntv$diff = res_low_icc_noIntv$mediaani-res_low_icc_noIntv$target_ICC
res_low_icc_noIntv$rel_bias = paste0( round(((res_low_icc_noIntv$mediaani-res_low_icc_noIntv$target_ICC)/res_low_icc_noIntv$target_ICC) *100), " %")
res_low_icc_noIntv$n_iter = 1000
res_low_icc_noIntv$conv_perc = paste0(round((1-res_low_icc_noIntv$n_na/res_low_icc_noIntv$n_iter)*100,1), " %")
res_low_icc_noIntv$m = c(50,5,50,50,50,
                         50,5,50,50,50,
                         50,5,50,50,50)
res_low_icc_noIntv$nClus = c(100,100,100,100,100,
                         50,50,50,50,50,
                         10,10,10,10,10)
res_low_icc_noIntv$prev = c(0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02)
res_low_icc_noIntv$trtEff = 0
res_low_icc_noIntv$scenario = c(1,2,3,4,5,
                                1,2,3,4,5,
                                1,2,3,4,5)
a = res_low_icc_noIntv[,c("target_ICC","nClus","m","prev","mediaani","Q1","Q3","diff","rel_bias","conv_perc","scenario")]



#######################
#######################
#######################
#######################

res_hi_icc_noIntv = hi_icc_noIntv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(identityICC,na.rm=T), 
                   n_na = sum(is.na(identityICC)),
                   "2.5%" = quantile(na.omit(identityICC), c(.025)),
                   "97.5%" = quantile(na.omit(identityICC), c(.975)),
                   "Q1" = quantile(na.omit(identityICC), c(.25)),
                   "Q3" = quantile(na.omit(identityICC), c(.75))) %>% 
  as.data.frame() %>% na.omit()

res_hi_icc_noIntv$target_ICC = 0.1
res_hi_icc_noIntv$diff = res_hi_icc_noIntv$mediaani-res_hi_icc_noIntv$target_ICC
res_hi_icc_noIntv$rel_bias = paste0( round(((res_hi_icc_noIntv$mediaani-res_hi_icc_noIntv$target_ICC)/res_hi_icc_noIntv$target_ICC) *100), " %")
res_hi_icc_noIntv$n_iter = 1000
res_hi_icc_noIntv$conv_perc = paste0(round((1-res_hi_icc_noIntv$n_na/res_hi_icc_noIntv$n_iter)*100,1), " %")
res_hi_icc_noIntv$m = c(50,5,50,50,50,
                         50,5,50,50,50,
                         50,5,50,50,50)
res_hi_icc_noIntv$nClus = c(100,100,100,100,100,
                             50,50,50,50,50,
                             10,10,10,10,10)
res_hi_icc_noIntv$prev = c(0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02)
res_hi_icc_noIntv$trtEff = 0
res_hi_icc_noIntv$scenario = c(1,2,3,4,5,
                                1,2,3,4,5,
                                1,2,3,4,5)
b = res_hi_icc_noIntv[,c("target_ICC","nClus","m","prev","mediaani","Q1","Q3","diff","rel_bias","conv_perc","scenario")]




###########################
###########################
###########################
###########################

######################################################
######################################################
######################################################
######################################################

###########################
###########################
###########################
###########################


res_low_icc_Intv = low_icc_Intv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(identityICC,na.rm=T), 
                   n_na = sum(is.na(identityICC)),
                   "2.5%" = quantile(na.omit(identityICC), c(.025)),
                   "97.5%" = quantile(na.omit(identityICC), c(.975)),
                   "Q1" = quantile(na.omit(identityICC), c(.25)),
                   "Q3" = quantile(na.omit(identityICC), c(.75))) %>% 
  as.data.frame() %>% na.omit()

res_low_icc_Intv$target_ICC = 0.01
res_low_icc_Intv$diff = res_low_icc_Intv$mediaani-res_low_icc_Intv$target_ICC
res_low_icc_Intv$rel_bias = paste0( round(((res_low_icc_Intv$mediaani-res_low_icc_Intv$target_ICC)/res_low_icc_Intv$target_ICC) *100), " %")
res_low_icc_Intv$n_iter = 1000
res_low_icc_Intv$conv_perc = paste0(round((1-res_low_icc_Intv$n_na/res_low_icc_Intv$n_iter)*100,1), " %")
res_low_icc_Intv$m = c(50,5,50,50,50,
                      50,5,50,50,50,
                      50,5,50,50,50)
res_low_icc_Intv$nClus = c(100,100,100,100,100,
                          50,50,50,50,50,
                          10,10,10,10,10)
res_low_icc_Intv$prev = c(0.5,0.5,0.9,0.2,0.02,
                         0.5,0.5,0.9,0.2,0.02,
                         0.5,0.5,0.9,0.2,0.02)
res_low_icc_Intv$trtEff = 0.1
res_low_icc_Intv$scenario = c(1,2,3,4,5,
                             1,2,3,4,5,
                             1,2,3,4,5)
c = res_low_icc_Intv[,c("target_ICC","nClus","m","prev","mediaani","Q1","Q3","diff","rel_bias","conv_perc","scenario")]

###########################
###########################
###########################
###########################


res_hi_icc_Intv = hi_icc_Intv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(identityICC,na.rm=T), 
                   n_na = sum(is.na(identityICC)),
                   "2.5%" = quantile(na.omit(identityICC), c(.025)),
                   "97.5%" = quantile(na.omit(identityICC), c(.975)),
                   "Q1" = quantile(na.omit(identityICC), c(.25)),
                   "Q3" = quantile(na.omit(identityICC), c(.75))) %>% 
  as.data.frame() %>% na.omit()

res_hi_icc_Intv$target_ICC = 0.1
res_hi_icc_Intv$diff = res_hi_icc_Intv$mediaani-res_hi_icc_Intv$target_ICC
res_hi_icc_Intv$rel_bias = paste0( round(((res_hi_icc_Intv$mediaani-res_hi_icc_Intv$target_ICC)/res_hi_icc_Intv$target_ICC) *100), " %")
res_hi_icc_Intv$n_iter = 1000
res_hi_icc_Intv$conv_perc = paste0(round((1-res_hi_icc_Intv$n_na/res_hi_icc_Intv$n_iter)*100,1), " %")
res_hi_icc_Intv$m = c(50,5,50,50,50,
                        50,5,50,50,50,
                        50,5,50,50,50)
res_hi_icc_Intv$nClus = c(100,100,100,100,100,
                            50,50,50,50,50,
                            10,10,10,10,10)
res_hi_icc_Intv$prev = c(0.5,0.5,0.9,0.2,0.02,
                           0.5,0.5,0.9,0.2,0.02,
                           0.5,0.5,0.9,0.2,0.02)
res_hi_icc_Intv$trtEff = 0.1
res_hi_icc_Intv$scenario = c(1,2,3,4,5,
                               1,2,3,4,5,
                               1,2,3,4,5)
d = res_hi_icc_Intv[,c("target_ICC","nClus","m","prev","mediaani","Q1","Q3","diff","rel_bias","conv_perc","scenario")]



ab = rbind.data.frame(a,b)
cd = rbind.data.frame(c,d)


ab$Q1_rel_bias = paste0( round(((ab$Q1-ab$target_ICC)/ab$target_ICC) *100), " %")
ab$Q3_rel_bias = paste0( round(((ab$Q3-ab$target_ICC)/ab$target_ICC) *100), " %")

cd$Q1_rel_bias = paste0( round(((cd$Q1-cd$target_ICC)/cd$target_ICC) *100), " %")
cd$Q3_rel_bias = paste0( round(((cd$Q3-cd$target_ICC)/cd$target_ICC) *100), " %")



ab$"mediaani (Q1, Q3)" = paste0(round(ab$mediaani,3)," (",round(ab$Q1,4),", ",round(ab$Q3,3),")" )
ab$"Bias (Q1, Q3)" = paste0(ab$rel_bias," (",ab$Q1_rel_bias,", ",ab$Q3_rel_bias,")" )

cd$"mediaani (Q1, Q3)" = paste0(round(cd$mediaani,3)," (",round(cd$Q1,4),", ",round(cd$Q3,3),")" )
cd$"Bias (Q1, Q3)" = paste0(cd$rel_bias," (",cd$Q1_rel_bias,", ",cd$Q3_rel_bias,")" )




agg1 = aggregate(low_icc_noIntv,by=list(low_icc_noIntv$scenario),FUN = function(x) c("mean" = mean(x,na.rm=T),
                                                                                     "sd" = sqrt(var(x,na.rm=T)),
                                                                                     "n" = length(na.omit(x)),
                                                                                     "2.5%" = quantile(na.omit(x), c(.025)),
                                                                                     "97.5%" = quantile(na.omit(x), c(.975))) )



agg2 = aggregate(hi_icc_noIntv,by=list(hi_icc_noIntv$scenario),FUN = function(x) c("mean" = mean(x,na.rm=T),
                                                                                   "sd" = sqrt(var(x,na.rm=T)),
                                                                                   "n" = length(na.omit(x)),
                                                                                   "2.5%" = quantile(na.omit(x), c(.025)),
                                                                                   "97.5%" = quantile(na.omit(x), c(.975))) )






labels <- c("n = 100, \np = 0.5, \nm = 50",
            "n = 100, \np = 0.5, \nm = 5",
            "n = 100, \np = 0.9, \nm = 50",
            "n = 100, \np = 0.2, \nm = 50",
            "n = 100, \np = 0.02, \nm = 50",
            
            "n = 50, \np = 0.5, \nm = 50",
            "n = 50, \np = 0.5, \nm = 5",
            "n = 50, \np = 0.9, \nm = 50",
            "n = 50, \np = 0.2, \nm = 50",
            "n = 50, \np = 0.02, \nm = 50",
            
            "n = 10, \np = 0.5, \nm = 50",
            "n = 10, \np = 0.5, \nm = 5",
            "n = 10, \np = 0.9, \nm = 50",
            "n = 10, \np = 0.2, \nm = 50",
            "n = 10, \np = 0.02, \nm = 50"
            )


par(mar = c(5, 7, 4, 2), mfrow=c(2,1)) 

plot(agg1$Group.1,agg1$identityICC[,1],ylim = c(-.0,0.5), pch = 18,xlim=c(0.5,15),
     xlab = substitute(paste(bold(""))), 
     ylab = substitute(paste(bold("Mean ICC estimate \n(Monte Carlo 95% CI )"))),
     main = "Simulation results, 1000 iteration rounds / scenario", 
     xaxt="n",cex.lab=2,cex.main=2,cex.axis = 2)
axis(1, at = agg$Group.1)
abline(h = 0.01,lty=2)
segments(agg1$Group.1,agg1$identityICC[,4],agg1$Group.1,agg1$identityICC[,5],col=1)
legend("topleft",legend = c("Target ICC"), lty=c(2),
       cex=1.4)


text(x = agg1$Group.1+0.05,
     y = 0.4,
     labels = labels,
     srt = 45,          
     adj = 1,           
     xpd = TRUE,
     cex=.88
)



labels <- c("n = 100, \np = 0.5, \nm = 50",
            "n = 100, \np = 0.5, \nm = 5",
            "n = 100, \np = 0.9, \nm = 50",
            "n = 100, \np = 0.2, \nm = 50",
            "n = 100, \np = 0.02, \nm = 50",
            
            "n = 50, \np = 0.5, \nm = 50",
            "n = 50, \np = 0.5, \nm = 5",
            "n = 50, \np = 0.9, \nm = 50",
            "n = 50, \np = 0.2, \nm = 50",
            "n = 50, \np = 0.02, \nm = 50",
            
            "n = 10, \np = 0.5, \nm = 50",
            "n = 10, \np = 0.5, \nm = 5",
            "n = 10, \np = 0.9, \nm = 50",
            "n = 10, \np = 0.2, \nm = 50",
            "n = 10, \np = 0.02, \nm = 50"
)
plot(agg2$Group.1,agg2$identityICC[,1],ylim = c(-.0,0.8), pch = 18,xlim=c(0.5,15),
     xlab = substitute(paste(bold("Simulation Scenario"))), 
     ylab = substitute(paste(bold("Mean ICC estimate \n(Monte Carlo 95% CI )"))),
     main = "", 
     xaxt = "n",cex.lab=2,cex.main=2,cex.axis = 2)
axis(1, at = agg$Group.1)
abline(h = 0.1,lty=2)
segments(agg2$Group.1,agg2$identityICC[,4],agg2$Group.1,agg2$identityICC[,5],col=1)
legend("topleft",legend = c("Target ICC"), lty=c(2),
       cex=1.4)

text(x = agg2$Group.1+0.05,
     y = 0.6,
     labels = labels,
     srt = 45,                 
     adj = 1,                  
     xpd = TRUE,
     cex=.88
)



abcd = rbind.data.frame(ab,cd)

median(
  abs(as.numeric(gsub("%","",abcd$rel_bias))));

quantile(
  abs(as.numeric(gsub("%","",abcd$rel_bias))),probs=.25);

quantile(
  abs(as.numeric(gsub("%","",abcd$rel_bias))),probs=.75)




#openxlsx::write.xlsx(ab,"ICC sim res, no intv 2026-07-10.xlsx")
#openxlsx::write.xlsx(cd,"ICC sim res, intv 2026-07-10.xlsx")

median(as.numeric(gsub("%","",ab$rel_bias)));
quantile(as.numeric(gsub("%","",ab$rel_bias)),probs=.25)
quantile(as.numeric(gsub("%","",ab$rel_bias)),probs=.75)

median(as.numeric(gsub("%","",cd$rel_bias)));
quantile(as.numeric(gsub("%","",cd$rel_bias)),probs=.25)
quantile(as.numeric(gsub("%","",cd$rel_bias)),probs=.75)







