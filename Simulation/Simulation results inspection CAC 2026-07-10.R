

## Simulation

### Simulation result inspection and tabulation

library(tidyverse)
library(plotrix)
library(Hmisc)

# Replace with your own paths
cac_no_intv_path_hi = "nointv/CAC = 0.8"
cac_no_intv_path_lo = "nointv/CAC = 0.3"
cac_intv_path_hi = "intv/CAC = 0.8"
cac_intv_path_lo = "intv/CAC = 0.3"




hi_cac_noIntv = haven::read_dta(paste0(cac_no_intv_path_hi,"/low_CAC results_intv.dta"))
low_cac_noIntv = haven::read_dta(paste0(cac_no_intv_path_lo,"/low_CAC results_intv.dta"))



low_cac_Intv = haven::read_dta(paste0(cac_intv_path_lo,"/low_CAC results_w_intv.dta"))
hi_cac_Intv = haven::read_dta(paste0(cac_intv_path_hi,"/low_CAC results_w_intv.dta"))




agg1 = aggregate(low_cac_noIntv,by=list(low_cac_noIntv$scenario),FUN = function(x) c("mean" = mean(x,na.rm=T),
                                                                                    "sd" = sqrt(var(x,na.rm=T)),
                                                                                    "n" = length(na.omit(x)),
                                                                                    "2.5%" = quantile(na.omit(x), c(.025)),
                                                                                    "97.5%" = quantile(na.omit(x), c(.975))) )



agg2 = aggregate(hi_cac_noIntv,by=list(hi_cac_noIntv$scenario),FUN = function(x) c("mean" = mean(x,na.rm=T),
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

plot(agg1$Group.1,agg1$CAC[,1],ylim = c(-.0,2.0), pch = 18,xlim=c(0.5,15),
     xlab = substitute(paste(bold(""))), 
     ylab = substitute(paste(bold("Mean CAC estimate \n(Monte Carlo 95% CI )"))),
     main = "Simulation results, 1000 iteration rounds / scenario", 
     xaxt="n",cex.lab=2,cex.main=2,cex.axis = 2, yaxt="n")
axis(1, at = agg1$Group.1)
axis(2, at = seq(0,1,by=0.1),cex.axis=1.5)
abline(h = 0.3,lty=2)
segments(agg1$Group.1,agg1$CAC[,4],agg1$Group.1,agg1$CAC[,5],col=1)
legend("topleft",legend = c("Target CAC"), lty=c(2),
       cex=1.4)


text(x = agg1$Group.1+0.05,
     y = 1.5,
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
plot(agg2$Group.1,agg2$CAC[,1],ylim = c(-.0,2.0), pch = 18,xlim=c(0.5,15),
     xlab = substitute(paste(bold("Simulation Scenario"))), 
     ylab = substitute(paste(bold("Mean CAC estimate \n(Monte Carlo 95% CI )"))),
     main = "", 
     xaxt = "n",cex.lab=2,cex.main=2,cex.axis = 2,yaxt="n")
axis(1, at = agg2$Group.1)
axis(2, at = seq(0,1,by=0.1),cex.axis=1.5)
abline(h = 0.8,lty=2)
segments(agg2$Group.1,agg2$CAC[,4],agg2$Group.1,agg2$CAC[,5],col=1)
legend("topleft",legend = c("Target CAC"), lty=c(2),
       cex=1.4)

text(x = agg2$Group.1+0.05,
     y = 1.5,
     labels = labels,
     srt = 45,                 
     adj = 1,                  
     xpd = TRUE,
     cex=.88
)




#######################
#######################
#######################
#######################


res_low_cac_noIntv = low_cac_noIntv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(CAC,na.rm=T), 
                   targetCAC = median(targetCAC,na.rm=T), 
                   n_na = sum(is.na(CAC)),
                   "2.5%" = quantile(na.omit(CAC), c(.025)),
                   "97.5%" = quantile(na.omit(CAC), c(.975)),
                   "Q1" = quantile(na.omit(CAC), c(.25)),
                   "Q3" = quantile(na.omit(CAC), c(.75)) ) %>% 
  as.data.frame() %>% na.omit()




res_low_cac_noIntv$diff = res_low_cac_noIntv$mediaani-res_low_cac_noIntv$targetCAC
res_low_cac_noIntv$rel_bias = paste0( round(((res_low_cac_noIntv$mediaani-res_low_cac_noIntv$targetCAC)/res_low_cac_noIntv$targetCAC) *100), " %")

res_low_cac_noIntv$Q1_rel_bias = paste0( round(((res_low_cac_noIntv$Q1-res_low_cac_noIntv$targetCAC)/res_low_cac_noIntv$targetCAC) *100), " %")
res_low_cac_noIntv$Q3_rel_bias = paste0( round(((res_low_cac_noIntv$Q3-res_low_cac_noIntv$targetCAC)/res_low_cac_noIntv$targetCAC) *100), " %")

res_low_cac_noIntv$n_iter = 1000
res_low_cac_noIntv$conv_perc = paste0(round((1-res_low_cac_noIntv$n_na/res_low_cac_noIntv$n_iter)*100,1), " %")
res_low_cac_noIntv$m = c(50,5,50,50,50,
                         50,5,50,50,50,
                         50,5,50,50,50)
res_low_cac_noIntv$nClus = c(100,100,100,100,100,
                             50,50,50,50,50,
                             10,10,10,10,10)
res_low_cac_noIntv$prev = c(0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02,
                            0.5,0.5,0.9,0.2,0.02)
res_low_cac_noIntv$trtEff = 0
res_low_cac_noIntv$scenario = c(1,2,3,4,5,
                                1,2,3,4,5,
                                1,2,3,4,5)
a = res_low_cac_noIntv[,c("targetCAC","nClus","m","prev","mediaani","Q1","Q3",
                          "Q1_rel_bias","Q3_rel_bias",
                          "diff","rel_bias","conv_perc","scenario")]

a$"mediaani (Q1, Q3)" = paste0(round(a$mediaani,4)," (",round(a$Q1,4),", ",round(a$Q3,4),")" )
a$"Bias (Q1, Q3)" = paste0(a$rel_bias," (",a$Q1_rel_bias,", ",a$Q3_rel_bias,")" )


#######################
#######################
#######################
#######################

res_hi_cac_noIntv = hi_cac_noIntv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(CAC,na.rm=T), 
                   targetCAC = median(targetCAC,na.rm=T), 
                   n_na = sum(is.na(CAC)),
                   "2.5%" = quantile(na.omit(CAC), c(.025)),
                   "97.5%" = quantile(na.omit(CAC), c(.975)),
                   "Q1" = quantile(na.omit(CAC), c(.25)),
                   "Q3" = quantile(na.omit(CAC), c(.75))) %>% 
  as.data.frame() %>% na.omit()


res_hi_cac_noIntv$diff = res_hi_cac_noIntv$mediaani-res_hi_cac_noIntv$targetCAC
res_hi_cac_noIntv$rel_bias = paste0( round(((res_hi_cac_noIntv$mediaani-res_hi_cac_noIntv$targetCAC)/res_hi_cac_noIntv$targetCAC) *100), " %")

res_hi_cac_noIntv$Q1_rel_bias = paste0( round(((res_hi_cac_noIntv$Q1-res_hi_cac_noIntv$targetCAC)/res_hi_cac_noIntv$targetCAC) *100), " %")
res_hi_cac_noIntv$Q3_rel_bias = paste0( round(((res_hi_cac_noIntv$Q3-res_hi_cac_noIntv$targetCAC)/res_hi_cac_noIntv$targetCAC) *100), " %")

res_hi_cac_noIntv$n_iter = 1000
res_hi_cac_noIntv$conv_perc = paste0(round((1-res_hi_cac_noIntv$n_na/res_hi_cac_noIntv$n_iter)*100,1), " %")
res_hi_cac_noIntv$m = c(50,5,50,50,50,
                        50,5,50,50,50,
                        50,5,50,50,50)
res_hi_cac_noIntv$nClus = c(100,100,100,100,100,
                            50,50,50,50,50,
                            10,10,10,10,10)
res_hi_cac_noIntv$prev = c(0.5,0.5,0.9,0.2,0.02,
                           0.5,0.5,0.9,0.2,0.02,
                           0.5,0.5,0.9,0.2,0.02)
res_hi_cac_noIntv$trtEff = 0
res_hi_cac_noIntv$scenario = c(1,2,3,4,5,
                               1,2,3,4,5,
                               1,2,3,4,5)
b = res_hi_cac_noIntv[,c("targetCAC","nClus","m","prev","mediaani","Q1","Q3",
                         "Q1_rel_bias","Q3_rel_bias",
                         "diff","rel_bias","conv_perc","scenario")]


b$"mediaani (Q1, Q3)" = paste0(round(b$mediaani,4)," (",round(b$Q1,4),", ",round(b$Q3,4),")" )
b$"Bias (Q1, Q3)" = paste0(b$rel_bias," (",b$Q1_rel_bias,", ",b$Q3_rel_bias,")" )

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


res_low_cac_Intv = low_cac_Intv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(CAC,na.rm=T), 
                   targetCAC = median(targetCAC,na.rm=T), 
                   n_na = sum(is.na(CAC)),
                   "2.5%" = quantile(na.omit(CAC), c(.025)),
                   "97.5%" = quantile(na.omit(CAC), c(.975)),
                   "Q1" = quantile(na.omit(CAC), c(.25)),
                   "Q3" = quantile(na.omit(CAC), c(.75))) %>% 
  as.data.frame() %>% na.omit()


res_low_cac_Intv$diff = res_low_cac_Intv$mediaani-res_low_cac_Intv$targetCAC
res_low_cac_Intv$rel_bias = paste0( round(((res_low_cac_Intv$mediaani-res_low_cac_Intv$targetCAC)/res_low_cac_Intv$targetCAC) *100), " %")

res_low_cac_Intv$Q1_rel_bias = paste0( round(((res_low_cac_Intv$Q1-res_low_cac_Intv$targetCAC)/res_low_cac_Intv$targetCAC) *100), " %")
res_low_cac_Intv$Q3_rel_bias = paste0( round(((res_low_cac_Intv$Q3-res_low_cac_Intv$targetCAC)/res_low_cac_Intv$targetCAC) *100), " %")

res_low_cac_Intv$n_iter = 1000
res_low_cac_Intv$conv_perc = paste0(round((1-res_low_cac_Intv$n_na/res_low_cac_Intv$n_iter)*100,1), " %")
res_low_cac_Intv$m = c(50,5,50,50,50,
                       50,5,50,50,50,
                       50,5,50,50,50)
res_low_cac_Intv$nClus = c(100,100,100,100,100,
                           50,50,50,50,50,
                           10,10,10,10,10)
res_low_cac_Intv$prev = c(0.5,0.5,0.9,0.2,0.02,
                          0.5,0.5,0.9,0.2,0.02,
                          0.5,0.5,0.9,0.2,0.02)

res_low_cac_Intv$scenario = c(1,2,3,4,5,
                              1,2,3,4,5,
                              1,2,3,4,5)
c = res_low_cac_Intv[,c("targetCAC","nClus","m","prev","mediaani","Q1","Q3",
                        "Q1_rel_bias","Q3_rel_bias",
                        "diff","rel_bias","conv_perc","scenario")]


c$"mediaani (Q1, Q3)" = paste0(round(c$mediaani,3)," (",round(c$Q1,3),", ",round(c$Q3,3),")" )
c$"Bias (Q1, Q3)" = paste0(c$rel_bias," (",c$Q1_rel_bias,", ",c$Q3_rel_bias,")" )

###########################
###########################
###########################
###########################


res_hi_cac_Intv = hi_cac_Intv %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(mediaani = median(CAC,na.rm=T), 
                   targetCAC = median(targetCAC,na.rm=T), 
                   n_na = sum(is.na(CAC)),
                   "2.5%" = quantile(na.omit(CAC), c(.025)),
                   "97.5%" = quantile(na.omit(CAC), c(.975)),
                   "Q1" = quantile(na.omit(CAC), c(.25)),
                   "Q3" = quantile(na.omit(CAC), c(.75))) %>% 
  as.data.frame() %>% na.omit()


res_hi_cac_Intv$diff = res_hi_cac_Intv$mediaani-res_hi_cac_Intv$targetCAC
res_hi_cac_Intv$rel_bias = paste0( round(((res_hi_cac_Intv$mediaani-res_hi_cac_Intv$targetCAC)/res_hi_cac_Intv$targetCAC) *100), " %")

res_hi_cac_Intv$Q1_rel_bias = paste0( round(((res_hi_cac_Intv$Q1-res_hi_cac_Intv$targetCAC)/res_hi_cac_Intv$targetCAC) *100), " %")
res_hi_cac_Intv$Q3_rel_bias = paste0( round(((res_hi_cac_Intv$Q3-res_hi_cac_Intv$targetCAC)/res_hi_cac_Intv$targetCAC) *100), " %")

res_hi_cac_Intv$n_iter = 1000
res_hi_cac_Intv$conv_perc = paste0(round((1-res_hi_cac_Intv$n_na/res_hi_cac_Intv$n_iter)*100,1), " %")
res_hi_cac_Intv$m = c(50,5,50,50,50,
                      50,5,50,50,50,
                      50,5,50,50,50)
res_hi_cac_Intv$nClus = c(100,100,100,100,100,
                          50,50,50,50,50,
                          10,10,10,10,10)
res_hi_cac_Intv$prev = c(0.5,0.5,0.9,0.2,0.02,
                         0.5,0.5,0.9,0.2,0.02,
                         0.5,0.5,0.9,0.2,0.02)
res_hi_cac_Intv$trtEff = 0.1
res_hi_cac_Intv$scenario = c(1,2,3,4,5,
                             1,2,3,4,5,
                             1,2,3,4,5)
d = res_hi_cac_Intv[,c("targetCAC","nClus","m","prev","mediaani","Q1","Q3",
                       "Q1_rel_bias","Q3_rel_bias",
                       "diff","rel_bias","conv_perc","scenario")]

d$"mediaani (Q1, Q3)" = paste0(round(d$mediaani,3)," (",round(d$Q1,3),", ",round(d$Q3,3),")" )
d$"Bias (Q1, Q3)" = paste0(d$rel_bias," (",d$Q1_rel_bias,", ",d$Q3_rel_bias,")" )

ab = rbind.data.frame(a,b)
cd = rbind.data.frame(c,d)


abcd = rbind.data.frame(ab,cd)

median(
  abs(as.numeric(gsub("%","",abcd$rel_bias))));

quantile(
  abs(as.numeric(gsub("%","",abcd$rel_bias))),probs=.25);

quantile(
  abs(as.numeric(gsub("%","",abcd$rel_bias))),probs=.75)

openxlsx::write.xlsx(ab,"CAC sim res, no intv 2026-07-10.xlsx")
openxlsx::write.xlsx(cd,"CAC sim res, intv 2026-07-10.xlsx")

median(as.numeric(gsub("%","",ab$rel_bias)));
quantile(as.numeric(gsub("%","",ab$rel_bias)),probs=.25)
quantile(as.numeric(gsub("%","",ab$rel_bias)),probs=.75)

median(as.numeric(gsub("%","",cd$rel_bias)));
quantile(as.numeric(gsub("%","",cd$rel_bias)),probs=.25)
quantile(as.numeric(gsub("%","",cd$rel_bias)),probs=.75)

