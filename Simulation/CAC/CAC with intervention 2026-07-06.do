/*

Simulation: Comparing CAC from a model with identity link 
to one obtained via delta methd

Parameters 

- Number of clusters: n 
- Treatment effect: trtEff
- Target within-period ICC: target_icc_within
- Target between-period ICC: target_icc_between
- Prevalence/prob: p

- Cluster size: m

*/
clear


*log using "cac_twoGroup_comparison", replace text

capture program drop cac_intv_comparison_lo
program define cac_intv_comparison_lo,rclass

version 19.0
 syntax [, n(integer 1000) periods(integer 2) prev(real 0.9) target_cac(real 0.5) target_icc_within(real 0.1) trtEff(real .9)  m(integer 150)  ] // target_icc(real .01)
qui{
drop _all

/*
clear
local prev = 0.9
local trtEff = 0.8
local target_cac = 0.9
local target_icc_within = 0.9
local periods = 3
local n = 100
local m = 150

*/















** Get Between-period ICC from known CAC and Within-period ICC
local target_icc_between = `target_cac' * `target_icc_within'




** Total variance: 
* Within-period ICC = (sigma2_cluster + sigma2_periodCluster) / (sigma2_cluster + sigma2_periodCluster + sigma2_residual)
* Within-period ICC = (sigma2_cluster + sigma2_periodCluster) / Total variance
** Residual variance: sigma2_residual = p(1-p)

local p_bar =  `trtEff' * `prev' //`prev'
local residual_var = `p_bar'*(1-`p_bar')

** Rearranging from within-period ICC formula: sigma2_cluster + sigma2_periodCluster = wpICC*(sigma2_cluster + sigma2_periodCluster + sigma2_residual)
* Express A = sigma2_cluster + sigma2_periodCluster
* ==> A = wpICC*(A + sigma2_residual)
* ==> A = wpICC*A + wpICC*sigma2_residual
* ==> A - wpICC*A = wpICC*sigma2_residual
* ==> (1-wpICC)*A = wpICC*sigma2_residual
* ==> A = (wpICC*sigma2_residual) / (1-wpICC)
* ==> total variance = ( (wpICC*sigma2_residual) / (1-wpICC) ) + sigma2_residual
* ==> total variance = ( (wpICC*sigma2_residual) + (1-wpICC)*(sigma2_residual) / (1-wpICC) ) 
* ==> total variance =  sigma2_residual / (1-wpICC) // when wpICC<1 and wpICC>0
 
local total_var =  `residual_var' / (1 - `target_icc_within')  // total variance by target within period ICC



** Between period level parameters
local sig2_u = `target_icc_between' * `total_var'
local sig_u = sqrt(`sig2_u')


local a_clus = min((1-`prev')/2,`prev'/2)


 return scalar trtEffect = `trtEff'

 return scalar targetCAC = `target_cac'

 return scalar targetICCwithin = `target_icc_within'
 return scalar targetICCbetween = `target_icc_between'




** Within period level parameter
* bpICC = sigma2_cluster / total variance ==> sigma2_cluster = bpICC* total variance
* wpICC = (sigma2_cluster + sigma2_periodCluster) / Total variance  ==> (sigma2_cluster + sigma2_periodCluster) = wpICC* total variance
* ==> (bpICC* total variance + sigma2_periodCluster) = wpICC* total variance
* ==> sigma2_periodCluster = (wpICC-bpICC)* total variance

local sig2_v = (`target_icc_within' - `target_icc_between') * `total_var'
local sig_v = sqrt(`sig2_v')


* Set number of periods (time points)
set obs `periods'

gen period = _n
gen n = `n'


* Create data set
expand n

sort period

bysort period: gen clusID = _n
gen intv=0
replace intv=1 if clusID>round(`n'/2,1)


local low = normal((-`a_clus'-0)/`sig_u')
local high = normal((`a_clus'-0)/`sig_u')

gen help_u = runiform(`low', `high')
bys clusID: gen trunc_u = 0 + `sig_u' * invnormal(help_u)
bys clusID: replace trunc_u = trunc_u[1]

bys clusID: gen u = rnormal(0,`sig_u') if _n==1
bys clusID: replace u = u[1]


gen clusSize = round(max(2,rpoisson(`m')))

expand clusSize
gen indID = _n

bys clusID: gen i = _n

egen timegroup = group(clusID period)

local low_v = normal((-`a_clus'-0)/`sig_v')
local high_v = normal((`a_clus'-0)/`sig_v')
gen help_v = runiform(`low_v', `high_v')

bys clusID period: gen trunc_v = 0 + `sig_v' * invnormal(help_v)
bys clusID period: replace trunc_v = trunc_v[1]


bys clusID period: gen v = rnormal(0,`sig_v') if _n==1
bys clusID period: replace v = v[1]




gen p_i = `prev' + u + v  if intv ==0
replace p_i = `trtEff' * `prev' + u + v if intv == 1






gen problem = . 
replace problem = 1 if p_i>=0.999 | p_i<=0

count if problem == 1
local probability_problem = r(N)

if `probability_problem' > 0{
local flag = 1 
} 
else {
local flag = 0
}

return scalar probability_problem = `flag'

count
return scalar scale_of_probability_problem = `probability_problem'/r(N)

return scalar tailSize = `probability_problem'
return scalar allN = r(N)



qui summ p_i
return scalar orig_p_SD = r(sd)

replace p_i = 0.9999 if p_i>0.9999
replace p_i = 0.0001 if p_i<0.0001
* * * * * * *
gen outcome = rbinomial(1, p_i)



mixed outcome i.intv i.period || clusID: || timegroup:, reml

estat sd, post coeflegend
	
nlcom ( ( _b[clusID:sd(_cons)] )^2 ) / ( ( ( _b[clusID:sd(_cons)] )^2 ) + ( ( _b[timegroup:sd(_cons)] )^2 ) + ( _b[Residual:sd(e)])^2 ) , iterate(100000) 
matrix mat = r(table)
return scalar betweenICC = mat[1,1]
return scalar betweenICC_se = mat[2,1]

nlcom ( ( _b[clusID:sd(_cons)] )^2 + ( _b[timegroup:sd(_cons)] )^2 ) / ( ( ( _b[clusID:sd(_cons)] )^2 ) + ( ( _b[timegroup:sd(_cons)] )^2 ) + ( _b[Residual:sd(e)])^2 ) , iterate(100000) 
matrix mat = r(table)
return scalar withinICC = mat[1,1]
return scalar withinICC_se = mat[2,1]

nlcom ( ( _b[clusID:sd(_cons)] )^2 ) / ( ( ( _b[clusID:sd(_cons)] )^2 ) + ( ( _b[timegroup:sd(_cons)] )^2 ) ) , iterate(100000) 
matrix mat = r(table)
return scalar CAC = mat[1,1]
return scalar CAC_se = mat[2,1]



* Estimate p 
count if outcome==1
local numerator = r(N)
count if outcome==1 | outcome==0
local denominator = r(N)

local EmpProb = `numerator'/`denominator'



return scalar empiricalPrevalence = `EmpProb'

return scalar givenPrevalence = `prev'

return scalar empiricalSD_bound = sqrt( `EmpProb' * ( 1- `EmpProb') )

qui summ p_i
return scalar p_SD_bound = r(sd)




eret clear
}
end




** Save

/*
******** Scenario 1
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(100) ///
							periods(4) prev(0.5) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 1
save "low_Row_1.dta",replace



******** Scenario 2
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(100) ///
							prev(0.5) periods(4) target_cac(0.3) /// 
							target_icc_within(0.1) m(5) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 2
save "low_Row_2.dta",replace



******** Scenario 3
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(100) ///
							periods(4) prev(0.9) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 3
save "low_Row_3.dta",replace



******** Scenario 4
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(100) ///
							periods(4) prev(0.2) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 4
save "low_Row_4.dta",replace



******** Scenario 5
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(100) ///
							periods(4) prev(0.02) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 5
save "low_Row_5.dta",replace



*************************************
*************************************
*************************************


******** Scenario 6
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(50) ///
							periods(4) prev(0.5) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 6
save "low_Row_6.dta",replace



******** Scenario 7
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(50) ///
							periods(4) prev(0.5) target_cac(0.3) /// 
							target_icc_within(0.1) m(5) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 7
save "low_Row_7.dta",replace

*/

******** Scenario 8
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(50) ///
							periods(4) prev(0.9) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 8
save "low_Row_8.dta",replace



******** Scenario 9
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(50) ///
							periods(4) prev(0.2) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 9
save "low_Row_9.dta",replace



******** Scenario 10
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(50) ///
							periods(4) prev(0.02) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 10
save "low_Row_10.dta",replace




*************************************
*************************************
*************************************


******** Scenario 11
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(10) ///
							periods(4) prev(0.5) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 11
save "low_Row_11.dta",replace



******** Scenario 12
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(10) ///
							periods(4) prev(0.5) target_cac(0.3) /// 
							target_icc_within(0.1) m(5) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 12
save "low_Row_12.dta",replace



******** Scenario 13
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(10) ///
							periods(4) prev(0.9) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 13
save "low_Row_13.dta",replace



******** Scenario 14
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(10) ///
							periods(4) prev(0.2) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 14
save "low_Row_14.dta",replace



******** Scenario 15
simulate, seed(123) reps(1000): cac_intv_comparison_lo, n(10) ///
							periods(4) prev(0.02) target_cac(0.3) /// 
							target_icc_within(0.1) m(50) 				
gen lb = CAC - 1.96* CAC_se
gen ub = CAC + 1.96* CAC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetCAC & targetCAC <= ub
gen scenario = 15
save "low_Row_15.dta",replace






clear

use "low_Row_1.dta"

append using "low_Row_2.dta" ///
"low_Row_3.dta" ///
"low_Row_4.dta" ///
"low_Row_5.dta" ///
"low_Row_6.dta" ///
"low_Row_7.dta" ///
"low_Row_8.dta" ///
"low_Row_9.dta" ///
"low_Row_10.dta" ///
"low_Row_11.dta" ///
"low_Row_12.dta" ///
"low_Row_13.dta" ///
"low_Row_14.dta" ///
"low_Row_15.dta" 

save "low_CAC results_w_intv.dta", replace
	
