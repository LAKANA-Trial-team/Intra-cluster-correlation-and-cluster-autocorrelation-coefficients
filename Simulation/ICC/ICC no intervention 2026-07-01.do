clear
capture program drop icc_oneGroup_comparison
program define icc_oneGroup_comparison,rclass
version 19.0
 syntax [, n(integer 1000) prev(real 0.9) target_icc(real 0.1)  m(integer 150)  ] 
qui{
drop _all
	return scalar targetICC = `target_icc'
	set obs  `n'
	gen clusID = _n
	gen clusSize = round(max(2,rpoisson(`m'))) 
// Between-cluster variance
	local sig2_u = ( `target_icc'* (`prev'*(1-`prev') ) )/(1-`target_icc') 
	local sig_u = sqrt(`sig2_u')
	// Instead rnormal, generate noise from beta
	local alpha = (`prev'*(1-`prev'))/(`sig2_u') - 1
	gen clusP = rbeta(`alpha'*`prev',`alpha'*(1-`prev')) // if intv==0
	expand clusSize
	* Generate outcome
	gen outcome = rbinomial(1, clusP)
	mixed outcome || clusID: , reml
	qui estat icc	
	return scalar identityICC = r(icc2)
	return scalar identityICC_se = r(se2)
eret clear
}
end


** Save


******** Scenario 1
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(100) ///
							prev(0.5) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 1
save "Row_1.dta",replace



******** Scenario 2
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(100) ///
							prev(0.5) target_icc(0.1) m(5) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 2
save "Row_2.dta",replace



******** Scenario 3
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(100) ///
							prev(0.9) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 3
save "Row_3.dta",replace



******** Scenario 4
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(100) ///
							prev(0.20) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 4
save "Row_4.dta",replace



******** Scenario 5
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(100) ///
							prev(0.02) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 5
save "Row_5.dta",replace



*************************************
*************************************
*************************************


******** Scenario 6
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(50) ///
							prev(0.5) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 6
save "Row_6.dta",replace



******** Scenario 7
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(50) ///
							prev(0.5) target_icc(0.1) m(5) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 7
save "Row_7.dta",replace



******** Scenario 8
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(50) ///
							prev(0.9) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 8
save "Row_8.dta",replace



******** Scenario 9
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(50) ///
							prev(0.20) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 9
save "Row_9.dta",replace



******** Scenario 10
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(50) ///
							prev(0.02) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 10
save "Row_10.dta",replace




*************************************
*************************************
*************************************


******** Scenario 11
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(10) ///
							prev(0.5) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 11
save "Row_11.dta",replace



******** Scenario 12
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(10) ///
							prev(0.5) target_icc(0.1) m(5) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 12
save "Row_12.dta",replace



******** Scenario 13
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(10) ///
							prev(0.9) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 13
save "Row_13.dta",replace



******** Scenario 14
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(10) ///
							prev(0.20) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 14
save "Row_14.dta",replace



******** Scenario 15
simulate, seed(123) reps(1000): icc_oneGroup_comparison, n(10) ///
							prev(0.02) target_icc(0.1) m(50) 				
gen lb = identityICC - 1.96* identityICC_se
gen ub = identityICC + 1.96* identityICC_se
gen withinCI = 0
replace withinCI = 1 if lb <= targetICC & targetICC <= ub
gen scenario = 15
save "Row_15.dta",replace











clear

use "Row_1.dta"

append using "Row_2.dta" ///
"Row_3.dta" ///
"Row_4.dta" ///
"Row_5.dta" ///
"Row_6.dta" ///
"Row_7.dta" ///
"Row_8.dta" ///
"Row_9.dta" ///
"Row_10.dta" ///
"Row_11.dta" ///
"Row_12.dta" ///
"Row_13.dta" ///
"Row_14.dta" ///
"Row_15.dta" 

save "onegroup ICC results_no intv.dta", replace