


*** Table 1: Numbers

* Cluster period sizes

clear

import delimited "LAKANA Primary outcome data_chptr4.csv"

keep anon_id villageid death_* waz_who_* calendarmonth_* intv eligible_*

reshape long death_ waz_who_ eligible_ calendarmonth_ , i(anon_id) j(visit)

gen one = .
replace one = 1 if eligible_==1 

collapse (sum) one, by(villageid visit)
sum one, detail





clear
import delimited "LAKANA Growth outcomes data_de-ident_HIPAA_20251003.csv"

gen one = .
replace one = 1 if waz!="NA" | laz!="NA" | wlz!="NA"

collapse (sum) one, by(anon_villid visit)
sum one, detail




clear

import excel "LAKANA-AMR data_de-ident_HIPAA_2026-06-03.xlsx", firstrow sheet("S. pneumoniae")

gen one = .
replace one = 1 if Spneumoniaecultureresults!="."

collapse (sum) one, by(anon_village Visit)
sum one, detail





clear

import delimited "LAKANA EED mechanistic data, with EE score.csv"

keep aatresults_* anon_id villageid

reshape long aatresults_ , i(anon_id) j(visit)

gen one = .
replace one = 1 if aatresults_ !="."

collapse (sum) one, by(villageid visit)
sum one, detail