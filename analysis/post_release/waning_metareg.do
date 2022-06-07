*cd "/Users/eh1415/Documents/covid-change-ve-over-time/release20220505"

*ssc install metareg

log using waning_metareg, replace

set more off
clear


import delimited data_metareg.csv

replace estimate="" if estimate=="NA"
destring estimate, replace

replace confhigh="" if confhigh=="NA"
destring confhigh, replace

replace conflow="" if conflow=="NA"
destring conflow, replace

gen loghr=estimate
gen lnu=confhigh
gen lnl=conflow

gen se1=(lnu-loghr)/1.96
gen se2=(loghr-lnl)/1.96

corr se1 se2

gen seloghr=(se1+se2)/2
// lines 26, 27, and 31 are usually simply coded as: gen seloghr = (lnu - lnl)/(2 * 1.96)
drop se1 se2 lnu lnl

* create strata from subgroup
gen stratum=1 if subgroup=="65+ years"
replace stratum=2 if subgroup=="18-64 years and clinically vulnerable"
replace stratum=3 if subgroup=="40-64 years"
replace stratum=4 if subgroup=="18-39 years"

* create sex
rename sex temp

gen sex=1 if temp=="Both"
replace sex=2 if temp=="Female"
replace sex=3 if temp=="Male"

label define sex 1 "Both" 2 "Female" 3 "Male" 
label values sex sex
drop temp

* create ageband
rename ageband temp
gen ageband=1 if temp=="all"
replace ageband=2 if temp=="65-74 years"
replace ageband=3 if temp=="75+ years"

drop temp

* create outcome
rename outcome temp

gen outcome=1 if temp=="covidadmitted"
replace outcome=2 if temp=="coviddeath"
replace outcome=3 if temp=="postest"
replace outcome=4 if temp=="noncoviddeath"
replace outcome=5 if temp=="anytest"

label define outcome 1 "COVID-19 hospitalisation" 2 "COVID-19 death" 3 "Positive test" ///
 4 "Non-COVID death" 5  "Any test"
label values outcome outcome
drop temp

* create vaccine from comparison
encode comparison, gen(vaccine)
tab vaccine
label list vaccine

replace k=k-1

sort outcome stratum sex vaccine
save waning_metareg.dta, replace

metareg loghr k if outcome==1 & stratum==1 & sex==1 & ageband==1 & vaccine==3, wsse(seloghr)
local a=_b[k]
local b=_se[k]
local c=_b[_cons]
local d=_se[_cons]

di `a'
di `b'
di `c'
di `d'

tempname memhold
postfile `memhold' outcome stratum sex ageband vaccine logrhr selogrhr loghr1 seloghr1 using results, replace

  forvalues i=1/5 {
  	forvalues v=1/3 {
		forvalues s=1/4 {
		forvalues x=1/3 {
			forvalues g=1/3 {
				di "A: " `i' `s' `g' `x' `v'

				count if outcome==`i' & stratum==`s' & sex==`g' & ageband==`x' & vaccine==`v' &loghr<.
				if r(N)>2 {
				di "B: " `i' `s' `g' `a' `v'
				metareg loghr k if outcome==`i' & stratum==`s' & sex==`g' & ageband==`x' & vaccine==`v', wsse(seloghr)
				local a=_b[k]
				local b=_se[k]
				local c=_b[_cons]
				local d=_se[_cons]
				post `memhold' (`i') (`s') (`g') (`x') (`v') (`a') (`b') (`c') (`d')
				}			
			}				
		}
		}
	}
  }

postclose `memhold'
di "`memhold'"

use results, clear
sort outcome stratum sex vaccine
save results, replace

log close
