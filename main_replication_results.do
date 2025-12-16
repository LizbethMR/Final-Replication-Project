/**************************************************************************
   main_replication_results.do
   Replicates: Table 2, Figure 1, Figure 3
*************************************************************************/
clear all
set more off

*============================*
* 1. Directory paths  *
*============================*

global root "/Users/lizbethmartinez/Desktop/Main Replications"

global output    "$root/output"
global figures   "$root/figures"
global tables    "$root/tables"

capture mkdir "$output"
capture mkdir "$figures"
capture mkdir "$tables"



*============================*
* 2. Table 2  Replication    *
*============================*

cd "$root"
use longrun_20002017acs_cleaned.dta, clear


* List of outcomes
local outcomes cpi_incwage cpi_incwage_no0 ln_cpi_income poverty100 employed hrs_worked

* Controls
local controls i.bpl i.birthyr i.ageblackfemale i.bpl_black i.bpl_female i.bpl_black_female black female

eststo clear

local i = 1
foreach y of local outcomes {
    qui reg `y' M12_exp_rate `controls' i.year, robust cluster(bplcohort)
    eststo model`i'
    local ++i
}

* COMPUTE OUTCOME MEANS (prevaccine cohorts)

matrix means = J(1,6,.)

local j = 1
foreach y of local outcomes {
    quietly summarize `y' if birthyr == 1947 | birthyr == 1948
    matrix means[1,`j'] = r(mean)
    local ++j
}

* AVERAGE 12-YEAR PREVACCINE INCIDENCE RATE
quietly summarize avg_12yr_measles_rate if birthyr == 1947 | birthyr == 1948
scalar mean_incidence_raw = r(mean)
scalar mean_incidence = mean_incidence_raw / 100000   // convert cases/100k -> population share

* IMPACT AT FULL EXPOSURE = β × 16 × incidence rate

matrix impact = J(1,6,.)
local k = 1

estimates dir
local models `r(names)'
foreach m of local models {
    estimates restore `m'
    scalar b = _b[M12_exp_rate]
    matrix impact[1,`k'] = b * 16 * mean_incidence
    local ++k
}

* EXPORT TABLE 

cd "$tables"

esttab model* using "Table2_clean.rtf", ///
    replace ///
    keep(M12_exp_rate) ///
    b(4) se(4) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    mtitle("Income" "Income>0" "ln Income" "Poverty" "Employed" "Hours Worked") ///
    scalars("r2 R-squared" "N Observations") ///
    nonumber noobs label ///
    title("Table 2 — Effects on Adult Labor Market Outcomes (Replication)") ///
    addnotes( ///
        "Outcome means (prevaccine cohorts): `=means[1,1]' `=means[1,2]' `=means[1,3]' `=means[1,4]' `=means[1,5]' `=means[1,6]'" ///
        "Average 12-year pre-vaccine incidence rate: `=mean_incidence'" ///
        "Impact at full exposure (16 years): `=impact[1,1]'  `=impact[1,2]'  `=impact[1,3]'  `=impact[1,4]'  `=impact[1,5]'  `=impact[1,6]'" ///
    )
	
* --- TXT version ---
esttab model* using "Table2_clean.txt", ///
    replace ///
    keep(M12_exp_rate) ///
    b(4) se(4)


display "Table 2 saved as: rtf, txt"

*============================*
* 3. Figure 1 Replication    *
*  (National disease counts) *
*============================*

cd "$root"
use case_counts_population.dta, clear

collapse (sum) measles (sum) pertussis (sum) chicken_pox ///
		 (sum) mumps (sum) rubella, by(year)

gen infectious_disease = pertussis + mumps + rubella + chicken_pox

keep if year >= 1956

twoway ///
    (line measles year,lcolor(blue)) ///
    (line infectious_disease year, lpattern(dash) lcolor(red)), ///
    xline(1964, lpattern(dash)) ///
    ytitle("Number of reported cases", margin(medlarge)) ///
    ylabel(, angle(horizontal)) ///
    xtitle("Year") ///
    legend(order(1 "Measles" 2 "Other infectious diseases")) ///
	title("Figure 1. Measles and Other Disease Incidence Over Time")

cd "$figures"
graph export "Figure1_replication.png", replace

*============================*
* 3. Figure 3 Replication    *
*  (Event Study on Measles)  *
*============================

cd "$root"
use inc_rate_ES.dta, clear

* interaction term for year and Measles pre-vac level in state
local i = 1
while `i' <= 5 {
    gen exp_Mpre_`i' = _Texp_`i' * avg_12yr_measles_rate
    local i = `i' + 1
}

local i = 7
while `i' <= 18 {
    gen exp_Mpre_`i' = _Texp_`i' * avg_12yr_measles_rate
    local i = `i' + 1
}

* run event-study regression
reg Measles exp_M* _Is* population _T* avg_12yr_measles_rate, cluster(statefip) robust

regsave, ci pval

* drop the unneeded coefficients and add in the 0 for omitted year
drop in 18/87

set obs 18
replace var      = "exp_Mpre_6" in 18
replace coef     = 0           in 18
replace stderr   = 0           in 18
replace N        = 1108        in 18
replace ci_lower = 0           in 18
replace ci_upper = 0           in 18

* event-time axis 
gen exp = 0
replace exp = -6 if var == "exp_Mpre_1"
replace exp = -5 if var == "exp_Mpre_2"
replace exp = -4 if var == "exp_Mpre_3"
replace exp = -3 if var == "exp_Mpre_4"
replace exp = -2 if var == "exp_Mpre_5"
replace exp = -1 if var == "exp_Mpre_6"
replace exp =  0 if var == "exp_Mpre_7"
replace exp =  1 if var == "exp_Mpre_8"
replace exp =  2 if var == "exp_Mpre_9"
replace exp =  3 if var == "exp_Mpre_10"
replace exp =  4 if var == "exp_Mpre_11"
replace exp =  5 if var == "exp_Mpre_12"
replace exp =  6 if var == "exp_Mpre_13"
replace exp =  7 if var == "exp_Mpre_14"
replace exp =  8 if var == "exp_Mpre_15"
replace exp =  9 if var == "exp_Mpre_16"
replace exp = 10 if var == "exp_Mpre_17"
replace exp = 11 if var == "exp_Mpre_18"

sort exp

* original-style event-study plot
scatter coef ci* exp if exp>-6 & exp<11, ///
    c(l l l) cmissing(y n n) ///
    msym(i i i) lcolor(gray gray gray) ///
    lpatter(solid dash dash) ///
    lwidth(thick medthick medthick) ///
    yline(0, lcolor(black)) xline(-1, lcolor(black)) ///
    ytitle("Measles rate by year (per 100,000)") ///
    xtitle("Years relative to measles vaccine availability") ///
    xlabel(-5(5)10, labsize(small)) ///
    ylabel(, nogrid angle(horizontal) labsize(small)) ///
    legend(off) ///
    title("Figure 3. Event Study – Effect of Measles Vaccine on Measles Incidence") ///
    graphregion(color(white))
	
cd "$figures"
graph export "Figure3_replication.png", replace

display "Figure 3"





