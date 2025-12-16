/**************************************************************************
   EXTENSION: Heterogeneous Effects by Race and Gender
   Outcome: ln_cpi_income
*************************************************************************/

clear all
set more off

global root "/Users/lizbethmartinez/Desktop/Main Replications"

global tables    "$root/tables"
global figures   "$root/figures"

capture mkdir "$tables"
capture mkdir "$figures"


cd "$root"
use longrun_20002017acs_cleaned.dta, clear

*============================*
* 1. Interaction Terms
*============================*

* Race interaction
gen exp_black  = M12_exp_rate * black

* Gender interaction
gen exp_female = M12_exp_rate * female


*============================*
* 2. Combined Regression
*============================*

* Full specification with baseline controls
reg ln_cpi_income ///
    M12_exp_rate black exp_black ///
    female exp_female ///
    i.bpl i.birthyr i.year ///
    i.ageblackfemale i.bpl_black i.bpl_female i.bpl_black_female ///
    , cluster(bplcohort)

* Store results
eststo combined_ext


*============================*
* 3. Export Combined Table
*============================*

cd "$tables"

esttab combined_ext using "Extension_Heterogeneity.rtf", ///
    replace ///
    b(4) se(4) ///
    keep(M12_exp_rate exp_black exp_female black female) ///
    label ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    title("Extension: Heterogeneous Effects by Race and Gender") ///
    addnotes( ///
    "M12_exp_rate coefficient = effect for non-Black men." ///
    "exp_black coefficient = difference in effect for Black individuals." ///
    "exp_female coefficient = difference in effect for women." ///
    "Effect for Black individuals = M12_exp_rate + exp_black." ///
    "Effect for women = M12_exp_rate + exp_female." ///
    )

display "Extension Table Saved: $tables/Extension_Heterogeneity.rtf"

