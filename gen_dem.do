/* generate number of individuals of age group and gender cell

This .do file generates state level population estimates according to parameters
specified by the user, including study period (year), age brackets, and gender.
It uses the County-Level Population Files - Single-year Age Groups data from
the SEER website (https://seer.cancer.gov/popdata/download.html) 
*/
args first_year last_year age_lower_limit age_upper_limit gender
local age_range "`age_lower_limit'_`age_upper_limit'"
local time_range "`first_year'_`last_year'"
local pop_range "`age_range'_`time_range'"

// specify the paths of population estimates and working directory
use "some_path\pop_singleage_raw.dta", clear
cd "some_other_path"

// mac
// use "/Volumes/Extreme SSD/pop_comp/pop_singleage_raw.dta", clear
// cd "/Volumes/Extreme SSD/p47_std/std rates"

keep if year >= `first_year' & year <= `last_year'

keep if age >= `age_lower_limit' & age <= `age_upper_limit'

if `gender' > 0 {
	keep if sex == `gender'
}

preserve
keep if race == 2 & hispanic == 0

sort year statefips

by year statefips: egen nb_pop_tot_`pop_range' = total(pop)

by year statefips: gen id = _n

keep if id == 1
keep year statefips nb_pop_tot_`pop_range'

save nb_temp.dta, replace
restore

preserve
keep if race == 1 & hispanic == 0

sort year statefips

by year statefips: egen nw_pop_tot_`pop_range' = total(pop)

by year statefips: gen id = _n

keep if id == 1
keep year statefips nw_pop_tot_`pop_range'

save nw_temp.dta, replace
restore

preserve
keep if hispanic == 1

sort year statefips

by year statefips: egen hispanic_pop_tot_`pop_range' = total(pop)

by year statefips: gen id = _n

keep if id == 1
keep year statefips hispanic_pop_tot_`pop_range'

save hispanic_temp.dta, replace
restore

sort year statefips

by year statefips: egen pop_tot_`pop_range' = total(pop)

by year statefips: gen id = _n

keep if id == 1
keep year statefips pop_tot_`pop_range'

merge 1:1 statefips year using nb_temp.dta
keep if _merge == 3
drop _merge

merge 1:1 statefips year using nw_temp.dta
keep if _merge == 3
drop _merge

merge 1:1 statefips year using hispanic_temp.dta
keep if _merge == 3
drop _merge

merge 1:1 statefips year using total_pop.dta
keep if _merge == 3
drop _merge

gen pct_nw = nw_pop_tot_/pop_tot_
gen pct_nb = nb_pop_tot_/pop_tot_
gen pct_hispanic = hispanic_pop_tot_/pop_tot_
gen pct_gr = pop_tot_/pop_total

save pop_tot_`pop_range'_all_races, replace

















