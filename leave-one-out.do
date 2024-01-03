// further robustness test
// leave-one-out distribution of event study point estimates and CI bounds
use "E:\pop_comp\state_abv_and_fips.dta", clear
cd "C:\Users\TanBW\Desktop\crime graphs\"

// number of iteration
local iter_num = 500
gen index = _n

// for reproducibility
set seed 339487731

// generate random indicators for deletion of one state
foreach i of numlist 1/`iter_num'  {
	gen chosen_states`i' = 0
	local a = runiformint(1,51)
	replace chosen_states`i' = 1 if index == `a'
}

rename statefip statefips
save random_states.dta, replace

use "E:\pop_comp\final_with_controls_min.dta", clear
cd "C:\Users\TanBW\Desktop\crime graphs\"

merge 1:1 statefips year using "E:\pop_comp\nvss_hr_age_adjusted.dta"
replace ar_r = hr1 if ar_r == -1 | ar_r == -2
drop _merge

// the basic specification

// state and year dummies
tab year, gen(yy) nof
tab statefips, gen(ss) nof

// log state population
gen txxpop=ln(pop_tot)

// interaction terms between pre-reformulation OxyContin misuse rates and year
// dummies, excluding year 2009 (following Alpert et al. (2018))
foreach num of numlist 2000/2008 2010/2017 {
	gen ee`num'=(initial_oxy)*(year==`num')
}

// interaction terms between pre-reformulation other pain relievers misuse rates
// and year dummies
foreach num of numlist 2000/2017 {
	gen cc`num'=initial_pr*(year==`num')
}


// states bordering Mexico
gen border_m = 0
replace border_m = 1 if statefips == 6 | statefips == 4 | statefips == 35 | statefips == 48

// Florida
gen border_f = 0
replace border_f = 1 if statefips == 12

// Hawaii and Alaska
gen offshore = 0
replace offshore = 1 if statefips == 2 | statefips == 15

// states that harbor major heroin markets
gen heroin_seizure = 0
replace heroin_seizure = 1 if statefips == 4 | statefips == 8 | statefips == 12 | statefips == 17 | statefips == 36 | statefips == 53 | statefips == 34

// time trends for states considered as sources and distribution hubs for heroin
gen trend_m = (year - 2000)*border_m
gen trend_f = (year - 2000)*border_f
gen trend_o = (year - 2000)*offshore
gen trend_h = (year - 2000)*heroin_seizure


merge m:1 statefips using random_states.dta
drop _merge

// initialize point estimates and their bonds
foreach i of numlist 1/`iter_num' {
	gen beta`i' = 0
	gen se`i' = .
	gen high`i' = 0
	gen low`i' = 0
}

// estimation 
foreach i of numlist 1/`iter_num' {
	reg ar_r ss* yy* ee* cc* trend_m trend_h trend_f trend_o [aw=pop_tot] if chosen_states`i' != 1, cluster(statefips) 
	foreach num of numlist 2000/2008 2010/2017 {
		replace beta`i' = _b[ee`num'] if year==`num'
		replace se`i' = _se[ee`num'] if year == `num'
	}
	replace high`i' = beta`i' + 1.96 * se`i'
	replace low`i' = beta`i' - 1.96 * se`i'
}


// calculating averages of parameters
gen beta_sum = 0
gen high_sum = 0
gen low_sum = 0


foreach i of numlist 1/`iter_num' {
	replace beta_sum = beta_sum + beta`i'
	replace high_sum = high_sum + high`i'
	replace low_sum = low_sum + low`i'
}

gen beta_avg = beta_sum/`iter_num'
gen high_avg = high_sum/`iter_num'
gen low_avg = low_sum/`iter_num'

foreach i of numlist 1/`iter_num' {
	replace high`i' = (high`i' - high_avg) * (high`i' - high_avg)
	replace low`i' = (low`i' - low_avg) * (low`i' - low_avg)
}

replace high_sum = 0
replace low_sum = 0

foreach i of numlist 1/`iter_num' {
	replace high_sum = high_sum + high`i'
	replace low_sum = low_sum + low`i'
}

gen high_se = sqrt(high_sum/(`iter_num' - 1))
gen low_se = sqrt(low_sum/(`iter_num' - 1))

gen high_upper = high_avg + 1.96 * high_se
gen low_lower = low_avg - 1.96 * low_se


// graph
	# delimit ;
twoway (rarea high_avg low_avg year if year>=2000 & year<2018& year~=2009 , sort fcolor(dknavy%50) fintensity(30)	lcolor(white) lwidth(vthin)) || (connected beta_avg year if year>1999 & year<2018,  sort yaxis(1) msymbol(O) msize(tiny) lcolor(dknavy) mcolor(dknavy) lwidth(vthin))  , yline(0, lcolor(gs0)) graphregion(color(white)) bgcolor(white) xline(2010, lcolor(gs0)) legend(label(2 "point estimates") label(1 "95% CI")) xlabel(2000(1)2017, angle(45) labsize(vsmall) nogrid)  ylabel(-3(1.5)3, labsize(vsmall) nogrid)  xtitle("Year", size(small)) title("(a) Average Bounds", size(small)) ytitle("Average Estimated Coefficient" , axis(1) size(small));
# delimit cr
graph save graphs/homicide_leave_one_out_avg, replace
drop high* low*