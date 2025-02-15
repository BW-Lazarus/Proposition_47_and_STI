// worker for estimation
args batch_num

sleep 5000

// load data
use processed_data.dta, clear

// load parameters
frame create "param"
frame change param
use tmp_param.dta, clear

local i = 0
foreach var of varlist _all {
	local `var' = `var'[1]
}

frame change default


// estimation
g treat = cond(statefip == `treated' & year >= 2015, 1, 0)

local start = (`batch_num' - 1)/`core_num'
di `start'
local step_size = 1/(`batch_count' * `core_num')
di `step_size'
local last = `start' + (`batch_count' - 1) * `step_size'
di `last'

// TODO doublecheck correctness
frame create "tmp_tn"
frame change tmp_tn
gen treated = .
gen alp_new = .
gen rmse = .
gen r2 = .
gen optimal_lambda = .
insobs `batch_count'
gen batch = `batch_num'
gen ord = .



frame change default

keep statefip year c_rate 
levelsof statefip, local(states)
greshape wide `oc', j(statefip) i(year)

// foreach state of local states {
// 	if `state' != `treated' {
// 		drop as_rate`state'
// 		drop homeless_rate`state'
// 	}
// }

tsset year
forvalues i = 1/`batch_count' {
	local pos = `i' - 1
	local alp_new = `start' + `pos' * `step_size'
	di "current alpha: `alp_new'"
// 	sculm `oc', ahead(1) q(0.16) treated(treat) lambda(lopt) cv(adaptive)
	
	cvlasso `oc'`treated' ///
	`oc'1-`oc'56 ///
	if year < 2015, ///
	lopt ///
	lglmnet ///
	roll ///
	h(1) adaptive postres alpha(`alp_new') ///
	prest maxiter(100000) tolopt(1e-6)  tolzero(1e-4) sklearn lminratio(1e-4) lmax(`lm') lcount(`lc') seed(`seed')
	
	predict double cf_`batch_num'_`pos', lopt
	
	frame change tmp_tn
	replace treated = `treated' if _n == `i'
	replace alp_new = `alp_new' if _n == `i'
	replace rmse = `e(rmse)' if _n == `i'
	replace r2 = `e(r2)' if _n == `i'
	replace optimal_lambda = `e(lambda)' if _n == `i'
	replace ord = `pos' if _n == `i'
	frame change default
}

if `batch_num' == `core_num' {
	di "current alpha: 1"
// 	sculm `oc', ahead(1) q(1) treated(treat) lambda(lopt) cv(adaptive)

	cvlasso `oc'`treated' ///
	`oc'1-`oc'56 ///
	if year < 2015, ///
	lopt ///
	lglmnet ///
	roll ///
	h(1) adaptive postres alpha(1) ///
	prest maxiter(200000) tolopt(1e-7)  tolzero(1e-4) sklearn lminratio(1e-4) lmax(`lm') lcount(`lc') seed(`seed')
	
	predict double cf_`batch_num'_`batch_count', lopt
	
	frame change tmp_tn
	insobs 1
	replace treated = `treated' if _n == _N
	replace alp_new = 1 if _n == _N
	replace rmse = `e(rmse)' if _n == _N
	replace r2 = `e(r2)' if _n == _N
	replace optimal_lambda = `e(lambda)' if _n == _N
	replace batch = `batch_num' if _n == _N
	replace ord = `batch_count' if _n == _N
	frame change default
}

local sa: label st `treated'
frame change tmp_tn
save `sa'_tn_`batch_num'.dta, replace

frame change default
save `sa'_predictions_`batch_num'.dta, replace

sleep 1000
// clear all 