/*runner for regularized synthetic control

This .do file automates model selection and other commonly applied placebo/robustness
tests in the synthetic control literature for the elastic net synthetic control 
estimator (Doudchenko and Imbens, 2017) with cross-validation in both unit and 
time (Hollingsworth and Wing, 2020).

The model selection procedure is slightly different from Doudchenko and Imbens (2017) in
that each donor pool unit is associated with its own selected model. The implementation
is similar to Jared Greathouse's scul command, but differs in 
(1) systematic way of choosing both alpha (relative strength of L_1 and L_2 penalty)
and beta (overall penalty intensity) for the elastic net regression and
(2) computational efficiency

In addition to minimizing mspe, two conditions are applied when choosing models:
(1) models with optimal lambda at the boundaries are dropped and
(2) A candidate model M needs to have a few "neighboring  models" and the performance of them
could not drop significantly compared to M's.

This .do file incorporates Python code. Therefore, the user needs to set up PyStata. To increase
computational efficiency, the worker .do file uses the package cvlasso with sklearn option.
However, this could lead the Python backend to report a TypeError exception. The cause is that
constructors of sklearn's ElasticNet() and Ridge() class expect a Boolean in fit_intercept=sk_cons 
(line 4557, 4569, 4689 of lassoutils.ado in the lassopack package), but STATA does not have a 
Boolean data type. Thus, the user needs to cast sk_cons to Boolean type in Python by adding the 
following line of code after line 4553 and 4663 in lassoutils.ado:

stata("python: sk_cons = bool(sk_cons)")

This .do file uses the batcher package by Jesse Wursten that runs several STATA sessions in paralell,
which only runs on windows. 

arguments:
	outcome := outcome name; used for selecting the correct .dta file for analysis
	age_bracket := age bracket of the study population
	gender := gender of the study population
	cause := underlying cause of current std status
	time_bracket := study period
	alp := default alpha
	rmspe_threshold := units with rmspe in the pre-treatment period larger than 
	rmspe_threshold * rmspe of the treated unit will be discarded
	treated := state fips code of the treated state
	batch_count := number of alphas computed by each STATA session
	core_num := number of parallel STATA sessions
	lambda_exl := number of models with max/min lambdas deleted sequentially
	radius := number of valid alphas = radius * 2
	rf := the threshold "large" drop in performance is range(rmspe)/rf
	lm := max lambda
	lc := lm/lc is the smallest lambda
	seed := ensures replicability
	
example:
	cd "some_path"
	log using hiv_sc, replace
	di "execution starts"
	log off
	do std_runner_ml.do "chl" "55-64" "male" "na" "2000-2019" 1 10 6 20 5 3 3 2 20 400 2345
*/
args outcome age_bracket gender cause time_bracket alp rmspe_threshold treated batch_count core_num lambda_exl radius rf lm lc seed

clear all

if substr("`age_bracket'", 1, 1) == "0" {
	local ll 0
	local ul = substr("`age_bracket'", 3, 2)
}
else {
	local ll = substr("`age_bracket'", 1, 2)
	local ul = substr("`age_bracket'", 4, 2)
}

local fy = substr("`time_bracket'", 1, 4)
local ly = substr("`time_bracket'", -4, 4)

local name "`outcome'"


if "`gender'" != "both" {
	local name "`name'_`gender'"
}

if `ll' >= 0 & `ul' > 0 {
	local name "`name'_`ll'_`ul'"
}



if "`cause'" != "na" {
	local name "`name'_`cause'"
}


// di "`name'"
use `name'.dta, clear

save std_temp.dta, replace

if "`gender'" == "female" {
	local sex 2
}
if "`gender'" == "male" {
	local sex 1
}
if "`gender'" == "both" {
	local sex 0
}

do gen_dem_all_races.do `fy' `ly' `ll' `ul' `sex'
save dem_temp.dta, replace
use dem_temp.dta, clear
merge 1:1 statefips year using std_temp.dta
keep if _merge == 3
drop _merge

merge m:1 statefips using "E:\pop_comp\state_hook.dta"
keep if _merge == 3
drop _merge

gen c_rate = cases/pop_tot_ * 100000
gen lr = log(c_rate)
gen lp = log(pop_tot_)

levelsof statefips if c_rate == ., local(missing)

foreach state of local missing {
	drop if statefips == `state'
}

foreach num of numlist 2 49 9 40 {
	drop if statefips == `num'
}

log on
di "starts processing `name'"
log off

local oc "c_rate"

xtset statefip year
save processed_data.dta, replace
local sa: label st `treated'
local fn_main "`sa'_`name'_r"
local fn_alpha "`sa'_`name'_alp"
local fn_weights "`sa'_`name'_w"
local fn_pweights "`sa'_`name'_pw"
local fn_inference "`sa'_`name'_in"
local fn_placebo_ip "`sa'_`name'_ip"
local fn_placebo_it "`sa'_`name'_it"
local fn_robustness_loo "`sa'_name_rl"
local fn_param "`name'_param"

// parameters
frame create "`fn_param'"
frame change `fn_param'
gen treated = .
gen oc = "" 
gen alp = .
gen batch_count = .
gen core_num = .
gen lambda_exl = .
gen radius = .
gen rf = .
gen lm = .
gen lc = .
gen seed = .
insobs 1
// replace treated = `treated' if _n == 1
replace alp = `alp' if _n == 1
replace oc = "`oc'" if _n == 1
replace batch_count = `batch_count' if _n == 1
replace core_num = `core_num' if _n == 1
replace lambda_exl = `lambda_exl' if _n == 1
replace radius = `radius' if _n == 1
replace rf = `rf' if _n == 1
replace lm = `lm' if _n == 1
replace lc = `lc' if _n == 1
replace seed = `seed' if _n == 1
save tmp_param.dta, replace
frame change default

levelsof statefip, local(states)

// select model and store results
frame create "`fn_alpha'"
frame change `fn_alpha'
gen statefip = .
gen alpopt = .
// gen found = .
insobs 1

// 	sfi.Scalar.setValue("found", found)
// 	sfi.Scalar.setValue("alp_opt", alp_opt)

python:
def get_alp() -> float:
	import pandas as pd
	import numpy as np
	from functools import partial
	import warnings
    
	warnings.filterwarnings("ignore")
	
	def check_neighborhood(rmse_range: float, models: pd.DataFrame, params: dict, current_model: int) -> float:
		ul = int(params['radius']) + 1
		for idx in range(0, ul):
			cond_1 = (abs(models['rmse'][current_model] - models['rmse'][current_model + idx]) > rmse_range/params['rf'])
			cond_2 = (abs(models['rmse'][current_model] - models['rmse'][current_model - idx]) > rmse_range/params['rf'])
			cond_3 = (models['optimal_lambda'][current_model + idx] == -1)
			cond_4 = (models['optimal_lambda'][current_model - idx] == -1)
			if cond_1 or cond_2 or cond_3 or cond_4:
				return False
		return True


	df_models = pd.read_stata(r'E:\p47_std\std rates\tmp_test.dta')
	df_params = pd.read_stata(r'E:\p47_std\std rates\tmp_param.dta')

	params = df_params.to_dict('dict')
	params = {k:params[k][0] for k in params.keys()}

	rmse_range = df_models[df_models['tested'] == 0]['rmse'].max() - df_models[df_models['tested'] == 0]['rmse'].min()
	check = partial(check_neighborhood, rmse_range, df_models, params)
	found = False
	alp_opt = 0.5

	if len(df_models[df_models['tested'] == 0]) > 1 + 2 * params['radius']:
		while (len(df_models[df_models['tested'] == 0]) >= 1 + 2 * params['radius']) & (not found):
			current_model = df_models[df_models['tested'] == 0]['rmse'].idxmin()
			try:
				if check(current_model):
					found = True
					alp_opt = df_models['alp_new'].iloc[current_model]
					break
				else:
					df_models['tested'].iloc[current_model] = 1
			except KeyError:
				print('at extreme values!')
				break

		if not found:
			current_model = df_models[df_models['tested'] == 0]['rmse'].idxmin()
			alp_opt = df_models['alp_new'].iloc[current_model] 
		
	return alp_opt
end

// the loop
local idx = 1
foreach state of local states {
	frame change default
	local sa: label st `state'
	frame change `fn_param'
	replace treated = `state'
	save tmp_param.dta, replace
	
	log on
	di "`sa' written as parameters"
	log off
	
	mkdir "E:\p47_std\std rates\log"
	local temp_loc "E:\p47_std\std rates\log"

	batcher "worker_alpha.do", i(1/`core_num') tempfolder("`temp_loc'") maxparallel(`core_num') tr(1) be(3) up(10)

	sleep 2000
	!rmdir "E:\p47_std\std rates\log"  /s /q
	
	log on
	di "cross validation of alphas and model predictions for `sa' completed"
	log off
	
	// merge predictions
	frame create "tmp_preds"
	frame change tmp_preds
	use `sa'_predictions_1.dta, clear
	keep year `oc'`state' cf_*
	frame create "tmp_pred"

	forvalues i = 1/`core_num' {
		frame change default
		local sa: label st `state'
		frame change tmp_pred
		use `sa'_predictions_`i'.dta, clear
		keep year `oc'`state' cf_*
		save tmp_pred.dta, replace
		frame change tmp_preds
		merge 1:1 year `oc'`state' using tmp_pred.dta
		keep if _merge == 3
		drop _merge
	}
	save `sa'_preds.dta, replace
	frame drop tmp_pred
	frame change default
	frame drop tmp_preds
	
	log on
	di "predictions for `sa' recorded"
	log off
	
	
	// merge alpha-rmse data
	frame create "tmp_alphas"
	frame create "tmp_alpha"

	forvalues i = 1/`core_num' {
		frame change default
		local sa: label st `state'
		frame change tmp_alpha
		use `sa'_tn_`i'.dta, clear
		frame change tmp_alphas
		frameappend tmp_alpha
	}
	frame drop tmp_alpha
	
	log on
	di "`sa''s alphas recorded"
	log off
	sort alp_new

// 	local lambda_exl = 2
	frame change tmp_alphas
	sort alp_new
	gen id = _n
	save `sa'_models_raw.dta, replace
	preserve
	local lambda_exl 2
	local lm 20
	forvalues i = 1/`lambda_exl' {
		local lm 20
		local lmin = `lm' * 0.001
		qui sum optimal_lambda
		if _N != 0 {
			drop if (optimal_lambda == `r(min)') | (optimal_lambda == `lmin')
		}
		qui sum optimal_lambda
		if _N != 0 {
			drop if (optimal_lambda == `r(max)') | (optimal_lambda == `lm')
		}
		qui sum optimal_lambda
		if _N == 0 {
			break
		}
	}
	
	qui sum optimal_lambda
	if _N == 0 {
		log on
		di "`sa' has no valid alpha"
		log off
		local dn = 0
	}
	else {
		local dn = 1
	}
	restore
	
	if  `dn' == 0 {
		frame change `fn_alpha'
		replace statefip = `state' if _n == `idx'
		replace alpopt = -1 if _n == `idx'
		
		insobs 1
		local ++idx
		frame change default
		frame drop tmp_alphas
		
		log on
		di "iteration updated"
		log off
		
		continue
	}
	
	log on
	di "`sa' precheck for lambdas completed"
	log off
	
	qui sum rmse
	local r_range = `r(max)' - `r(min)'
	
	log on
	di "`sa' rmse range caculated"
	log off

	// exclude non-convergent/singular models
	noisily di `lambda_exl'
	noisily di `lm'
	forvalues i = 1/`lambda_exl' {
		local lmin = `lm' * 0.001
		qui sum optimal_lambda
		replace optimal_lambda = -1 if optimal_lambda == r(min) | (optimal_lambda == `lmin')
		
		qui sum optimal_lambda
		replace optimal_lambda = -1 if (optimal_lambda == r(max)) | (optimal_lambda == `lm')
	}
	
	log on
	di "modified alphas for `sa'"
	log off
	
	// model uncertainty
	gen tested = 0
	replace tested = 1 if optimal_lambda == -1
	save tmp_test.dta, replace
	frame change default
	frame drop tmp_alphas
	
	log on
	di "`sa' alphas saved"
	log off
	
	python: import sfi
	python: sfi.Scalar.setValue("alp_opt", get_alp())
	
	log on
	di "`sa' optimal alpha picked"
	log off
	
	frame change `fn_alpha'
	replace statefip = `state' if _n == `idx'
	replace alpopt = alp_opt if _n == `idx'
	
	log on
	di "optimal alpha for `sa' recorded"
	log off
	
	insobs 1
	local ++idx
	frame change default
	
	log on
	di "iteration updated"
	log off
}

log on
di "model selection finished"
log off

frame change `fn_alpha'
drop if statefip == .
save optimal_alpha.dta, replace

log on
di "optimal alphas recorded"
log off


frame change default
local sa: label st `treated'
frame change `fn_alpha'
levelsof alpopt if statefip == `treated'
local this_alp = `r(levels)'

frame create "tmp_ml"
frame change tmp_ml
use `sa'_models_raw.dta, clear

gen alp_diff = alp_new - `this_alp'
replace alp_diff = alp_diff * alp_diff
qui sum alp_diff
replace alp_diff = alp_diff - 2 * `r(min)'
preserve
keep if alp_diff <= 0
local this_model = id[1]
restore

levelsof batch if id == `this_model'
local this_batch = `r(levels)'

levelsof ord if id == `this_model'
local this_ord = `r(levels)'

use `sa'_preds.dta, clear
gen tr = `oc'`treated' - cf_`this_batch'_`this_ord'

log on
di "level comparison and treatment effects preparation done"
log off

graph twoway (line `oc'`treated' year, lcolor(black) lwidth(medium) lpattern(solid)) (line cf_`this_batch'_`this_ord' year, lcolor(black) lwidth(medium) lpattern(dash)), legend(order(1 "`sa'" 2 "Synthetic `sa'") size(small) position(6) rows(1)) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2014, lcolor(gs10%50)) yline(0, lcolor(gs10%50)) title("`sa', Actual vs. Synthetic", size(medium)) 
// text(60 2018.5 "ATT = `att'", place(a))


sleep 5000
graph export std_graph\pred_`name'.pdf, as(pdf) replace
log on
di "level comparison graph finished."
log off

//
// treatment effect
// local treated 6
graph twoway (line tr year, lcolor(black) lwidth(medium) lpattern(solid)), legend(off) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2014, lcolor(gs10%50)) yline(0, lcolor(gs10%50)) title("Treatment Effects", size(medium))

sleep 5000
graph export std_graph\te_`name'.pdf, as(pdf) replace
frame change default
frame drop tmp_ml

log on
di "treatment effect graph finished."
log off

// merge predictions
frame change default

levelsof statefip, local(states)

frame create "`fn_placebo_ip'"
frame create "tmp_models"
frame create "tmp_preds"

frame change `fn_alpha'
gen ladopt = .
gen rmse_pre = .
gen rmse_post = .
gen rmse_ratio = .
gen att = .
gen adj_att = .

log on
di "model performance updated"
log off

foreach state of local states {
	frame change default
	local sa: label st `state'
	
	frame change `fn_alpha'
	levelsof alpopt if statefip == `state'
	local this_alp = `r(levels)'
	
	log on
	di "`sa' alpha got"
	log off
	
	if `this_alp' == -1 {
		frame change default
		log on
		di "`sa' has invalid alpha; better use other models"
		log off
		continue
	}

	frame change tmp_models
	use `sa'_models_raw.dta, clear
	
	gen alp_diff = alp_new - `this_alp'
	replace alp_diff = alp_diff * alp_diff
	qui sum alp_diff
	replace alp_diff = alp_diff - 2 * `r(min)'
	preserve
	keep if alp_diff <= 0
	local this_model = id[1]
	restore
	
	drop alp_diff
	levelsof batch if id == `this_model'
	local this_batch = `r(levels)'

	levelsof ord if id == `this_model'
	local this_ord = `r(levels)'
	
	levelsof optimal_lambda if id == `this_model'
	local this_lambda = `r(levels)'
	
	log on
	di "`sa' model found"
	log off
	
	frame change tmp_preds
	use `sa'_preds.dta, clear
	keep year `oc'`state' cf_`this_batch'_`this_ord'
	rename `oc'`state' `oc'
	rename cf_`this_batch'_`this_ord' pred
	gen diff = `oc' - pred
	gen time = year - 2015
	gen post = (year >= 2015)
	
	log on
	di "`sa''s difference calculated"
	log off
	
	// att
	sort post
	by post: egen total_effect = total(diff)
	count if post == 1
	gen att = total_effect/r(N)
	
	log on
	di "`sa''s att calculated"
	log off
	
	// rmse
	sort post
	gen diff2 = diff * diff
	by post: egen mse = total(diff2)
	count if post == 1
	replace mse = mse/r(N)
	gen rmse = sqrt(mse)
	
	log on
	di "`sa''s rmse calculated"
	log off
	
	sort year
	local rmse_pre = rmse[1]
	qui sum
	local rmse_post = rmse[r(N)]
	qui sum
	local att = att[r(N)]
	
	log on
	di "`sa''s rmse and att saved to locals"
	log off
	
	keep year `oc' pred diff time 
	gen statefip = `state'
	
	// fill model performance data
	frame change `fn_alpha'
	replace ladopt = `this_lambda' if statefip == `state'
	replace rmse_pre = `rmse_pre' if statefip == `state'
	replace rmse_post = `rmse_post' if statefip == `state'
	replace rmse_ratio = `rmse_post'/`rmse_pre' if statefip == `state'
	replace adj_att = `att'/`rmse_pre' if statefip == `state'
	replace att = `att' if statefip == `state'
	
	log on
	di "`sa''s model performance updated"
	log off
	
	// append predictions
	frame change `fn_placebo_ip'
	frameappend tmp_preds
	
	log on
	di "`sa''s predictions updated"
	log off
}

// save the aggreate results
frame change `fn_alpha'
save model_performance.dta, replace
log on
di "model performance recorded"
log off

frame change `fn_placebo_ip'
save model_predictions.dta, replace
log on
di "model predictions recorded"
log off

frame change default
frame drop tmp_models
frame drop tmp_preds

// plot in-place placebo test
frame change `fn_placebo_ip'
local draw_list 
local treated 6
local sa: label st `treated'

levelsof statefip, local(states)

foreach state of local states {
	if `state' == `treated' {
		continue
	}
	local draw_list "`draw_list' (line diff year if statefip == `state', lcolor(gs10%60) lwidth(vthin))"
}

local draw_list "`draw_list' (line diff year if statefip == `treated', lcolor(black) lwidth(medium))"


graph twoway `draw_list', legend(off) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Treatment Effects", size(small)) xtitle("Time", size(small)) xline(2014, lcolor(gs10%50)) yline(0, lcolor(gs10%60)) title("In-Place Placebo Test", size(medium))

sleep 5000
graph export std_graph\placebo_ip_`name'.pdf, as(pdf) replace

log on
di "in-place placebo test graph finished."
log off

// in-time placebo test using synth2
// due to small number of pre-treatment periods
frame change default

synth2 `oc' `oc'(2008) `oc'(2009) `oc'(2010) `oc'(2011) `oc'(2012) `oc'(2013) `oc'(2014) , trunit(`treated')  trperiod(2015)  xperiod(2008(1)2014) maxiter(50000) sigf(8) nested allopt savegraph(`treated', replace) frame("`fn_placebo_it'") placebo(period(2011) show(10))

frame change `fn_placebo_it'

local att: di %4.3f `e(att)'
graph twoway (line `oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(solid)) (line pred·`oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(dash)), legend(order(1 "`sa'" 2 "Synthetic `sa'") size(small) position(6) rows(1)) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2015, lcolor(gs10%50)) yline(0, lcolor(gs10%60)) title("`sa', Actual vs. Synthetic", size(medium)) subtitle("Alternative Estimator", size(small)) text(60 2018.5 "ATT = `att'", place(a))

sleep 5000
graph export std_graph\placebo_am_`name'.pdf, as(pdf) replace

graph twoway (line `oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(solid)) (line pred·`oc'·2011 year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(dash)), legend(order(1 "`sa'" 2 "Synthetic `sa'") size(small) position(6) rows(1)) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2011, lcolor(gs10%50)) xline(2015, lcolor(gs10%50)) yline(0, lcolor(gs10%60)) title("In-Time Placebo Test, `sa'", size(medium)) subtitle("Alternative Estimator", size(small))

sleep 5000
graph export std_graph\placebo_it_`name'.pdf, as(pdf) replace

frame change default
// local fn_placebo_it  "CA_hivpre_male_55_64_idu_alp"
frame drop `fn_placebo_it'

frame change default
local treated 6
local oc c_rate

// g treat = cond(statefip == `treated' & year >= 2015, 1, 0)

// scul `oc', ahead(1) q(0.5) treated(treat) lambda(lopt)  obscol(black) cfcol(blue) legpos(7) cv(adaptive) transform(`oc')

local sa: label st `treated'

synth2 `oc' `oc'(2008) `oc'(2009) `oc'(2010) `oc'(2011) `oc'(2012) `oc'(2013) `oc'(2014) , trunit(`treated')  trperiod(2015)  xperiod(2008(1)2014) maxiter(50000) sigf(8) nested allopt savegraph(`treated', replace) frame("`fn_placebo_it'") placebo(period(2011) show(10))

frame change `fn_placebo_it'

local att: di %4.3f `e(att)'
graph twoway (line `oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(solid)) (line pred·`oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(dash)), legend(order(1 "`sa'" 2 "Synthetic `sa'") size(small) position(6) rows(1)) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2015, lcolor(gs10%50)) yline(0, lcolor(gs10%60)) title("`sa', Actual vs. Synthetic", size(medium)) subtitle("Alternative Estimator", size(small)) text(60 2018.5 "ATT = `att'", place(a))

sleep 5000
graph export std_graph\placebo_am_CA.pdf, as(pdf) replace

graph twoway (line `oc' year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(solid)) (line pred·`oc'·2011 year if statefip == `treated', lcolor(black) lwidth(medium) lpattern(dash)), legend(order(1 "`sa'" 2 "Synthetic `sa'") size(small) position(6) rows(1)) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2011, lcolor(gs10%50)) xline(2015, lcolor(gs10%50)) yline(0, lcolor(gs10%60)) title("In-Time Placebo Test, `sa'", size(medium)) subtitle("Alternative Estimator", size(small))

sleep 5000
graph export std_graph\placebo_it_CA.pdf, as(pdf) replace
frame change default
log on
di "ADH(2015) in-time placebo test graph finished."
log off

// plot rmspe ratio ranking
frame change default
local treated 6
local sa: label st `treated'
frame change `fn_alpha'
preserve
keep if alpopt >= 0
levelsof rmse_pre if statefip == `treated'
local rmsep_treated = `r(levels)'
keep if rmse_pre < `rmsep_treated' * `rmspe_threshold'
gen srmser = rmse_ratio * sign(att)
sort adj_att
gen rank = _N - _n + 1
keep if rank <= 10

graph hbar adj_att, over(statefip, sort(rank)) b1title("Adjusted ATT", size(small)) title("ATT/Pre RMSE Ranking, Top 10 States", size(medium)) bar(1, color(black)) ytitle("")

sleep 5000
graph export std_graph\placebo_rr1_`name'.pdf, as(pdf) replace
restore

local treated 6
local sa: label st `treated'
preserve
keep if alpopt >= 0
levelsof rmse_pre if statefip == `treated'
local rmsep_treated = `r(levels)'
keep if rmse_pre < `rmsep_treated' * `rmspe_threshold'
gen srmser = rmse_ratio * sign(att)
sort srmser
gen rank = _N - _n + 1
keep if rank <= 10

graph hbar srmser, over(statefip, sort(rank)) b1title("RMSPE Ratio * Sign(ATT)", size(small)) title("Post/Pre Treatment Signed RMSPE Ratio Ranking, Top 10 States", size(medium)) bar(1, color(black)) ytitle("")

sleep 5000
graph export std_graph\placebo_rr2_`name'.pdf, as(pdf) replace
restore

frame change default

log on
di "rmspe ranking graph finished."
log off
dd


// leave-one-out robustness checks
// get donor states with positive weights
frame change `fn_weights'
local pw
local prefix_len = 2 + strlen("`oc'") + 1

foreach var of varlist _all {
	if substr("`var'", -1, 1) == "s" {
		continue
	}
	else {
		local temp = substr("`var'", `prefix_len' , .)
		local pw "`pw' `temp'"
	}
}

frame change default


mkdir "E:\p47_std\std rates\log"
local temp_loc "E:\p47_std\std rates\log"

frame change `fn_param'
gen treated = .
replace treated = `treated' if _n == 1
save tmp_param.dta, replace
frame change default

batcher "worker_loo.do", i("`pw'") tempfolder("`temp_loc'") maxparallel(4) tr(5) be(10) up(30)
sleep 5000
!rmdir "E:\p47_std\std rates\log"  /s /q

// foreach d of local pw {
// 	local fn_loo_temp "`fn_main'_`d'"
// 	frame create "`fn_loo_temp'"
// 	frame change `fn_loo_temp'
// 	use scul_`sa'.dta, clear
// 	// label each iteration
// 	gen fid = `d'
// 	frame change default
// }

frame create `"fn_robustness_loo"'
frame create tmp_rb
frame change `fn_main'

frame copy `fn_main' `fn_robustness_loo', replace
frame change `fn_robustness_loo'
gen fid = `treated'


foreach d of local pw {
	frame change tmp_rb
	use scul_`sa'.dta, clear
	gen fid = `d'
	frame change `fn_robustness_loo'
	frameappend `fn_loo_temp'
}


// local oc "c_rate"
// local treated 6
local draw_list "(line `oc'`treated' year if fid == `treated', lcolor(black) lwidth(medium) lpattern(solid)) (line cf year if fid == `treated', lcolor(black) lwidth(medium) lpattern(dash))"

levelsof fid, local(states)

foreach state of local states {
	if `state' == `treated' {
		continue
	}
	local draw_list "`draw_list' (line cf year if fid == `state', lcolor(gs10%50) lwidth(vthin))"
}

graph twoway `draw_list', legend(off) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2014, lcolor(gs8%60)) yline(0, lcolor(gs10%50)) title("Leave-One-Out Robustness Test", size(medium))

sleep 5000
graph export std_graph\placebo_rb_`name'.pdf, as(pdf) replace
frame change default

log on
di "leave-one-out robustness test graph finished."
log off
// save frames to disk and clear memory
local fn_main "`sa'_`name'_r"
local fn_weights "`sa'_`name'_w"
local fn_pweights "`sa'_`name'_pw"
local fn_inference "`sa'_`name'_in"
local fn_placebo_ip "`sa'_`name'_ip"
local fn_placebo_it "`sa'_`name'_it"
local fn_robustness_loo "`sa'_name_rl"
local fn_param "`name'_param"
local frame_list `fn_main' `fn_weights' `fn_pweights' `fn_placebo_ip' `fn_placebo_it' `fn_param' `fn_robustness_loo'
foreach fn of local frame_list {
	frame change `fn'
	save `fn'.dta, replace
	frame change default
	frame drop `fn'
}

log on
di "`name' estimation finished."
log off

graph twoway line tr·c_rate·2011 year, lcolor(black) lwidth(medium) lpattern(solid)  legend(off) xlabel(2008(4)2020, labsize(small) nogrid) ylabel(, labsize(small) nogrid) ytitle("Number of Cases per 100,000", size(small)) xtitle("Time", size(small)) xline(2014, lcolor(gs10%50)) title("In-Time Placebo Test", size(medium)) xline(2011, lcolor(gs10%50)) yline(0, lcolor(gs10%50))

graph export graphs/placebo_it_`name'.png, replace

save processed_`name'.dta, replace

