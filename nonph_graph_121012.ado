
********************************************************************************************************************************
********************************************************************************************************************************
* This .ado file implements the program `nonph_binary',
* which is used for generating graphs for variables that reject the proportional hazards assumption

* See the bottom of the file for examples of implementing the code


capture program drop nonph_graph
program define nonph_graph
local j = 5000
matrix V = e(V)
scalar names = e(cmdline) 
quietly sum _t, d

scalar min = r(min)
scalar max = r(max)
scalar scale = max-min
       label variable combcoef "Combined Coefficient"
local z = -invnormal((100-`ci')/200)
       label variable combcoef_hi "`ci'% Confidence Interval"
} 
else {      
       label variable combcoef_hi "95% Confidence Interval"
       
/*chart effect with confidence bounds*/
       ytitle("Combined Coefficient") xtitle("Time") title("Combined Coefficient of `name2'") legend(order(1 2))
       
drop se_combcoef combcoef*
}

*     ****************************************************************  *

if ("`graph'" == "hazard_ratio"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){
gen stuff = exp((`hi'-`lo')*(b_1 + b_2*ln(`a')))
else{
gen stuff = exp(b_1 + b_2*ln(`a'))
}

if ("`graph'" == "first_difference"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){
}
}
else{
 gen stuff = (exp(b_1 + b_2*ln(`a')) - 1)*100
}
}
local x1 = (100-`ci')/2
local x2 = (100-(100-`ci')/2)
    _pctile stuff, p(`x1',`x2')		/*for a ci percent ci*/
}






use sim, clear
quietly drop if _n > 501



if (`counter_ci' > 0){
label variable diff_lo "`ci'% Confidence Interval"
label variable diff_lo "95% Confidence Interval"



if ("`graph'" == "hazard_ratio"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){
label variable diff_hat "Hazard Ratio"

graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///
			xlabel(minmax)
}
else{
label variable diff_hat "Relative Hazard"

graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///
			xlabel(minmax)
}







if ("`graph'" == "first_difference"){ 
label variable diff_hat "First Difference"



graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///
			xlabel(minmax)
}



/*
**************************************************************************************
* An example of how to use the program
webuse stan3, clear
tab surgery
rename age Age
rename surgery Surgery
stcox Age posttran Surgery year, nohr
estat phtest, rank detail 
gen t_surgery = ln(_t)*Surgery
gen t_age = ln(_t)*Age






* To produce the combined coefficient of the binary variable `surgery':
stcox surgery t_surgery age posttran year
nonph_graph combined_coefficient
* Or, with a 90% confidence interval 
stcox surgery t_surgery age posttran year
nonph_graph combined_coefficient 90

* To produce the relative hazard of the variable `surgery' with a 90\% confidence interval
stcox surgery t_surgery age posttran year
nonph_graph hazard_ratio 90

* To produce the hazard ratio as age changes from 40 to 43 with a 94% confidence interval
stcox age t_age surgery posttran year
nonph_graph hazard_ratio 94 40 43

* To graph the first differences as age changes from 40 to 44 with a 99\% confidence interval
stcox age t_age surgery posttran year
nonph_graph first_difference 99 40 44


 
*/