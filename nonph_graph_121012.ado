
********************************************************************************************************************************
********************************************************************************************************************************
* This .ado file implements the program `nonph_binary',
* which is used for generating graphs for variables that reject the proportional hazards assumption

* See the bottom of the file for examples of implementing the code


capture program drop nonph_graph
program define nonph_graphversion 11args graph ci lo hi* `j' defines the number of draws in each iteration
local j = 5000matrix beta = e(b)local b = colsof(beta)
matrix V = e(V)
scalar names = e(cmdline) * Define the lenght of time
quietly sum _t, d

scalar min = r(min)
scalar max = r(max)
scalar scale = max-minlocal counter_ci: word count `ci'local counter_min: word count `lo'local counter_max: word count `hi'/***********************************************//*** Combined Coefficient ***//***********************************************/if ("`graph'" == "combined_coefficient"){/*generate combined coef*/       gen combcoef = beta[1,1] + ln(_t)*beta[1,2]
       label variable combcoef "Combined Coefficient"/*calculate standard error of combined coef*/       gen se_combcoef = sqrt(V[1,1] + (ln(_t))^2*V[2,2] +2*ln(_t)*V[1,2])/*calculate confidence intervals*/if (`counter_ci' > 0){
local z = -invnormal((100-`ci')/200)       gen combcoef_hi = combcoef + `z'*se_combcoef       gen combcoef_lo = combcoef - `z'*se_combcoef
       label variable combcoef_hi "`ci'% Confidence Interval"       label variable combcoef_lo "`ci'% Confidence Interval"
} 
else {             gen combcoef_lo = combcoef - 1.96*se_combcoef       gen combcoef_hi = combcoef + 1.96*se_combcoef       
       label variable combcoef_hi "95% Confidence Interval"       label variable combcoef_lo "95% Confidence Interval"}       
       
/*chart effect with confidence bounds*/local name2 =word(names,2)       twoway line combcoef combcoef_lo combcoef_hi _t, sort yline(0) ///
       ytitle("Combined Coefficient") xtitle("Time") title("Combined Coefficient of `name2'") legend(order(1 2))
       
drop se_combcoef combcoef*
}/***********************************************//*** All Other Graphs ***//***********************************************/

*     ****************************************************************  **       Take `j' draws from the estimated coefficient vector and     **       variance-covariance matrix.                                     **     ****************************************************************  *preserve*local b = columnsquietly drawnorm b_1-b_`b', n(`j') means(e(b)) cov(e(V)) clear*     ****************************************************************  **       To calculate the desired quantities of interest we need to set  **       up a loop.  This is what we do here.                            **     ****************************************************************  **       First, specify what you quantities should be saved and what     **       these quantities should be called.                              **     ****************************************************************  *postutil clearpostfile mypost  diff_hat diff_lo diff_hi using sim , replace            *     ****************************************************************  **       Start loop.  Let this  											 **       run from min to max of time.   			   	                  **     ****************************************************************  *                                  /*Set the size of the matrix for stata */quietly set mat `j'/*Create a matrix of the simulated values from the betas and call it 'X'*/mkmat b_1-b_`b', matrix(X)/*Begin the loop*/local a=min while `a' <= max + scale/500 { /***********************************************//*** Hazard Ratio***//***********************************************/

if ("`graph'" == "hazard_ratio"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){
gen stuff = exp((`hi'-`lo')*(b_1 + b_2*ln(`a')))}}
else{
gen stuff = exp(b_1 + b_2*ln(`a'))
}}/***********************************************//*** First Difference ***//***********************************************/

if ("`graph'" == "first_difference"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){ gen stuff = (exp((`hi'-`lo')*(b_1 + b_2*ln(`a'))) - 1)*100
}
}
else{
 gen stuff = (exp(b_1 + b_2*ln(`a')) - 1)*100
}
}* Save the median and CI for each value of a        egen diff=median(stuff)        tempname  diff_hat diff_lo diff_hi   if (`counter_ci' > 0){
local x1 = (100-`ci')/2
local x2 = (100-(100-`ci')/2)
    _pctile stuff, p(`x1',`x2')		/*for a ci percent ci*/
}else{    _pctile stuff, p(2.5,97.5)		/*for a 95 percent ci*/}



    scalar `diff_lo'= r(r1)    scalar `diff_hi'= r(r2)      scalar `diff_hat'=diff * Post the values that have been generated           post mypost (`diff_hat') (`diff_lo') (`diff_hi')                   drop   diff  stuff 
* We are going to do the loop 500 times to get a pretty graph    local a=`a'+ (scale/500)     display "." _c    } postclose mypost* Use simulated values 

use sim, cleargen MV = ((_n-1)/500)*scale + min
quietly drop if _n > 501

*Generate the names of the variables to be used in the graphlocal name2 =word(names,2)* Generate the label for the graph

if (`counter_ci' > 0){
label variable diff_lo "`ci'% Confidence Interval"}else{
label variable diff_lo "95% Confidence Interval"}



if ("`graph'" == "hazard_ratio"){ 
if (`counter_min' > 0){
if (`counter_max' > 0){
label variable diff_hat "Hazard Ratio"

graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///        ||   line diff_lo  MV, clwidth(medium) clcolor(red) ///        ||   line diff_hi  MV, clwidth(medium) clcolor(red) ///        ||  ,   ///            yline(1) ///            legend(order(1 2)) ///            title("`name2' Increases From `lo' to `hi'") ///            xtitle("Time") ///            ytitle("Hazard Ratio") ///
			xlabel(minmax)
}}
else{
label variable diff_hat "Relative Hazard"

graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///        ||   line diff_lo  MV, clwidth(medium) clcolor(red) ///        ||   line diff_hi  MV, clwidth(medium) clcolor(red) ///        ||  ,   ///            yline(1) ///            legend(order(1 2)) ///            title("Relative Hazard of `name2'") ///            xtitle("Time") ///            ytitle("Relative Hazard") ///
			xlabel(minmax)
}}







if ("`graph'" == "first_difference"){ 
label variable diff_hat "First Difference"



graph twoway line diff_hat MV, clwidth(medium) clcolor(blue) ///        ||   line diff_lo  MV, clwidth(medium) clcolor(red) ///        ||   line diff_hi  MV, clwidth(medium) clcolor(red) ///        ||  ,   ///            yline(0) ///            legend(order(1 2)) ///            title("`name2'  Increases From `lo' to `hi'") ///            xtitle("Time") ///            ytitle("First Difference") ///
			xlabel(minmax)
}end



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
