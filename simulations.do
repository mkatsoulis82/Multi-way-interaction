* This script is written by Manuel Gomes


*******************************************************************************
**** Simulation study - Multi - way interactions 

**** Katsoulis et al, 2020


clear all
set seed 110582

capture program drop MGsim
program MGsim, rclass
clear 

* Generate binary indicators with a certain prevalence and correlation (this is for the 2nd scenario)
* using the Stata command "rbinary" presented in the following paper
* Minxing Chen. Generating nonnegatively correlated binary random variates. The Stata Journal (2015) 15, Number 1, pp. 301–308

rbinary x1 x2 x3, means(0.34,0.39,0.19) corr(1,0.027,-0.008\0.027,1,-0.152\-0.008,-0.152,1) n(16000) 

/*
for the 1st scenario, we ran the same command, but we assumed that the binary variables where independent, i.e.

rbinary x1 x2 x3, means(0.34,0.39,0.19) corr(1,0,0\0,1,0\0,0,1) n(16000) 
*/


* Generate interactions
gen x12=x1*x2
gen x13=x1*x3
gen x23=x2*x3
gen x123=x1*x2*x3

* simulate data using Weibull
survsim stime, distribution(weibull) lambdas(.0066) gamma(2) ///
covariates(x1 .36 x2 .29 x3 .41 x12 -.27 x13 -.23 x23 -.24 x123 .92)
gen cens = runiform()*6
gen died = cens>stime
replace stime = cens if cens<stime
stset stime, failure(died = 1)

* analyse data using Cox regression
stcox x1 x2 x3 x12 x13 x23 x123, nohr

* Compute RERI with 3-way interactions
nlcom TotRERI3: exp(_b[x1]+_b[x2]+_b[x3]+_b[x12]+_b[x13]+_b[x23]+_b[x123])-exp(_b[x1])-exp(_b[x2])-exp(_b[x3])+2
return scalar TotRERI3=el(r(b),1,1)
return scalar seTotRERI3=sqrt(el(r(V),1,1))

nlcom RERI3:    exp(_b[x1]+_b[x2]+_b[x3]+_b[x12]+_b[x13]+_b[x23]+_b[x123])-exp(_b[x1]+_b[x2]+_b[x12])-exp(_b[x1]+_b[x3]+_b[x13])-exp(_b[x2]+_b[x3]+_b[x23])+exp(_b[x1])+exp(_b[x2])+exp(_b[x3])-1
return scalar RERI3=el(r(b),1,1)
return scalar seRERI3=sqrt(el(r(V),1,1))

*We compute 2-way interactions, given the 3rd risk factor is absent
*RERI(x1,x2/x3=0)
nlcom RERI2_x12: exp(_b[x1]+_b[x2]+_b[x12])-exp(_b[x1])-exp(_b[x2])+1
return scalar RERI2_x12=el(r(b),1,1)
return scalar seRERI2_x12=sqrt(el(r(V),1,1))
*RERI(x1,x3/x2=0)
nlcom RERI2_x13: exp(_b[x1]+_b[x3]+_b[x13])-exp(_b[x1])-exp(_b[x3])+1
return scalar RERI2_x13=el(r(b),1,1)
return scalar seRERI2_x13=sqrt(el(r(V),1,1))
*RERI(x2,x3/x1=0)
nlcom RERI2_x23: exp(_b[x2]+_b[x3]+_b[x23])-exp(_b[x2])-exp(_b[x3])+1
return scalar RERI2_x23=el(r(b),1,1)
return scalar seRERI2_x23=sqrt(el(r(V),1,1))

*We compute 2-way interactions, given the 3rd risk factor is present
*RERI(x1,x2/x3=1)
nlcom RERI2_x12e: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x12]+_b[x13]+_b[x23]+_b[x123])-exp(_b[x1]+_b[x3]+_b[x13])-exp(_b[x2]+_b[x3]+_b[x23])+exp(_b[x3]))/exp(_b[x3])
return scalar RERI2_x12e=el(r(b),1,1)
return scalar seRERI2_x12e=sqrt(el(r(V),1,1))
*RERI(x1,x3/x2=1)
nlcom RERI2_x13e: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x12]+_b[x13]+_b[x23]+_b[x123])-exp(_b[x1]+_b[x2]+_b[x12])-exp(_b[x2]+_b[x3]+_b[x23])+exp(_b[x2]))/exp(_b[x2])
return scalar RERI2_x13e=el(r(b),1,1)
return scalar seRERI2_x13e=sqrt(el(r(V),1,1))
*RERI(x2,x3/x1=1)
nlcom RERI2_x23e: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x12]+_b[x13]+_b[x23]+_b[x123])-exp(_b[x1]+_b[x2]+_b[x12])-exp(_b[x1]+_b[x3]+_b[x13])+exp(_b[x1]))/exp(_b[x1])
return scalar RERI2_x23e=el(r(b),1,1)
return scalar seRERI2_x23e=sqrt(el(r(V),1,1))
end

*run simulations and collect both point estimates and SEs.
simulate TotRERI3=r(TotRERI3) RERI3=r(RERI3) RERI2_x12=r(RERI2_x12) RERI2_x13=r(RERI2_x13) ///
RERI2_x23=r(RERI2_x23) RERI2_x12e=r(RERI2_x12e) RERI2_x13e=r(RERI2_x13e) RERI2_x23e=r(RERI2_x23e) ///
seTotRERI3=r(seTotRERI3) seRERI3=r(seRERI3) seRERI2_x12=r(seRERI2_x12) seRERI2_x13=r(seRERI2_x13) ///
seRERI2_x23=r(seRERI2_x23) seRERI2_x12e=r(seRERI2_x12e) seRERI2_x13e=r(seRERI2_x13e) ///
seRERI2_x23e=r(seRERI2_x23e), dots(1) reps(1000): MGsim

summarize


/* Bias */
gen true_TotRERI3=exp(.36 +.29 + .41 -.27 -.23 -.24 + .92)-exp(.36)-exp(.29)-exp(.41)+2
generate bias_TotRERI3 = TotRERI3 - true_TotRERI3

gen true_RERI3=exp(.36 +.29 + .41 -.27 -.23 -.24 + .92)-exp(.36+.29-.27)-exp(.36+.41-.23)-exp(.29+.41-.24)+exp(.36)+exp(.29)+exp(.41)-1
generate bias_RERI3 = RERI3 - true_RERI3

gen true_RERI2_x12=exp(.36+.29-.27)-exp(.36)-exp(.29)+1
generate bias_RERI2_x12 = RERI2_x12 - true_RERI2_x12

gen true_RERI2_x13=exp(.36+.41-.23)-exp(.36)-exp(.41)+1
generate bias_RERI2_x13 = RERI2_x13 - true_RERI2_x13

gen true_RERI2_x23=exp(.29 + .41 -.24)-exp(.29)-exp(.41)+1
generate bias_RERI2_x23 = RERI2_x23 - true_RERI2_x23

gen true_RERI2_x12e=(exp(.36 +.29 + .41 -.27 -.23 -.24 + .92)-exp(.36+.41-.23)-exp(.29+.41-.24)+exp(.41))/exp(.41)
generate bias_RERI2_x12e = RERI2_x12e - true_RERI2_x12e

gen true_RERI2_x13e=(exp(.36 +.29 + .41 -.27 -.23 -.24 + .92)-exp(.36+.29-.27)-exp(.29+.41-.24)+exp(.29))/exp(.29)
generate bias_RERI2_x13e = RERI2_x13e - true_RERI2_x13e

gen true_RERI2_x23e=(exp(.36 +.29 + .41 -.27 -.23 -.24 + .92)-exp(.36 +.29 -.27)-exp(.36+.41-.23)+exp(.36))/exp(.36)
generate bias_RERI2_x23e = RERI2_x23e - true_RERI2_x23e

summarize true_TotRERI3 true_RERI3 true_RERI2_x12 true_RERI2_x13 true_RERI2_x23 true_RERI2_x12e true_RERI2_x13e true_RERI2_x23e
summarize TotRERI3 RERI3 RERI2_x12 RERI2_x13 RERI2_x23 RERI2_x12e RERI2_x13e RERI2_x23e
summarize bias_TotRERI3 bias_RERI3 bias_RERI2_x12 bias_RERI2_x13 bias_RERI2_x23 bias_RERI2_x12e bias_RERI2_x13e bias_RERI2_x23e

summarize bias_TotRERI3
di r(mean)/true_TotRERI3
summarize bias_RERI3
di r(mean)/true_RERI3
summarize bias_RERI2_x12
di r(mean)/true_RERI2_x12
summarize bias_RERI2_x13
di r(mean)/true_RERI2_x13
summarize bias_RERI2_x23
di r(mean)/true_RERI2_x23
summarize bias_RERI2_x12e
di r(mean)/true_RERI2_x12e
summarize bias_RERI2_x13e
di r(mean)/true_RERI2_x13e
summarize bias_RERI2_x23e
di r(mean)/true_RERI2_x23e

/* rMSE */

generate rmse_TotRERI3 = (TotRERI3 - true_TotRERI3)^2
generate rmse_RERI3 = (RERI3 - true_RERI3)^2
generate rmse_RERI2_x12 = (RERI2_x12 - true_RERI2_x12)^2
generate rmse_RERI2_x13 = (RERI2_x13 - true_RERI2_x13)^2
generate rmse_RERI2_x23 = (RERI2_x23 - true_RERI2_x23)^2
generate rmse_RERI2_x12e = (RERI2_x12e - true_RERI2_x12e)^2
generate rmse_RERI2_x13e = (RERI2_x13e - true_RERI2_x13e)^2
generate rmse_RERI2_x23e = (RERI2_x23e - true_RERI2_x23e)^2

summarize rmse_TotRERI3 rmse_RERI3 rmse_RERI2_x12 rmse_RERI2_x13 rmse_RERI2_x23 rmse_RERI2_x12e rmse_RERI2_x13e rmse_RERI2_x23e



/* Coverage */
generate cov_TotRERI3 = (TotRERI3 + invnorm(0.975)*seTotRERI3>true_TotRERI3 & TotRERI3 - invnorm(0.975)*seTotRERI3<true_TotRERI3)
generate cov_RERI3 = (RERI3 + invnorm(0.975)*seRERI3>true_RERI3 & RERI3 - invnorm(0.975)*seRERI3<true_RERI3)
generate cov_RERI2_x12 = (RERI2_x12 + invnorm(0.975)*seRERI2_x12>true_RERI2_x12 & RERI2_x12 - invnorm(0.975)*seRERI2_x12<true_RERI2_x12)
generate cov_RERI2_x13 = (RERI2_x13 + invnorm(0.975)*seRERI2_x13>true_RERI2_x13 & RERI2_x13 - invnorm(0.975)*seRERI2_x13<true_RERI2_x13)
generate cov_RERI2_x23 = (RERI2_x23 + invnorm(0.975)*seRERI2_x23>true_RERI2_x23 & RERI2_x23 - invnorm(0.975)*seRERI2_x23<true_RERI2_x23)
generate cov_RERI2_x12e = (RERI2_x12e + invnorm(0.975)*seRERI2_x12e>true_RERI2_x12e & RERI2_x12e - invnorm(0.975)*seRERI2_x12e<true_RERI2_x12e)
generate cov_RERI2_x13e = (RERI2_x13e + invnorm(0.975)*seRERI2_x13e>true_RERI2_x13e & RERI2_x13e - invnorm(0.975)*seRERI2_x13e<true_RERI2_x13e)
generate cov_RERI2_x23e = (RERI2_x23e + invnorm(0.975)*seRERI2_x23e>true_RERI2_x23e & RERI2_x23e - invnorm(0.975)*seRERI2_x23e<true_RERI2_x23e)

summarize cov_TotRERI3 cov_RERI3 cov_RERI2_x12 cov_RERI2_x13 cov_RERI2_x23 cov_RERI2_x12e cov_RERI2_x13e cov_RERI2_x23e
