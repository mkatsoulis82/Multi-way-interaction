*  Implementation in STATA
/*
DESCRIPTION
time; survival time
death (outcome); 0-->alive, 1-->dead
x1 (non adherence to Mediterranean diet - 1st risk factor); 
0--> high adherence to Mediterranean diet, 1--> low adherence to Mediterranean diet

x2 (being obese - 2nd risk factor); 
0--> not being obese (BMI<30), 1--> being obese (BMI>=30)

x3 (smoking status - 3rd risk factor); 
0--> never/former smoker, 1--> current smoker

u1 (age             - 1st confounder); continuous in years
u2 (education level - 2nd confounder); categorical in 4 levels
*/

* IMPLEMENTATION
* At first we compute the product terms
gen x1x2=x1*x2
gen x1x3=x1*x3
gen x2x3=x2*x3
gen x1x2x3=x1*x2*x3

*We run the Cox model
stset time, failure(death)

stcox x1 x2 x3 x1x2 x1x3 x2x3 x1x2x3 u1 i.u2


*We compute TotRERI3
nlcom TotRERI3: exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1])-exp(_b[x2])-exp(_b[x3])+2

*We compute RERI3
nlcom RERI3:    exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x1]+_b[x3]+_b[x1x3])-exp(_b[x2]+_b[x3]+_b[x2x3])+exp(_b[x1])+exp(_b[x2])+exp(_b[x3])-1

*We compute 2-way interactions, given the 3rd risk factor is absent
*RERI(x1,x2/x3=0)
nlcom RERI2_x1_x2_given_x3is0: exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x1])-exp(_b[x2])+1
*RERI(x1,x3/x2=0)
nlcom RERI2_x1_x3_given_x2is0: exp(_b[x1]+_b[x3]+_b[x1x3])-exp(_b[x1])-exp(_b[x3])+1
*RERI(x2,x3/x1=0)
nlcom RERI2_x2_x3_given_x1is0: exp(_b[x2]+_b[x3]+_b[x2x3])-exp(_b[x2])-exp(_b[x3])+1

*We compute 2-way interactions, given the 3rd risk factor is present

*RERI(x1,x2/x3=1)
nlcom RERI2_x1_x2_given_x3is1: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x3]+_b[x1x3])-exp(_b[x2]+_b[x3]+_b[x2x3])+exp(_b[x3]))/exp(_b[x3])

*RERI(x1,x3/x2=1)
nlcom RERI2_x1_x3_given_x2is1: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x2]+_b[x3]+_b[x2x3])+exp(_b[x2]))/exp(_b[x2])

*RERI(x2,x3/x1=1)
nlcom RERI2_x2_x3_given_x1is1: (exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x1]+_b[x3]+_b[x1x3])+exp(_b[x1]))/exp(_b[x1])

*The same formulae for all these RERIs are used when running logistic regression
* CHECK FOR QUALITATIVE INTERACTION

*We run again the Cox model
stset time, failure(death)

stcox x1 x2 x3 x1x2 x1x3 x2x3 x1x2x3 u1 i.u2

* To check whether the risk for x1 is increasing across strata of x2,x3, we have to examine whether the following quantities are positive (i.e. >0)
* 1a) to see if RR100>RR000, we check whether RR100-RR000>0
disp exp(_b[x1])-1
* 1b) to see if RR110>RR010, we check whether RR110-RR010>0
disp exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x2])
* 1c) to see if RR101>RR001, we check whether RR101-RR001>0
disp exp(_b[x1]+_b[x3]+_b[x1x3])-exp(_b[x3])
* 1d) to see if RR111>RR011, we check whether RR111-RR011>0
disp exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x2]+_b[x3]+_b[x2x3])


* To check whether the risk for x2 is increasing across strata of x1,x3
* 2a) to see if RR010>RR000, we check whether RR010-RR000>0
disp exp(_b[x2])-1
* 2b) to see if RR110>RR100, we check whether RR110-RR100>0
disp exp(_b[x1]+_b[x2]+_b[x1x2])-exp(_b[x1])
* 2c) to see if RR011>RR001, we check whether RR011-RR001>0
disp exp(_b[x2]+_b[x3]+_b[x2x3])-exp(_b[x3])
* 2d) to see if RR111>RR101, we check whether RR111-RR101>0
disp exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x3]+_b[x1x3])


* To check whether the risk for x3 is increasing across strata of x1,x2
* 2a) to see if RR001>RR000, we check whether RR001-RR000>0
disp exp(_b[x3])-1
* 2b) to see if RR101>RR100, we check whether RR101-RR010>0
disp exp(_b[x1]+_b[x3]+_b[x1x3])-exp(_b[x1])
* 2c) to see if RR011>RR010, we check whether RR011-RR100>0
disp exp(_b[x2]+_b[x3]+_b[x2x3])-exp(_b[x2])
* 2d) to see if RR111>RR110, we check whether RR111-RR110>0
disp exp(_b[x1]+_b[x2]+_b[x3]+_b[x1x2]+_b[x1x3]+_b[x2x3]+_b[x1x2x3])-exp(_b[x1]+_b[x2]+_b[x1x2])


