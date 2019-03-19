*--------------------------------------------------------------------------------------------;
*-This is sample codes for estimating adjusted risk difference, risk ratio and hazard ratio--;
*-of structural mean models using G_ESTIMATION_PSEUDO_OBSERVATIONS macro---------------------;
*--------------------------------------------------------------------------------------------;
%let path   = C:\Users\Shiro\Dropbox\(ÅB•_•°)…\2016\Pseudo-observations\GitHub;
%inc "&path.\G_ESTIMATION_PSEUDO_OBSERVATIONS.sas";

*--------------------------------------------------------------------------------------------;
*-DIABETES dataset, a simulated dataset of 1000 patients from the JDCS/J-EDIT----------------;
*--------------------------------------------------------------------------------------------;
proc format;
value DF 0='Censored'
         1='Event';
value EPSILONF 0='Censored'
               1='Stroke'
               2='Death';
value AF 0='High LTPA'
         1='Low LTPA';
data DIABETES;
infile "&path.\data_diabetes.csv" dlm="," dsd lrecl=32767;
input ID T D EPSILON A L1 L2 @@;
label T='Time to event in days';
label D='Censoring';
label EPSILON='Type of event';
label A='Exposure';
label L1='Age';
label L2='Woman';
format D DF. EPSILON EPSILONF. A AF.;
run;

*--------------------------------------------------------------------------------------------;
*-Estimation of causal effects of A adjusted for L1 and L2-----------------------------------;
*-using of G_ESTIMATION_PSEUDO_OBSERVATIONS G_ESTIMATION macro-------------------------------;
*--------------------------------------------------------------------------------------------;
%G_ESTIMATION_PSEUDO_OBSERVATIONS(
	DATA_IN=DIABETES, 		/* dataset for analysis */
	DATA_OUT=DIABETES_OUT, 		/* dataset which includes pseudo-observations*/
	ID=ID, 				/* unique identifier for subjects */
	TIMEPOINT=8,			/* maximum time point for analysis */
	N_TIMEPOINT=5,			/* number of time points for analysis */
	T=T, 				/* outcome variable: time to event */
	D=D, 				/* outcome variable: 1 uncensored, 0 censored */
	EPSILON=EPSILON,	 	/* outcome variable: 1 event of interest, 2 competing risks */
	A=A, 				/* exposure variable: 1 exposed, 0 not exposed */
	L_LIST=L1 L2
);

*--------------------------------------------------------------------------------------------;
*-Comparisons with generalized linear models and the Fine-Gray model-------------------------;
*--------------------------------------------------------------------------------------------;
%PSUEDO_OBS_AJ(
	DATA_IN=DIABETES,
	DATA_OUT=DIABETES_OUT,
	ID=ID,
	T=T,
	D=D,
	EPSILON=EPSILON,
	TIMEPOINT=8
);
proc sort data=DIABETES;
by ID;
run;
proc sort data=DIABETES_OUT;
by ID;
run;
data DIABETES;
merge DIABETES DIABETES_OUT;
by ID;
run;
proc genmod data=DIABETES;			/*risk difference*/
class ID;
model pseudo_obs=A L1 L2 / dist=normal link=identity;
repeated subject=ID / corr=unstr;
run;
proc genmod data=DIABETES;			/*risk ratio*/
class ID;
model pseudo_obs=A L1 L2 / dist=normal link=log;
repeated subject=ID / corr=unstr;
run;
proc phreg data=DIABETES;			/*hazard ratio*/
model T*EPSILON(0)=A L1 L2 / eventcode=1;
run;
