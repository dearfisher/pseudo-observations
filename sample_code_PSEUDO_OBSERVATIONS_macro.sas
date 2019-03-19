*--------------------------------------------------------------------------------------------;
*-This is sample codes for calculating pseudo-observations based on the Kaplan-Meier---------;
*-estimator or the Aalen-Johansen estimator using PSUEDO_OBS_KM or PSUEDO_OBS_AJ macro-------;
*--------------------------------------------------------------------------------------------;
%let path   = C:\Users\Shiro\Dropbox\(ÅB•_•°)…\2016\Pseudo-observations\GitHub;
%inc "&path.\PSEUDO_OBSERVATIONS.sas";

*--------------------------------------------------------------------------------------------;
*-BMT dataset, originally used in SAS/STAT(R) 13.1 User's Guide------------------------------;
*-Example 71.15 Analysis of Competing-Risks Data---------------------------------------------;
*--------------------------------------------------------------------------------------------;
%inc "&path.\data_bmt.sas";

*--------------------------------------------------------------------------------------------;
*-Analysis in Example 71.15 Analysis of Competing-Risks Data---------------------------------;
*--------------------------------------------------------------------------------------------;
data Risk;
   Disease=1; output;
   Disease=2; output;
   Disease=3; output;
   format Disease DiseaseGroup.;
   run;
ods graphics on;
proc phreg data=BMT plots(overlay=stratum)=cif;
   class Disease (order=internal ref=first);
   model T*Status(0)=Disease / eventcode=1;
   Hazardratio 'Pairwise' Disease / diff=pairwise;
   baseline covariates=Risk out=out1 cif=_all_ / seed=191;
run;

*--------------------------------------------------------------------------------------------;
*-Estimation of a risk ratio using PSEUDO_OBSERVATONS macro----------------------------------;
*--------------------------------------------------------------------------------------------;
data BMT;
set BMT;
if Status=0 then D=0;		/*Creating a censoring variable*/
else D=1;
ID=_n_;				/*Creating an id variable for merge*/
run;
%PSUEDO_OBS_AJ(
	DATA_IN=BMT, 		/* dataset for analysis */
	DATA_OUT=BMT_OUT, 	/* dataset which includes pseudo-observations*/
	ID=ID, 			/* unique identifier for subjects */
	T=T, 			/* outcome variable: time to event */
	D=D, 			/* outcome variable: 1 uncensored, 0 censored */
	EPSILON=Status, 	/* outcome variable: 1 event of interest, 2 competing risks */
	TIMEPOINT=1825		/* time point for analysis*/
);
proc sort data=BMT;
by ID;
run;
proc sort data=BMT_OUT;
by ID;
run;
data BMT;
merge BMT BMT_OUT;
by ID;
run;
proc genmod data=BMT;
class ID;
model pseudo_obs=Disease / dist=normal link=log;
repeated subject=ID / corr=unstr;
run;
