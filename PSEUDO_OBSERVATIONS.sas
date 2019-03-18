*------------------------------------------------------------------------------;
*-SAS macro for calculation of pseudo-observations from Kaplan-Meier estimates-;
*------------------------------------------------------------------------------;
%MACRO PSUEDO_OBS_KM(
	DATA_IN=, 	/* dataset for analysis */
	DATA_OUT=, 	/* dataset which includes pseudo-observations*/
	ID=, 		/* unique identifier for subjects */
	T=, 		/* outcome variable: time to event */
	D=, 		/* outcome variable: 1 uncensored, 0 censored */
	TIMEPOINT=	/* time point for analysis*/
);
proc sort data=&DATA_IN.;
by &T.;
run;
data &DATA_IN.;
set &DATA_IN.;
seq=_n_;
call symput ('samplesize',_n_);
run;
proc iml symsize= 370000000 worksize= 280000000;
use &DATA_IN.;
read all var{&ID.} 		into id;
read all var{&T.} 		into t;
read all var{&D.} 		into d;

start km(t,d);
	n	  = nrow(t);
	t1 	  = shape(t,n,n);
	t2 	  = t1`;
	atrisk 	  = (t2 >= t1);
	n_atrisk  = atrisk[+,]`;
	if d[n,1] =0 then do;
		one = j(n,1,1);
		s = one-d/n_atrisk;
		log_s = log(s);
		km_total = exp(cusum(log_s));
	end;
	else do;
		one = j(n-1,1,1);
		s = one-d[1:n-1,]/n_atrisk[1:n-1,];
		log_s = log(s);
		km_ = exp(cusum(log_s));
		km_total = (km_//0);
	end;
	use 	  = (t <= &TIMEPOINT.);
	one 	  = j(n,1,1);
	usekm 	  = use#km_total + 2*(one-use);
	km	  = min(usekm);
	return(km);
finish;

start jack_knife(i) global(t,d,km_total);
	n	= nrow(t);
	t_i	= remove(t,i)`;
	d_i	= remove(d,i)`;
	pseudo  = n*km_total-(n-1)*km(t_i,d_i);
	return(pseudo);
finish;

km_total	= km(t,d);
n		= nrow(t);
seq 		= t(shape(1:n, 1)); 
pseudo 		= apply("jack_knife", seq);

read all var{&T.} 		into t;
read all var{&D.} 		into d;
result	  =(id||t||d||pseudo);
create &DATA_OUT. from result [colname={"&ID.","&T.","&D.",'PSEUDO_OBS'}];
append from result;
quit;
%MEND;

*-----------------------------------------------------------------------------------------;
*-SAS macro for calculation of pseudo-observations from estimates of cumulative incidence-;
*-----------------------------------------------------------------------------------------;
%MACRO PSUEDO_OBS_AJ(
	DATA_IN=, 	/* dataset for analysis */
	DATA_OUT=, 	/* dataset which includes pseudo-observations*/
	ID=, 		/* unique identifier for subjects */
	T=, 		/* outcome variable: time to event */
	D=, 		/* outcome variable: 1 uncensored, 0 censored */
	EPSILON=, 	/* outcome variable: 1 event of interest, 2 competing risks */
	TIMEPOINT=	/* time point for analysis*/
);
proc sort data=&DATA_IN.;
by &T.;
run;
data &DATA_IN.;
set &DATA_IN.;
seq=_n_;
call symput ('samplesize',_n_);
run;
proc iml symsize= 370000000 worksize= 280000000;
use &DATA_IN.;
read all var{&ID.} 		into id;
read all var{&T.} 		into t;
read all var{&D.} 		into d;
read all var{&EPSILON.} 	into epsilon;

start ci(t,d,epsilon);
	n	  = nrow(t);
	t1 	  = shape(t,n,n);
	t2 	  = t1`;
	atrisk 	  = (t2 >= t1);
	n_atrisk  = atrisk[+,]`;
	if d[n,1] =0 then do;
		one = j(n,1,1);
		s = one-d/n_atrisk;
		log_s = log(s);
		km_total = exp(cusum(log_s));
	end;
	else do;
		one = j(n-1,1,1);
		s = one-d[1:n-1,]/n_atrisk[1:n-1,];
		log_s = log(s);
		km_ = exp(cusum(log_s));
		km_total = (km_//0);
	end;
	d1 	  = d#(epsilon=1);
	km_total_ = km_total[1:n-1,];
	km_total_ = (1//km_total_);
	ci_total_ = km_total_#(d1/n_atrisk);
	ci_total  = cusum(ci_total_);
	one 	  = j(n,1,1);
	use 	  = (t <= &TIMEPOINT.);
	useci 	  = use#ci_total;
	ci	  = max(useci);
	return(ci);
finish;

start jack_knife(i) global(t,d,epsilon,ci_total);
	n		= nrow(t);
	t_i		= remove(t,i)`;
	d_i		= remove(d,i)`;
	epsilon_i	= remove(epsilon,i)`;
	pseudo  	= n*ci_total-(n-1)*ci(t_i,d_i,epsilon_i);
	return(pseudo);
finish;

ci_total	= ci(t,d,epsilon);
n		= nrow(t);
seq 		= t(shape(1:n, 1)); 
pseudo 		= apply("jack_knife", seq);

read all var{&T.} 		into t;
read all var{&D.} 		into d;
read all var{&EPSILON.} 	into e;
result	  =(id||t||d||e||pseudo);
create &DATA_OUT. from result [colname={"&ID.","&T.","&D.","&EPSILON.",'PSEUDO_OBS'}];
append from result;
quit;
%MEND;

*-------------------------------------------------------;
*-SAS macro for g-estimation using pseudo-observations--;
*-of a risk difference, a risk ratio and a hazard ratio-;
*-------------------------------------------------------;
%MACRO G_ESTIMATION_PSEUDO_OBSERVATIONS(
	DATA_IN=, 	/* dataset for analysis */
	DATA_OUT=, 	/* dataset which includes estimates and 95% confidence intervals*/
	ID=, 		/* unique identifier for subjects */
	TIMEPOINT=,	/* maximum time point for analysis*/
	N_TIMEPOINT=,	/* number of time points for analysis*/
	T=, 		/* outcome variable: time to event */
	D=, 		/* outcome variable: 1 uncensored, 0 censored */
	EPSILON=, 	/* outcome variable: 1 event of interest, 2 competing risks */
	A=, 		/* exposure variable: 1 exposed, 0 not exposed */
	L_LIST=		/* full list of covariates to include in the exposure and outcome models */
);

%if &EPSILON ne %then %do;
	%put G-ESTIMATION OF EFFECTS ON CUMULATIVE INCIDENCE FUNCTION AT &TIMEPOINT. OR SUBDISTRIBUTION HAZARD FUNCTION;
	%put OUTCOME: &T., &D., AND &EPSILON.;
	%put EXPOSURE: &A.;
	%put ADJUSTED FOR &L_LIST.;
	option nonotes;
	ods exclude all;

	%do m=1 %to &N_TIMEPOINT.;
		%let tp=%sysevalf(&m*&TIMEPOINT./&N_TIMEPOINT.);
		%PSUEDO_OBS_AJ(DATA_IN=&DATA_IN., DATA_OUT=TMP&m., ID=&ID., T=&T., D=&D., EPSILON=&EPSILON., TIMEPOINT=&tp.);
		data TMP&m.;set TMP&m.;y=PSEUDO_OBS;run;
		proc sort data=TMP&m.;by &ID.;run;
		proc sort data=&DATA_IN.;by &ID.;run;
		data TMP&m.;merge TMP&m. &DATA_IN.;by &ID.;run;
	%end;
%end;
%else %do;
	%put G-ESTIMATION OF EFFECTS ON SURVIVAL FUNCTION AT &TIMEPOINT. OR HAZARD FUNCTION;
	%put OUTCOME: &T., &D.;
	%put EXPOSURE: &A.;
	%put ADJUSTED FOR &L_LIST.;
	option nonotes;
	ods exclude all;

	%do m=1 %to &N_TIMEPOINT.;
		%let tp=%sysevalf(&m*&TIMEPOINT./&N_TIMEPOINT.);
		%PSUEDO_OBS_KM(DATA_IN=&DATA_IN., DATA_OUT=TMP&m., ID=&ID., T=&T., D=&D., TIMEPOINT=&tp.);
		data TMP&m.;set TMP&m.;y=1-PSEUDO_OBS;run;
		proc sort data=TMP&m.;by &ID.;run;
		proc sort data=&DATA_IN.;by &ID.;run;
		data TMP&m.;merge TMP&m. &DATA_IN.;by &ID.;run;
	%end;
%end;

proc iml symsize=370000000 worksize= 280000000;
%do m=1 %to &N_TIMEPOINT. %by 1;
	use TMP&m.;
	read all var{y} 		into y;
	read all var{y} 		into y_0  where(&A.=0);
	read all var{y} 		into y_1  where(&A.=1);
	read all var{&A.} 		into a;

	%let cnt = 0;
	%let unique =;
	%let count = 0;

	%if &L_LIST ne %then %do;
	    %let word = %scan(&L_LIST,1);
	    %do %while (&word ne);
	        %let cnt = %eval(&cnt+1);
	        %let word = %scan(&L_LIST,&cnt+1);
	    %end;
	%end;

	%do loop1=1 %to &cnt;
	    %let word = %scan(&L_LIST,&loop1);
	    %if &count = 0 %then %do;
	        %let unique = &word;
	        %let count = 1;
			read all var{&word.}		into l;
			read all var{&word.} 		into l_0 where(&A.=0);
			read all var{&word.}		into l_1 where(&A.=1);
	    %end;
	    %else %do;
	        %do loop2=1 %to &count;
	            %if &word = %scan(&unique,&loop2) %then %let word=;
	        %end;
	        %if &word ne %then %do;
	            %let unique = &unique &word;
	            %let count = %eval(&count+1);
			read all var{&word.}		into tmp;
			l  =(l ||tmp);
			read all var{&word.} 		into tmp where(&A.=0);
			l_0 =(l_0 ||tmp);
			read all var{&word.}		into tmp where(&A.=1);
			l_1 =(l_1 ||tmp);
	        %end;
	    %end;
	%end;

	%EXPECTATION(A,L,Y_0,L_0,Y_1,L_1);
	Y_repeated=(Y_repeated//Y);
	A_repeated=(A_repeated//A);
	EU_repeated=(EU_repeated//EU);
	EA_repeated=(EA_repeated//EA);
	XITA_repeated=(XITA_repeated//XITA);
%end;

%GESTIMATION_RD(Y,A,EU,EA);
%GESTIMATION_RR(Y,A,EU,EA);
%GESTIMATION_HR(Y_repeated,A_repeated,EU_repeated,EA_repeated,XITA_repeated);
tmp1=(risk_difference//risk_ratio//hazard_ratio);
tmp2=(cil_risk_difference//cil_risk_ratio//cil_hazard_ratio);
tmp3=(ciu_risk_difference//ciu_risk_ratio//ciu_hazard_ratio);
tmp4=(p_risk_difference//p_risk_ratio//p_hazard_ratio);
res=(tmp1||tmp2||tmp3||tmp4);
create &DATA_OUT. from res [colname={'ESTIMATE','LOWER_95%CI','UPPER_95%CI','P'}];
append from res;
quit;
data &DATA_OUT.;
set &DATA_OUT.;
if _n_=1 then LABEL="RISK DIFFERENCE";
if _n_=2 then LABEL="RISK RATIO";
if _n_=3 then LABEL="HAZARD RATIO";
run;
option notes;
ods exclude none;
proc print data=&DATA_OUT.;
run;
%MEND;

*------------------------------------------------------;
*------------------------------------------------------;
*------------------------------------------------------;
%MACRO EXPECTATION(A,L,Y_0,L_0,Y_1,L_1);
start ml_alpha(alpha) global(&A.,x_a,one);
e_x_alpha = exp(x_a*alpha`);
ea	  = e_x_alpha/(one+e_x_alpha);
ret	  = &A.#log(ea) + (one-&A.)#log(one-ea);
ret	  = ret[+,+];
return(ret);
finish ml_alpha;

one	  = j(nrow(&A.),1,1);
x_a	  = (one||&L.);
p_a	  = ncol(x_a);
alpha_0	  = j(1,p_a,0);
optn	  = (p_a||0);
call nlpnra(rc,alpha,"ml_alpha",alpha_0,optn);
e_x_alpha = exp(x_a*alpha`);
ea	  = e_x_alpha/(one+e_x_alpha);

*-------------------------------------------;
one_	  = j(nrow(&L_0.),1,1);
x_g_	  = (one_||&L_0.);
gamma     = ginv(x_g_`*x_g_)*x_g_`*&Y_0.;
one	  = j(nrow(&A.),1,1);
x_g	  = (one||&L.);
eu        = x_g*gamma;

*-------------------------------------------;
start ilink(x);
one	  = j(nrow(x),1,1);
y	  = exp(-exp(-x));
return y;
finish ilink;

start ml_ita(ita) global(&Y_1.,x_y1_,one_);
xita	  = x_y1_*ita`;
ey1	  = ilink(xita);
one_y1	  = one_-&Y_1.;
ret	  = -(one_y1-ey1)#(one_y1-ey1);
ret	  = ret[+,+];
return(ret);
finish ml_ita;

one_	  = j(nrow(&L_1.),1,1);
x_y1_	  = (one_||&L_1.);
p_y1_	  = ncol(x_y1_);
ita_0	  = j(1,p_y1_,0);
optn	  = (p_y1_||0);
call nlpnra(rc,ita,"ml_ita",ita_0,optn);
one	  = j(nrow(&L.),1,1);
x_y1	  = (one||&L.);
xita	  = x_y1*ita`;
%MEND;

*------------------------------------------------------;
*------------------------------------------------------;
*------------------------------------------------------;
%MACRO GESTIMATION_RD(Y,A,EU,EA); 
risk_difference	  = ginv((&A.-&EA.)`*(&A.-&EA.))*(&A.-&EA.)`*(&Y.-&EU.);
score	  = (&A.-&EA.)#(&Y.-risk_difference*&A.-&EU.);
var_score = score`*score/nrow(score);
hessian_  = (&A.-&EA.)#(-&A.);
hessian	  = hessian_[+,+]/nrow(score);
var	  = ginv(hessian)`*var_score*ginv(hessian)/nrow(score);
se	  = sqrt(var);
ciu_risk_difference = risk_difference+1.96*se;
cil_risk_difference = risk_difference-1.96*se;
chisq	  = risk_difference#risk_difference/var;
p_risk_difference	  = 1-probchi(chisq,1);
%MEND;

*------------------------------------------------------;
*------------------------------------------------------;
*------------------------------------------------------;
%MACRO GESTIMATION_RR(Y,A,EU,EA); 
start ee_risk_ratio(beta) global(&Y.,&EU.,&A.,&EA.);
score	  =(&A.-&EA.)`*(&Y.#exp(-beta*&A.)-&EU.);
return(score);
finish ee_risk_ratio;

beta_0	  = 0;
p_b	  = ncol(beta_0);
optn	  = (p_b||2);
call nlphqn(rc,beta,"ee_risk_ratio",beta_0,optn);
score	  = (&A.-&EA.)#(&Y.#exp(-beta*&A.)-&EU.);
var_score = score`*score/nrow(score);
hessian_  = (&A.-&EA.)#&Y.#(-&A.)#exp(-beta*&A.);
hessian	  = hessian_[+,+]/nrow(score);
risk_ratio = exp(beta);
var	  = ginv(hessian)`*var_score*ginv(hessian)/nrow(score);
se	  = sqrt(var);
ciu_risk_ratio  = exp(beta+(1.96*se));
cil_risk_ratio  = exp(beta-(1.96*se));
chisq	  = beta#beta/var;
p_risk_ratio	  = 1-probchi(chisq,1);
%MEND;

*------------------------------------------------------;
*------------------------------------------------------;
*------------------------------------------------------;
%MACRO GESTIMATION_HR(Y,A,EU,EA,XITA); 
start ilink(x);
one	  = j(nrow(x),1,1);
y	  = exp(-exp(-x));
return y;
finish ilink;

start ee_hazard_ratio(m_beta) global(&Y.,&EU.,&A.,&EA.,&XITA.);
tmp	  = &XITA. - m_beta*&A.;
one	  = j(nrow(&A.),1,1);
one_y	  = one-&Y.;
one_eu	  = one-&EU.;
h	  = (one-&A.)#one_y + &A.#ilink(tmp);
score	  = (&A.-&EA.)`*(h-one_eu);
return(score);
finish ee_hazard_ratio;

m_beta_0	  = 0;
p_b	  = ncol(m_beta_0);
optn	  = (p_b||2);
call nlphqn(rc,m_beta,"ee_hazard_ratio",m_beta_0,optn);
beta	  = -m_beta;

one	  = j(nrow(&A.),1,1);
one_y	  = one-&Y.;
one_eu	  = one-&EU.;
tmp	  = &XITA. - m_beta*&A.;
h	  = (one-&A.)#one_y + &A.#ilink(tmp);
score	  = (&A.-&EA.)#(h - one_eu);
var_score = score`*score/nrow(score);
call nlpfdd(f,g,hessian_,"ee_hazard_ratio",beta);
hessian  = hessian_/nrow(score);

hazard_ratio 		  = exp(beta);
var	  = ginv(hessian)`*var_score*ginv(hessian)/nrow(score);
se	  = sqrt(var);
ciu_hazard_ratio = exp(beta+1.96*se);
cil_hazard_ratio = exp(beta-1.96*se);
chisq	  = beta#beta/var;
p_hazard_ratio	  = 1-probchi(chisq,1);
%MEND;
