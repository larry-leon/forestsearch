/*******************************************
2022-7 with HELEN
overfitting
sensitivity analysis
**************************************************/
data _eff_cough; set lptda.adeff;
where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
	  paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (0, 4,8,12,16,20,24);
keep usubjid subjid trtp trtpn ADT fascnfl ontrtfl paramcd aval sex  avisitn base ablfl chg;
run;

proc sql;
create table __ot_24hr_pervisit as select 
	distinct avisitn, trtpn, count(distinct usubjid) as n, exp(mean(base)) as base_gm, exp(mean(aval)) as post_gm, exp(mean(chg)) as GMR
		from _eff_cough 
			where chg ne .
				group by avisitn, trtpn
					order by avisitn, trtpn;
quit;

data __outdesc; set __ot_24hr_pervisit;
	N =n;
	base_gm = round(base_gm, 0.01);
	post_gm = round(post_gm, 0.01);
	GMR = ROUND(GMR, 0.01);
	KEEP AVISITN TRTPN N BASE_GM POST_GM GMR;
	RUN;


proc mixed data=_eff_cough ;
	where avisitn in (4 8 12 16 20 24) ;
		class usubjid trtpn avisitn sex ;
		model chg= base trtpn avisitn trtpn*avisitn  sex /ddfm=kr solution ;
		repeated avisitn/subject=usubjid type=TOEP R rcorr;

		lsmestimate trtpn*avisitn 'Placebo, at Week 4'              [1, 1 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 4'        [1, 2 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 4'        [1, 2 1] [-1, 1 1]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 8'              [1, 1 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 8'        [1, 2 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 8'        [1, 2 2] [-1, 1 2]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 12'              [1, 1 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 12'        [1, 2 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 12'        [1, 2 3] [-1, 1 3]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 16'              [1, 1 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 16'        [1, 2 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 16'        [1, 2 4] [-1, 1 4]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 20'              [1, 1 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 20'        [1, 2 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 20'        [1, 2 5] [-1, 1 5]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 24'              [1, 1 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 24'        [1, 2 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 24'        [1, 2 6] [-1, 1 6]/exp cl;
		ods output LSMEstimates=_lsm;
run;



DATA __ot_outm; set _lsm;
where stmtno in (1,2, 4,5, 7,8, 10,11, 13,14, 16,17 );
if stmtno in (1,2) then avisitn=4;
if stmtno in (4,5) then avisitn=8;
if stmtno in (7,8) then avisitn=12;
if stmtno in (10,11) then avisitn=16;
if stmtno in (13,14) then avisitn=20;
if stmtno in (16,17) then avisitn=24;

if stmtno in (1,4,7,10,13,16) then trtpn=1;
if stmtno in (2,5,8,11,14,17) then trtpn=3;

keep label avisitn trtpn expestimate lowerexp upperexp ;
run;

DATA __ot_outdiff; set _lsm;
where stmtno in (3, 6, 9, 12, 15, 18 );
if stmtno in (3) then do; avisitn=4; trtpn=3; end;
if stmtno in (6) then do; avisitn=8;  trtpn=3; end;
if stmtno in (9) then do; avisitn=12;  trtpn=3; end;
if stmtno in (12) then do; avisitn=16;  trtpn=3; end;
if stmtno in (15) then do; avisitn=20; trtpn=3; end;
if stmtno in (18) then do; avisitn=24; trtpn=3; end;

modelrr = (expestimate-1)*100;
modelrrL = (lowerexp-1)*100;
modelrrU = (upperexp-1)*100;

keep label avisitn trtpn modelrr modelrrL modelrrU  ;
run;

data __part1; 
	merge __OUTDESC(IN=A) __ot_outm(in=b ) __ot_outdiff(in=c);
	if a;
	by avisitn trtpn;
	
	__modelGMR = round(expestimate, .01);
	__modelgmrL = round(lowerexp, .01);
	__modelgmrU = round(upperexp, .01);

	__modelrr = round(modelrr, .01);
	__modelrrL = round(modelrrL, .01);
	__modelrrU = round(modelrrU, .01);

	keep avisitn trtpn n base_gm post_gm gmr 
		__modelGMR __modelgmrL __modelgmrU 
		__modelrr __modelrrL __modelrrU;
	run;



proc sgplot data=__OUTDESC;
		   title1 'descriptive GMR over time';
		   scatter x=avisitn y=GMR / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=GMR / group=TRTPN groupdisplay=cluster clusterwidth=0.5;
  
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

proc sgplot data=__ot_outm;
		   title1 'model-based GMR over time -  trt avisit trt*avisit sex TOEP';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;
  
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;






/************************
MI
step 1: mcmc to get monotone missing pattern
Using a Markov Chain Monte Carlo method, make the observed dataset monotone-missing. 
This will be accomplished for each treatment group using "proc mi" within SAS 9.3 by utilizing the options "mcmc
chain=multiple impute=monotone;", in conjunction with all of the covariates (excluding
treatment) included in the primary analysis model. The random seed will be set equal to 653832. This step
will generate m monotone-missing datasets. 


step 2: MI
Applying parametric regression to the monotone-missing datasets, impute all the missing values. 
This will be accomplished for each treatment group using "proc mi" within SAS 9.3 utilizing the option
"monotone reg", in conjunction with all of the covariates (excluding treatment) included in the primary
analysis model. The random seed will be set equal to 653832. This step will generate m complete datasets.

step 3: ANOCOVA model for each imputed dataset

step 4: combine all using RUBIN's rule
***********************		*/

/*hypothetical*/
data _ontrt_eff_cough_abl; set _eff_cough;
if ablfl = 'Y' then output;
else if ontrtfl ='Y' and avisitn in (4,8,12,16,20,24) then output;
run;
proc sort data= _ontrt_eff_cough_abl;
  by usubjid trtpn sex avisitn;
run;
proc transpose data= _ontrt_eff_cough_abl out= _ontrt_eff_cough_wide (drop=_name_ _label_) prefix=y;
  by uSUBJID trtpn sex;
  var AVAL;
  id avisitn;
run;
proc sort data= _ontrt_eff_cough_wide;
  by TRTPN sex;
run;

* STEP 1: to get monotone missing pattern;
proc mi data=_ontrt_eff_cough_wide out=_ontrt_eff_mono nimpute=100 seed=7264030 noprint;
    by TRTPN sex;  		*model will be run for each combiniation of the varaibles in the BY STATEMENT;
    var  y0 y4 y8 y12 y16 y20 y24 ; 
    mcmc chain=multiple impute=monotone;
run; 

data _ontrt_eff_mono1; set _ontrt_eff_mono;
if y24=. then flag=1;
run;
proc sort data=_ontrt_eff_mono1 out=_ontrt_eff_miin;
    by  _Imputation_ trtpn;
run;

* STEP 2: to get imputed datasets;
proc mi data=_ontrt_eff_miin  out=_ontrt_eff_miout nimpute=1 seed=7264030  noprint ; 
	by  _Imputation_ trtpn; 	*model will be run for each combiniation of the varaibles in the BY STATEMENT;
							    * note: sparse data for treatment by strata combination, need to comment out trtpn for program to run properly;
	class sex;
	monotone reg;   			*linear regression, each variable in the VAR is imputed using all the preceeding variables in the VAR statement as covariate;
	var sex y0 y4 y8 y12 y16 y20 y24  ; 
run;

data _ontrt_temp_base; set _ontrt_eff_miout;
avisitn=0; base=y0; output;
run;
data _ontrt_temp_post; set _ontrt_eff_miout;
	avisitn=4; logaval=y4; chg=logaval-y0; output;
	avisitn=8; logaval=y8; chg=logaval-y0; output;
	avisitn=12; logaval=y12; chg=logaval-y0; output;
	avisitn=16; logaval=y16; chg=logaval-y0; output;
	avisitn=20; logaval=y20; chg=logaval-y0; output;
	avisitn=24; logaval=y24; chg=logaval-y0; output;
	run;

proc sort data=_ontrt_temp_base;
by _imputation_ usubjid;
proc sort data=_ontrt_temp_post;
by _imputation_ usubjid;
run;

	data _ontrt_temp; merge _ontrt_temp_base(keep=usubjid  _imputation_ base trtpn) 
							_ontrt_temp_post(keep=usubjid _imputation_ avisitn logaval chg sex);
	by  _imputation_ usubjid;
	run;

** STEP 3: ANCOVA;
*** note: before running the cLDA model, do the checking and create macro varaibles accordinly as done for the primary analysis;
proc mixed data=_ontrt_temp ;
	where avisitn in (4 8 12 16 20 24) ;
		class usubjid trtpn avisitn sex ;
		model chg= base trtpn avisitn trtpn*avisitn sex /ddfm=kr solution ;
		by _imputation_;
		repeated avisitn/subject=usubjid type=un R rcorr;

		lsmestimate trtpn*avisitn 'Placebo, at Week 4'              [1, 1 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 4'        [1, 2 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 4'        [1, 2 1] [-1, 1 1]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 8'              [1, 1 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 8'        [1, 2 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 8'        [1, 2 2] [-1, 1 2]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 12'              [1, 1 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 12'        [1, 2 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 12'        [1, 2 3] [-1, 1 3]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 16'              [1, 1 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 16'        [1, 2 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 16'        [1, 2 4] [-1, 1 4]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 20'              [1, 1 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 20'        [1, 2 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 20'        [1, 2 5] [-1, 1 5]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 24'              [1, 1 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 24'        [1, 2 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 24'        [1, 2 6] [-1, 1 6]/exp cl;
		ods output LSMEstimates=_ontrt_lsm;
run;


DATA _ontrt__ot_outm; set _ontrt_lsm;
where stmtno in (1,2, 4,5, 7,8, 10,11, 13,14, 16,17 );
if stmtno in (1,2) then avisitn=4;
if stmtno in (4,5) then avisitn=8;
if stmtno in (7,8) then avisitn=12;
if stmtno in (10,11) then avisitn=16;
if stmtno in (13,14) then avisitn=20;
if stmtno in (16,17) then avisitn=24;

if stmtno in (1,4,7,10,13,16) then trtpn=1;
if stmtno in (2,5,8,11,14,17) then trtpn=3;

*keep label avisitn trtpn expestimate lowerexp upperexp ;
run;

DATA _ontrt__ot_outdiff; set _ontrt_lsm;
where stmtno in (3, 6, 9, 12, 15, 18 );
if stmtno in (3) then do; avisitn=4; trtpn=3; end;
if stmtno in (6) then do; avisitn=8;  trtpn=3; end;
if stmtno in (9) then do; avisitn=12;  trtpn=3; end;
if stmtno in (12) then do; avisitn=16;  trtpn=3; end;
if stmtno in (15) then do; avisitn=20; trtpn=3; end;
if stmtno in (18) then do; avisitn=24; trtpn=3; end;

modelrr = (expestimate-1)*100;
modelrrL = (lowerexp-1)*100;
modelrrU = (upperexp-1)*100;
withinvar = stderr**2;
*keep label _label_ avisitn trtpn modelrr modelrrL modelrrU  ;
run;

** STEP 4: COMBINE BY RUBIN;
proc sort data=_ontrt__ot_outm out=_ontrt__ot_outm1;
by stmtno label _imputation_;
run;

proc sort data=_ontrt__ot_outdiff out=_ontrt__ot_outdiff1;
by stmtno label _imputation_;
run;

ODS OUTPUT PARAMETERESTIMATES=_ontrtmeans_rubin;
	PROC MIANALYZE DATA=_ontrt__ot_outm1 ;
	BY stmtno label;
	MODELEFFECTS estimate;
	STDERR stderr;
	RUN;
	ODS OUTPUT CLOSE;
	data _ontrtmeans_rubin1; set _ontrtmeans_rubin;
	expestimate = round(exp(estimate), .01);
	expupper = round(exp(uclmean), .01);
	explower = round(exp(lclmean), .01);
	run;

	data _ontrtmeans_rubin2;set _ontrtmeans_rubin1;
	if stmtno in (1,4,7,10,13,16) the trtpn=1;
	else trtpn=3;
	if stmtno in (1,2) then avisitn=4;
	else if stmtno in (4,5) then avisitn=8;
	else if stmtno in (7,8) then avisitn=12;
	else if stmtno in (10,11) then avisitn=16;
	else if stmtno in (13,14) then avisitn=20;
	else if stmtno in (16,17) then avisitn=24;
	keep stmtno label avisitn trtpn expestimate expupper explower;
run;
proc sgplot data=_ontrtmeans_rubin2;
		   title1 'model-based GMR over time - HYPOTHETICAL';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;
 
			xaxistable expestimate / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='GMR';
	
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

ODS OUTPUT PARAMETERESTIMATES=_ontrtdiffs_rubin;
	PROC MIANALYZE DATA=_ontrt__ot_outdiff1 ;
	BY stmtno label;
	MODELEFFECTS estimate;
	STDERR stderr;
	RUN;
	ODS OUTPUT CLOSE;

	data _ontrtdiffs_rubin1; set _ontrtdiffs_rubin;
	expestimate = exp(estimate);
	expupper = exp(uclmean);
	explower = exp(lclmean);
	modelrr = (expestimate-1)*100;
	modelrrL = (explower-1)*100;
	modelrrU = (expupper-1)*100;
	run;



/*******************
treatment policy
***********************/

data _eff_cough_abl; set _eff_cough;
if ablfl = 'Y' then output;
else if avisitn in (4,8,12,16,20,24) then output;
run;

proc sort data= _eff_cough_abl;
  by usubjid trtpn sex avisitn;
run;
proc transpose data= _eff_cough_abl out= _eff_cough_wide (drop=_name_ _label_) prefix=y;
  by uSUBJID trtpn sex;
  var AVAL;
  id avisitn;
run;
proc sort data= _ontrt_eff_cough_wide;
  by TRTPN sex;
run;

proc transpose data= _eff_cough_abl out= _eff_cough_wide_flag (drop=_name_ _label_) prefix=flagy;
  by uSUBJID trtpn sex;
  var ontrtfl;
  id avisitn;
run;

data _eff_cough_wide1; merge
	_eff_cough_wide _eff_cough_wide_flag;
	by usubjid;
	IF flagy24='Y' then do; 
		flagy4='Y';flagy8='Y';flagy12='Y'; flagy16='Y';flagy20='Y'; 
	end;
	else if flagy20='Y' then do;
		flagy4='Y';flagy8='Y';flagy12='Y'; flagy16='Y';
	end;
	else if flagy16='Y' then do;
		flagy4='Y';flagy8='Y';flagy12='Y';
	end;
	else if flagy12='Y' then do;
		flagy4='Y';flagy8='Y';
	end;
	else if flagy8='Y' then do;
		flagy4='Y';
	end;
	run;


* STEP 1: to get monotone missing pattern;
proc mi data=_ontrt_eff_cough_wide out=_ontrt_eff_mono nimpute=100 seed=7264030 noprint;
    by TRTPN sex;  		*model will be run for each combiniation of the varaibles in the BY STATEMENT;
    var  y0 y4 y8 y12 y16 y20 y24 ; 
    mcmc chain=multiple impute=monotone;
run; 

data _ontrt_eff_mono1; set _ontrt_eff_mono;
if y24=. then flag=1;
run;
proc sort data=_ontrt_eff_mono1 out=_ontrt_eff_miin;
    by  _Imputation_ trtpn;
run;

* STEP 2: to get imputed datasets;
proc mi data=_ontrt_eff_miin  out=_ontrt_eff_miout nimpute=1 seed=7264030  noprint ; 
	by  _Imputation_ trtpn; 	*model will be run for each combiniation of the varaibles in the BY STATEMENT;
							    * note: sparse data for treatment by strata combination, need to comment out trtpn for program to run properly;
	class sex;
	monotone reg;   			*linear regression, each variable in the VAR is imputed using all the preceeding variables in the VAR statement as covariate;
	var sex y0 y4 y8 y12 y16 y20 y24  ; 
run;

data _ontrt_temp_base; set _ontrt_eff_miout;
avisitn=0; base=y0; output;
run;
data _ontrt_temp_post; set _ontrt_eff_miout;
	avisitn=4; logaval=y4; chg=logaval-y0; output;
	avisitn=8; logaval=y8; chg=logaval-y0; output;
	avisitn=12; logaval=y12; chg=logaval-y0; output;
	avisitn=16; logaval=y16; chg=logaval-y0; output;
	avisitn=20; logaval=y20; chg=logaval-y0; output;
	avisitn=24; logaval=y24; chg=logaval-y0; output;
	run;

proc sort data=_ontrt_temp_base;
by _imputation_ usubjid;
proc sort data=_ontrt_temp_post;
by _imputation_ usubjid;
run;

	data _ontrt_temp; merge _ontrt_temp_base(keep=usubjid  _imputation_ base trtpn) 
							_ontrt_temp_post(keep=usubjid _imputation_ avisitn logaval chg sex);
	by  _imputation_ usubjid;
	run;



%macro computeMI(mcmcby=, mivar=, covar=, rtype=);

data _eff_cough_abl; set _eff_cough;
if ablfl = 'Y' then output;
else if avisitn in (4,8,12,16,20,24) then output;
run;
proc sort data= _eff_cough_abl;
  by usubjid trtpn sex avisitn;
run;
proc transpose data= _eff_cough_abl out= _eff_cough_wide (drop=_name_ _label_) prefix=y;
  by uSUBJID trtpn sex;
  var AVAL;
  id avisitn;
run;
proc sort data= _eff_cough_wide;
  by TRTPN sex;
run;

* to get monotone missing pattern;
proc mi data=_eff_cough_wide out=_eff_mono nimpute=100 seed=7264030 noprint;
    by TRTPN &mcmcby.;  		*model will be run for each combiniation of the varaibles in the BY STATEMENT;
    var  y0 y4 y8 y12 y16 y20 y24 ; 
    mcmc chain=multiple impute=monotone;
run; 

data _eff_mono1; set _eff_mono;
if y24=. then flag=1;
run;
proc sort data=_eff_mono1 out=_eff_miin;
    by  _Imputation_ trtpn;
run;

* to get imputed datasets;
proc mi data=_eff_miin out=_eff_miout nimpute=1 seed=7264030  noprint ; 
	by  _Imputation_ trtpn; 	*model will be run for each combiniation of the varaibles in the BY STATEMENT;
							    * note: sparse data for treatment by strata combination, need to comment out trtpn for program to run properly;
	class &mivar.;
	monotone reg;   			*linear regression, each variable in the VAR is imputed using all the preceeding variables in the VAR statement as covariate;
	var &mivar. y0 y4 y8 y12 y16 y20 y24  ; 
run;

data _temp_base; set _eff_miout;
avisitn=0; base=y0; output;
run;
data _temp_post; set _eff_miout;
	avisitn=4; logaval=y4; chg=logaval-y0; output;
	avisitn=8; logaval=y8; chg=logaval-y0; output;
	avisitn=12; logaval=y12; chg=logaval-y0; output;
	avisitn=16; logaval=y16; chg=logaval-y0; output;
	avisitn=20; logaval=y20; chg=logaval-y0; output;
	avisitn=24; logaval=y24; chg=logaval-y0; output;
	run;

proc sort data=_temp_base;
by _imputation_ usubjid;
proc sort data=_temp_post;
by _imputation_ usubjid;
run;

	data _temp; merge _temp_base(keep=usubjid  _imputation_ base trtpn) _temp_post(keep=usubjid _imputation_ avisitn logaval chg sex);
	by  _imputation_ usubjid;
	run;


*** note: before running the cLDA model, do the checking and create macro varaibles accordinly as done for the primary analysis;
proc mixed data=_temp ;
	where avisitn in (4 8 12 16 20 24) ;
		class usubjid trtpn avisitn sex ;
		model chg= &covar. /ddfm=kr solution ;
		by _imputation_;
		repeated avisitn/subject=usubjid type=&rtype. R rcorr;

		lsmestimate trtpn*avisitn 'Placebo, at Week 4'              [1, 1 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 4'        [1, 2 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 4'        [1, 2 1] [-1, 1 1]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 8'              [1, 1 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 8'        [1, 2 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 8'        [1, 2 2] [-1, 1 2]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 12'              [1, 1 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 12'        [1, 2 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 12'        [1, 2 3] [-1, 1 3]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 16'              [1, 1 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 16'        [1, 2 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 16'        [1, 2 4] [-1, 1 4]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 20'              [1, 1 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 20'        [1, 2 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 20'        [1, 2 5] [-1, 1 5]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 24'              [1, 1 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 24'        [1, 2 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 24'        [1, 2 6] [-1, 1 6]/exp cl;
		ods output LSMEstimates=_lsm;
run;


DATA __ot_outm; set _lsm;
where stmtno in (1,2, 4,5, 7,8, 10,11, 13,14, 16,17 );
if stmtno in (1,2) then avisitn=4;
if stmtno in (4,5) then avisitn=8;
if stmtno in (7,8) then avisitn=12;
if stmtno in (10,11) then avisitn=16;
if stmtno in (13,14) then avisitn=20;
if stmtno in (16,17) then avisitn=24;

if stmtno in (1,4,7,10,13,16) then trtpn=1;
if stmtno in (2,5,8,11,14,17) then trtpn=3;

*keep label avisitn trtpn expestimate lowerexp upperexp ;
run;

DATA __ot_outdiff; set _lsm;
where stmtno in (3, 6, 9, 12, 15, 18 );
if stmtno in (3) then do; avisitn=4; trtpn=3; end;
if stmtno in (6) then do; avisitn=8;  trtpn=3; end;
if stmtno in (9) then do; avisitn=12;  trtpn=3; end;
if stmtno in (12) then do; avisitn=16;  trtpn=3; end;
if stmtno in (15) then do; avisitn=20; trtpn=3; end;
if stmtno in (18) then do; avisitn=24; trtpn=3; end;

modelrr = (expestimate-1)*100;
modelrrL = (lowerexp-1)*100;
modelrrU = (upperexp-1)*100;
withinvar = stderr**2;
*keep label _label_ avisitn trtpn modelrr modelrrL modelrrU  ;
run;

proc sort data=__ot_outm out=__ot_outm1;
by stmtno label _imputation_;
run;

proc sort data=__ot_outdiff out=__ot_outdiff1;
by stmtno label _imputation_;
run;

ODS OUTPUT PARAMETERESTIMATES=means_rubin;
	PROC MIANALYZE DATA=__ot_outm1 ;
	BY stmtno label;
	MODELEFFECTS estimate;
	STDERR stderr;
	RUN;
	ODS OUTPUT CLOSE;
data means_rubin1; set means_rubin;
	expestimate = round(exp(estimate), .01);
	expupper = round(exp(uclmean), .01);
	explower = round(exp(lclmean), .01);
	run;

data means_rubin2;set means_rubin1;
	if stmtno in (1,4,7,10,13,16) then trtpn=1;
	else trtpn=3;
	if stmtno in (1,2) then avisitn=4;
	else if stmtno in (4,5) then avisitn=8;
	else if stmtno in (7,8) then avisitn=12;
	else if stmtno in (10,11) then avisitn=16;
	else if stmtno in (13,14) then avisitn=20;
	else if stmtno in (16,17) then avisitn=24;
	keep stmtno label avisitn trtpn expestimate expupper explower;
run;

ODS OUTPUT PARAMETERESTIMATES=diffs_rubin;
	PROC MIANALYZE DATA=__ot_outdiff1 ;
	BY stmtno label;
	MODELEFFECTS estimate;
	STDERR stderr;
	RUN;
	ODS OUTPUT CLOSE;

	data diffs_rubin1; set diffs_rubin;
	expestimate = exp(estimate);
	expupper = exp(uclmean);
	explower = exp(lclmean);
	modelrr = (expestimate-1)*100;
	modelrrL = (explower-1)*100;
	modelrrU = (expupper-1)*100;
	run;

%mend;

%computeMI(mcmcby=, mivar=, covar=%str( base trtpn avisitn trtpn*avisitn sex  ), rtype=UN);
proc sgplot data=means_rubin2;
		   title1 'model-based GMR over time - MI';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;
  
		   xaxistable expestimate / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='GMR';
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;
%computeMI(mcmcby=, mivar=%str(sex), covar=%str( base trtpn avisitn trtpn*avisitn sex  ), rtype=UN);
proc sgplot data=means_rubin2;
		   title1 'model-based GMR over time - MI regvar sex';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;

		   xaxistable expestimate / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='GMR';

		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

%computeMI(mcmcby=%str(sex), mivar=%str(sex), covar=%str( base trtpn avisitn trtpn*avisitn sex  ), rtype=UN);
proc sgplot data=means_rubin2;
		   title1 'model-based GMR over time - mcmc by sex, MI reg var sex';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;

			xaxistable expestimate / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='GMR';
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

%computeMI(mcmcby=%str(sex), mivar=%str(sex), covar=%str( base trtpn avisitn trtpn*avisitn sex  base*avisitn), rtype=UN);
proc sgplot data=means_rubin2;
		   title1 'model-based GMR over time - mcmc by sex, MI reg var sex full model';
		   scatter x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=expestimate / group=TRTPN groupdisplay=cluster clusterwidth=0.5;

			xaxistable expestimate / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='GMR';
		   yaxis  type=linear label = 'GMR'  Values=(0.3 to 1 by 0.1) grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

