/**************************************************************
* step 0: identify completers defined as having week 24 data;
****************************************************************/

data _eff_cough; set lptda.adeff;
where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
	  paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (0, 4,8,12,16,20,24);
keep usubjid subjid trtp trtpn ADT fascnfl paramcd aval sex  avisitn base ablfl chg;
run;

	proc sort data=_eff_cough out=adeff_Ltotfeq;
	by trtpn trtp usubjid;
	run;
	proc transpose data=adeff_Ltotfeq out=adeff_Ltotfeqt(drop=_name_ _label_) prefix = v;
	by trtpn trtp usubjid fascnfl;
	id avisitn;
	var chg;
	run;   * wide format per subject per treatment;
	proc sort data=adeff_Ltotfeqt out=adeff_Ltotfeqtt;  * wide format per subject per treatment;
	by trtpn v0 v4 v8 v12 v16 v20 v24;
	run;


    

proc sql;
create table _eff_coughc as select * from _eff_cough 
	where usubjid in select distinct usubjid from _eff_cough where avisitn=24;
quit;

/**************************************************************
* step 1: big model on completer and normality test;
****************************************************************/
	proc mixed data=_eff_coughc;
	where avisitn in (4 8 12 16 20 24);
		class usubjid trtpn avisitn sex /*reg*/;
		model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex /ddfm=kr solution OUTPM=normres RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;

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
		ods output LSMEstimates=_lsmc;
run;

data _lsmc1; set _lsmc;
where stmtno in (3,6,9,12,15,18);
rr = round((expestimate-1)*100, .01);
rrl = round((lowerexp-1)*100, .01);
rru = round((upperexp-1)*100, .01);
run;


proc mixed data=_eff_coughc;
	where avisitn in (4 8 12 16 20 24);
		class usubjid trtpn avisitn sex /*reg*/;
		model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex /ddfm=kr solution OUTPM=normres RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;
		run;
	PROC UNIVARIATE DATA = normres NORMAL NOPRINT;
					VAR scaledresid;
					OUTPUT OUT = normtest PROBN = testnorm;
	RUN;


	ods select histogram qqplot;
	ods output TestsForNormality = work.normal_test;
	ods output BasicMeasures = work.measures;
		proc univariate data=_eff_coughc normal;
		WHERE AVISITN IN (4,8,12,16,20,24);
		class /*trtpn*/ trtp AVISITN;
		var chg;
		histogram chg / NORMAL ;
		qqplot chg /Normal(mu=est sigma=est color=red l=1);
		run;

	proc sort data=measures; by avisitn  trtp ;run;
	proc transpose data=measures out=measures_loc;
	by avisitn trtp;
	id locmeasure;
	var locvalue;
	run;
	proc transpose data=measures out=measures_var;
	by avisitn  trtp;
	id varmeasure;
	var varvalue;
	run;
	data measures_all; merge measures_loc measures_var;
	by avisitn  trtp;
	AVISITN_=INPUT(strip(AVISITN), BEST.);
	hi = mean+std_deviation;
	lo = mean-std_deviation;
	drop _name_ _label_ mode;
	run;

	proc sort data=measures_all; by trtp avisitn_ ;run;
	
	data normal_test_pval; set normal_test;
	where test = 'Shapiro-Wilk';
	if pvalue <0.001 then flag=1;
	AVISITN_=INPUT(strip(AVISITN), BEST.);
	run;

/**************************************************************
* step 2: baseline char;
****************************************************************/

proc sql;
create table _eff_comp0 as select 
	a.*, b.tr01sdt, b.tr01edt, b.age, b.sex, b.eotstt1, b.dctreas1, b.TRCMPS,
	b.BCSEVAS, B.BCSEVASN, BVASSCR, B.SMOKE, B.CHRO, b.CCDURG1, BCO24HR, BCO24SCR
		from adeff_Ltotfeqt as a left join lptda.adbase as b
			on a.usubjid =b.usubjid;
quit;

DATA _eff_comp0; SET _eff_comp0;
length agegr1 $10;
IF V24 NE . THEN COMPLETER='Y'; ELSE COMPLETER = 'N';
if age<60 then agegr1 = '<60'; else agegr1 = '>=60';
RUN;
	libname adhoc '\\bardsar-prod\mk7264-cough\prot030-china-ext-main\adhoc\datasets';
	data adhoc.fascn_completer; set _eff_comp0;
	keep usubjid fascnfl completer;
	run;
 
    
PROC FREQ DATA=_eff_comp0;
table TRTPN/MISSING nopercent norow;
TABLE COMPLETER*TRTPN/MISSING nopercent norow;
RUN;

PROC FREQ DATA=_eff_comp0 NOPRINT;
TABLE completer*SEX*TRTPN /OUT=T_SEX;
TABLE completer*agegr1*TRTPN/out=t_age;
TABLE completer*chro*TRTPN/out=t_chro;
TABLE completer*CCDURG1*TRTPN/out=t_dur;
TABLE completer*BCO24HR*TRTPN/out=t_24hc;

TABLE completer*BCSEVAS*TRTPN/out=t_vas;
TABLE completer*SMOKE*TRTPN/out=t_smoke;
table completer*eotstt1*trtpn/out=t_eot;
table completer*eotstt1*dctreas1*trtpn/out=t_dctreas1;
run;

** treatment discontinuation and compliance;
proc transpose data=t_eot out=t_eott;
by completer eotstt1 ;
id trtpn;
var count;
run;

proc transpose data=t_dctreas1 out=t_dctreas1t;
*where eotstt1 = 'Discontinued';
by completer eotstt1 dctreas1;
id trtpn;
var count;
run;

proc means data=_eff_comp0;
class completer trtpn;
var TRCMPS;
run;

proc univariate data=_eff_comp0;
		class completer trtpn;
		var TRCMPS;
		histogram TRCMPS ;
		run;



/*MVN BASED ON SCALED RESIDUAL FROM MODEL*/
	proc mixed data=_eff_coughc;
	where avisitn in (4 8 12 16 20 24);
		class usubjid trtpn avisitn sex;
		model chg= trtpn avisitn trtpn*avisitn  /ddfm=kr solution OUTPM=__normres_C RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;
		
run;
	PROC UNIVARIATE DATA = __normres_C NORMAL NOPRINT;
					VAR scaledresid;
					OUTPUT OUT = __normtest_C PROBN = testnorm;
	RUN;


	proc mixed data=_eff_cough;
	where avisitn in (4 8 12 16 20 24);
		class usubjid trtpn avisitn sex;
		model chg= trtpn avisitn trtpn*avisitn  /ddfm=kr solution OUTPM=__normres_A RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;
run;

	PROC UNIVARIATE DATA = __normres_A NORMAL NOPRINT;
					VAR scaledresid;
					OUTPUT OUT = __normtest_A PROBN = testnorm;
	RUN;


/*****/

	DATA ADEX; SET LPTDA.ADEX;
	WHERE PRDFLAG  = 'MAIN STUDY';
	FORMAT DATE YYMMDD10.;
	DO I = ASTDT TO AENDT;
		DATE = I;
		OUTPUT;
	END;
	KEEP USUBJID DATE EXNUMDOS ETCD prdflag;
	RUN;

	PROC SQL;  /*total daily dose*/
		CREATE TABLE ADEX1 AS SELECT DISTINCT usubjid, prdflag, eTCD, date, sum(exnumdos) as tdd
			from adex group by usubjid,  etcd, date; 

	QUIT;

PROC SQL;  /*per subject per date of exposure with avisitn date as landmark*/
	CREATE TABLE _EFF_COMP1 AS SELECT A.USUBJID, A.TR01SDT, A.TR01EDT, A.TRTPN, A.COMPLETER,  B.*
		FROM _EFF_COMP0 AS A LEFT JOIN ADEX1 AS B
			ON A.USUBJID = B.USUBJID;
	CREATE TABLE _EFF_COMP2 AS SELECT A.* , B.AVISITN, B.ADT
		FROM _EFF_COMP1 AS A LEFT JOIN _EFF_COUGH AS B
			ON A.USUBJID=B.USUBJID AND A.DATE = B.ADT;
QUIT;

proc sort data=_eff_comp2 out=_eff_comp2_mk;
where trtpn=3;
by usubjid date avisitn adt;
run;

data _eff_comp3_mk; set _eff_comp2_mk;
by usubjid date avisitn adt;
retain v adt_;
if avisitn ne . then do;v=avisitn; adt_=adt; output;end;
else output;
format adt_ yymmdd10.;
run;

proc sql;
create table _eff_comp4_mk as select distinct usubjid, trtpn, completer, tr01sdt, tr01edt, v, adt_, count(usubjid) as ontrtday
	from _eff_comp3_mk
		group by usubjid, trtpn, completer, tr01sdt, tr01edt, v, adt_
			order by completer, usubjid, v;
quit;
	proc transpose data=_eff_comp4_mk  out=_eff_comp4_mkt;
	by  completer usubjid;
	id v;
	var ontrtday;
	run;
	data _eff_comp4_mkt1(drop=_name_); 
		retain completer usubjid _0 _4 _8 _12 _16 _20 _24; 
		set _eff_comp4_mkt;
		
		label 	_0 = '[base, week 4)'
				_4 = '[week 4, 8)'
				_8 = '[week 8, 12)'
				_12 = '[week 12, 16)'
				_16 = '[week 16, 20)'
				_20 = '[week 20, 24)'
				_24 = '[week 24, -';
	run;

%LET DIR = \\bardsar-prod\mk7264-cough\prot030-china-ext-main\adhoc\outputs;
%PUT &DIR;
	PROC EXPORT DATA= WORK._eff_comp4_mkt1 
			            OUTFILE= "&dir.\adhoc0mk450exposure0by0completer.xls" 
			            DBMS=xls LABEL REPLACE;
			     NEWFILE=YES;
			RUN;
proc sql;
create table _eff_comp5_ex as select distinct completer, trtpn, v, 
	count(distinct usubjid) as n, round(mean(ontrtday), .01) as mean, round(std(ontrtday), .01) as sd, median(ontrtday)as median
		from _eff_comp4_mk
			group by completer, trtpn, v;
quit;



/*2022-4-22*******************************
completer 1: with 24 week cough data available
completer 2: EOTSTT1 = 'COMPLETED'
*******************************************/

	** STEP 422_1: COMPLETER 1 VS COMPLETER 2 
	HOW DO THE TWO POPULATION OVERLAP? 
	DATA POINTS EXCLUDED?
	;

DATA SUBJ_comp12; SET _eff_comp0;
if upcase(eotstt1) = 'COMPLETED' THEN COMPLETER2 = 'Y'; ELSE COMPLETER2 = 'N';
rename completer=completer1;
RUN;

proc freq data=SUBJ_comp12;
table completer1*trtpn
	completer2*trtpn
	completer1*completer2
	/ missing nocol nopercent norow;
	run;


proc sql;
create table _eff_coughc12 as select a.*, b.completer1, b.completer2
	from _eff_cough as a left join subj_comp12 as b
		on a.usubjid = b.usubjid;
quit;

ods select histogram qqplot;
	ods output TestsForNormality = work.normal_test;
	ods output BasicMeasures = work.measures;
		proc univariate data=_eff_coughc12 normal;
		WHERE AVISITN IN (4,8,12,16,20,24) and completer2='Y';
		class /*trtpn*/ trtp AVISITN;
		var chg;
		*histogram chg / NORMAL ;
		*qqplot chg /Normal(mu=est sigma=est color=red l=1);
		run;
	
	data normal_test_pval; set normal_test;
	where test = 'Shapiro-Wilk';
	if pvalue <0.001 then flag=1;
	AVISITN_=INPUT(strip(AVISITN), BEST.);
	run;

** data points excluded;
proc sql;
create table __coughc1_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase
	from _eff_coughc12
		where completer1='Y' and ablfl = 'Y'
			group by avisitn, trtpn;

create table __coughc1_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase, round(exp(mean(aval)),.01) as geopost, round(exp(mean(chg)),.01) as geochg
	from _eff_coughc12
		where completer1='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
			group by avisitn, trtpn;


create table __coughc2_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase
	from _eff_coughc12
		where completer2='Y' and ablfl = 'Y'
			group by avisitn, trtpn;

create table __coughc2_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase, round(exp(mean(aval)),.01) as geopost, round(exp(mean(chg)),.01) as geochg
	from _eff_coughc12
		where completer2='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
			group by avisitn, trtpn;
quit;

data __coughc1_all; set __coughc1_base __coughc1_post;run;
data __coughc2_all; set __coughc2_base __coughc2_post;run;

	proc sgplot data=__coughc1_all;
		   title1 'Summary of cough freq ratio(post/baseline) by treatment COMPLETER 1';
		   scatter x=avisitn y=geochg / group=trtpn groupdisplay=cluster clusterwidth=0.5markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=geochg / group=trtpn groupdisplay=cluster clusterwidth=0.5;

			xaxistable n / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='N' ;
			xaxistable geobase / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='baseline GM of subjects at each post-baseline timepoint' ;
  
		   yaxis  type=linear label = 'Geometric Mean'  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;
	proc sgplot data=__coughc2_all;
		   title1 'Summary of cough freq ratio(post/baseline) by treatment COMPLETER 2 ';
		   scatter x=avisitn y=geochg / group=trtpn groupdisplay=cluster clusterwidth=0.5 markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=geochg / group=trtpn groupdisplay=cluster clusterwidth=0.5;

			xaxistable n / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='N' ;
			xaxistable geobase / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='baseline GM of subjects at each post-baseline timepoint' ;

		   yaxis  type=linear label = 'Geometric Mean'  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;
		title1;

	proc mixed data=_eff_coughc12;
	where avisitn in (4 8 12 16 20 24) and completer1='Y';
		class usubjid trtpn avisitn sex /*reg*/;
		model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex /ddfm=kr solution OUTPM=normres RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;

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
		ods output LSMEstimates=_lsmc1;
run;
		data _lsmc1_; set _lsmc1;
		where stmtno in (3,6,9,12,15,18);
		rr = round((expestimate-1)*100, .01);
		rrl = round((lowerexp-1)*100, .01);
		rru = round((upperexp-1)*100, .01);
		run;

	proc mixed data=_eff_coughc12;
	where avisitn in (4 8 12 16 20 24) and completer2='Y';
		class usubjid trtpn avisitn sex /*reg*/;
		model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex /ddfm=kr solution OUTPM=normres RESIDUAL VCIRY;
		repeated avisitn/subject=usubjid type=un;

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
		ods output LSMEstimates=_lsmc2;
run;

		data _lsmc2_; set _lsmc2;
		where stmtno in (3,6,9,12,15,18);
		rr = round((expestimate-1)*100, .01);
		rrl = round((lowerexp-1)*100, .01);
		rru = round((upperexp-1)*100, .01);
		run;

 *awake cough;

	data _eff_awc; set lptda.adeff;
	where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
		  paramcd = 'LAWFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
		  avisitn in (0, 4,8,12,16,20,24);
	keep usubjid subjid trtp trtpn ADT fascnfl paramcd aval sex  avisitn base ablfl chg;
	run;
	proc sql;
	create table _eff_awcc12 as select a.*, b.completer1, b.completer2
		from _eff_awc as a left join subj_comp12 as b
			on a.usubjid = b.usubjid;
	quit;

	** data points excluded;
	proc sql;
	create table __awcc1_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase
		from _eff_awcc12
			where completer1='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table __awcc1_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase, round(exp(mean(aval)),.01) as geopost, round(exp(mean(chg)),.01) as geochg
		from _eff_awcc12
			where completer1='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;


	create table __awcc2_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase
		from _eff_awcc12
			where completer2='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table __awcc2_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase, round(exp(mean(aval)),.01) as geopost, round(exp(mean(chg)),.01) as geochg
		from _eff_awcc12
			where completer2='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;
	quit;

	data __awcc1_all; set __awcc1_base __awcc1_post;run;
	data __awcc2_all; set __awcc2_base __awcc2_post;run;

%macro lda(in=, out=, cond=);
		proc mixed data=&in.;
		where avisitn in (4 8 12 16 20 24) and &cond.='Y';
			class usubjid trtpn avisitn sex /;
			model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex /ddfm=kr solution ;
			repeated avisitn/subject=usubjid type=un;

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
			ods output LSMEstimates=&out.;
	run;
			data &out._; set &out.;
			where stmtno in (3,6,9,12,15,18);
			rr = round((expestimate-1)*100, .01);
			rrl = round((lowerexp-1)*100, .01);
			rru = round((upperexp-1)*100, .01);
			run;
%mend;
%lda(in=_eff_awcc12, out=_lsmawc1, cond =completer1 );
%lda(in=_eff_awcc12, out=_lsmawc2, cond =completer2);

****** LCQ;


	data _eff_lcq; set lptda.adqs;
		where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
			paramcd = 'LCQTOT'  and ANL01FL = 'Y' and dtype ='' and 
			avisitn in (0, 4,8,12,16,20,24);
		keep usubjid subjid trtpn paramcd aval base chg pchg sex avisitn ablfl CSDRES1 CSDRES2 LCQRES VASRES PGICRES ;
		run;

	proc sql;
	create table _eff_lcqc12 as select a.*, b.completer1, b.completer2
		from _eff_lcq as a left join subj_comp12 as b
			on a.usubjid = b.usubjid;
	quit;

	** data points excluded;

	proc sql;
	create table __lcq_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from _eff_lcqc12
			where  ablfl = 'Y'
				group by avisitn, trtpn;

	create table __lcq_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from _eff_lcqc12
			where  avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;

	create table __lcqc1_base as select distinct avisitn, trtpn, count(distinct usubjid) as N,  round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from _eff_lcqc12
			where completer1='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table __lcqc1_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from _eff_lcqc12
			where completer1='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;


	create table __lcqc2_base as select distinct avisitn, trtpn, count(distinct usubjid) as N,  round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from _eff_lcqc12
			where completer2='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table __lcqc2_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from _eff_lcqc12
			where completer2='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;
	quit;

			data __lcq_all; set __lcq_base __lcq_post;run;
	data __lcqc1_all; set __lcqc1_base __lcqc1_post;run;
	data __lcqc2_all; set __lcqc2_base __lcqc2_post;run;

%macro logit(in=, resp=, out1=, out2=, cond=);
			PROC SQL;
			* BIGN: subject with data available at week 24;
			CREATE TABLE _24COUNT AS SELECT DISTINCT TRTPN, COUNT(DISTINCT USUBJID) as bign
				FROM &in. 
					WHERE AVISITN=24 AND chg NE . and &cond.='Y'
						GROUP BY TRTPN;
			*smalln: responder at week 24 as observed;					
			create table _24resp as select distinct trtpn, count(distinct usubjid)  as smalln
				from &in.
					where avisitn=24 and &resp. =1 and &cond.='Y'
						group by trtpn;

			QUIT;

			data _part1;
				merge _24count _24resp;
				by trtpn;
				rate = round(smalln/bign*100, .1);
				run;


			PROC SQL;

			* Nstar: subject in analysis;
			CREATE TABLE _24anal AS SELECT DISTINCT TRTPN, COUNT(DISTINCT USUBJID) as nstar
				FROM &in.
					WHERE  CHG NE . and &cond.='Y'
						GROUP BY TRTPN;

			QUIT;

			proc glimmix data=&in. maxopt=50 pconv=1E-5; 
			where  &cond.='Y' and avisitn in (4,8,12,16,20,24);
			    class subjid TRTPN avisitn sex ;
			    model &resp. (event='1')= TRTPN avisitn TRTPN*avisitn base base*avisitn sex / noint dist=binary link=logit ddfm=sat; 	

			    lsmestimate trtpn*avisitn "Placebo, at Week 24"                         [1, 1 6] / exp cl;
			    lsmestimate trtpn*avisitn "MK-7264 45 mg, at Week 24"                  [1, 2 6] / exp cl;
			    lsmestimate trtpn*avisitn "MK-7264 45 mg vs. Placebo, at Week 24"       [1, 2 6] [-1, 1 6] / exp cl;

			    ods output LSMEstimates=_lsmestimate;

				random avisitn / sub=subjid type=un residual; 
			run;

			data _outp; set _lsmestimate;
				where statement in (1,2);
				_p_est = expestimate/(1+expestimate);
				 p_est  = round(_p_est*100, .1);
			    

				if statement=1 then trtpn=1;
				if statement=2 then trtpn=3;
				run;

			data _outPdiff; 
				merge _outp(where=(trtpn=1) rename=( _p_est = _p_est1  )) 
					  _outp(where=(trtpn=3) rename=( _p_est = _p_est3  ));

				_pdiff_est = _p_est3 - _p_est1;
				 pdiff_est = round(_pdiff_est*100, .01);

				trtpn=3;
				run;

			data _outORdiff; set _lsmestimate;
				where index(label, 'vs. Placebo')>0;

				or = round(expestimate, .01);
				orlower = round(explower, .01);
				orupper = round(expupper, .01);
				pvalue = round(probt, 0.001);

				trtpn=3;
				run;


			/*MN CI for est p difference*/

			data lsm2; set _outp;
			p_est = _p_est;
			run;
			data _outmn1; set _outp;
				varp = (_p_est *(1- _p_est))**2 * stderr**2;
				nadj = int((_p_est *(1- _p_est))/varp) - 1;   * adjusted n;
				cadj = int(nadj* _p_est);				* adjusted responder count;
			run;

			data _outmn1_1; set _outmn1(where=(trtpn=3));
				n1 = nadj;
				s1 = cadj; 

			run;

			data _outmn1_2; set _outmn1(where=(trtpn=1));
				n2 = nadj;
				s2 = cadj; 
			run;
				
			data _outmn2; merge _outmn1_1 _outmn1_2;
				stratum = 1;
				examparm = '1';
			keep stratum examparm n1 s1 n2 s2;
			run;

			%include '\\bardsar-test\mk7264-cough\prot030-china-ext-main\documents\stats\xinru\rate0compare.sas';

			%rate0compare(
				 	/*Macro Setup Parameters*/
				  input_dataset 	= _outmn2
				 ,output_dataset 	= __outmn_observed
				 ,output_dataset2	= __outmn_adjusted
				 ,stratum_var 		= stratum
				 ,no_treatment 		= 2
				 ,by_var 			= examparm
				 ,debug 			= N
				 	/*Statistical Analysis Parameters*/
				 ,method 			= MN 	/*Choose: MN or BW */
				 ,continuity_corr	= 0 	/*Choose: 0 (None), 1 (Mehrotra and Railkar), or 2 (Mantel-Haenszel) */
				 ,delta 			= 0 	/* -1 <= DELTA < +1 */
				 ,weight_schm 		= 3 	/*Choose: 1 (Equal), 2 (Sample Size), 3 (CMH) , 4 (MR) , 5 (Pooled) 6 (Userspecified)*/
				 ,ci 				= 95 	/*Confidence Level */
				 ,exposure_adj 		= 0 	/*Choose: 0=No exposure adjustment, 1= Exposure Adjustment*/
				 ,testing 			= 2 	/*Choose: 0=No P-values, 1=1 sided P-value, 2=2 sided P-value*/
				); 

				data _outmnci; set __outmn_adjusted;
					trtpn=3;
					L_MNCI = round(l_bound, .01);
					u_MNCI = round(u_bound, .01);
					keep trtpn u_bound l_bound L_MNCI u_MNCI;
					run;

			data _part2; 
				merge _24anal(in=a) _outp(in=b) _outpdiff(in=c) _outmnci(in=d) _outORdiff(in=ce);;
				by trtpn;
				keep trtpn nstar p_est pdiff_est or L_MNCI u_MNCI orlower orupper pvalue;
				run;

			data &out1.; set _part1;run;
			data &out2.; set _part2;run;
%mend;
%logit(in=_eff_lcqc12, resp=LCQRES, out1=__p1_lcqc1, out2=__p2_lcqc1, cond=completer1);
%logit(in=_eff_lcqc12, resp=LCQRES, out1=__p1_lcqc2, out2=__p2_lcqc2, cond=completer2);


*** >=30% reduction in 24hr cough;

data _eff_cough_org; set lptda.adeff;
	where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
		paramcd = 'TOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
		avisitn in (4,8,12,16,20,24);
	keep usubjid subjid trtpn paramcd aval base chg pchg sex avisitn TOTRES30 TOTRES50 TOTRES70;
	run;

	proc sql;
	create table _eff_cough_orgc12 as select a.*, b.completer1, b.completer2
		from _eff_cough_org as a left join subj_comp12 as b
			on a.usubjid = b.usubjid;
	quit;

%logit(in=_eff_cough_orgc12, resp=TOTRES30, out1=__p1_totres30c1, out2=__p2_totres30c1, cond=completer1);
%logit(in=_eff_cough_orgc12, resp=TOTRES30, out1=__p1_totres30c2, out2=__p2_totres30c2, cond=completer2);

/** ************************
********VAS
*****************************/

	data _eff_VAS; set lptda.adqs;
		where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
			paramcd = 'VASWK'  and ANL01FL = 'Y' and dtype ='' and 
			avisitn in (0, 4,8,12,16,20,24);
		keep usubjid subjid trtpn paramcd aval base chg pchg sex avisitn ablfl CSDRES1 CSDRES2 LCQRES VASRES PGICRES ;
		run;

	proc sql;
	create table _eff_VASc12 as select a.*, b.completer1, b.completer2
		from _eff_VAS as a left join subj_comp12 as b
			on a.usubjid = b.usubjid;
	quit;

	** data points excluded;
%macro desc(in=, out=, );
	proc sql;
	create table &out._base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from &in.
			where  ablfl = 'Y'
				group by avisitn, trtpn;

	create table &out._post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from &in.
			where  avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;

	create table &out.c1_base as select distinct avisitn, trtpn, count(distinct usubjid) as N,  round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from &in.
			where completer1='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table &out.c1_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from &in.
			where completer1='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;


	create table &out.c2_base as select distinct avisitn, trtpn, count(distinct usubjid) as N,  round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase
		from &in.
			where completer2='Y' and ablfl = 'Y'
				group by avisitn, trtpn;

	create table &out.c2_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(mean(base),.1) as meanbase, round(std(base),.01) as sdbase, round(mean(aval),.1) as meanpost, round(std(aval),.01) as sdpost, round(mean(chg),.1) as meanchg, round(std(chg),.01) as sdchg 
		from &in.
			where completer2='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
				group by avisitn, trtpn;
	quit;

%mend;
%desc(in=_eff_VASc12, out=__VAS );

			data __VAS_all; set __VAS_base __VAS_post;run;
	data __VASc1_all; set __VASc1_base __VASc1_post;run;
	data __VASc2_all; set __VASc2_base __VASc2_post;run;

%logit(in=_eff_VASc12, resp=VASRES, out1=__p1_VASRESc1, out2=__p2_VASRESc1, cond=completer1);
%logit(in=_eff_VASc12, resp=VASRES, out1=__p1_VASRESc2, out2=__p2_VASRESc2, cond=completer2);

/************************
*** CSD
*/
	data _eff_csd; set lptda.adqs;
		where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
			paramcd = 'CSDTOTWK'  and ANL01FL = 'Y' and dtype ='' and 
			avisitn in (0, 4,8,12,16,20,24);
		keep usubjid subjid trtpn paramcd aval base chg pchg sex avisitn ablfl CSDRES1 CSDRES2 LCQRES VASRES PGICRES ;
		run;

	proc sql;
	create table _eff_CSDc12 as select a.*, b.completer1, b.completer2
		from _eff_CSD as a left join subj_comp12 as b
			on a.usubjid = b.usubjid;
	quit;

%desc(in=_eff_CSDc12, out=__CSD );
	data __CSD_all; set __CSD_base __CSD_post;run;
	data __CSDc1_all; set __CSDc1_base __CSDc1_post;run;
	data __CSDc2_all; set __CSDc2_base __CSDc2_post;run;
%logit(in=_eff_CSDc12, resp=CSDRES1, out1=__p1_CSDRES1c1, out2=__p2_CSDRES1c1, cond=completer1);
%logit(in=_eff_CSDc12, resp=CSDRES1, out1=__p1_CSDRES1c2, out2=__p2_CSDRES1c2, cond=completer2);
%logit(in=_eff_CSDc12, resp=CSDRES2, out1=__p1_CSDRES2c1, out2=__p2_CSDRES2c1, cond=completer1);
%logit(in=_eff_CSDc12, resp=CSDRES2, out1=__p1_CSDRES2c2, out2=__p2_CSDRES2c2, cond=completer2);
