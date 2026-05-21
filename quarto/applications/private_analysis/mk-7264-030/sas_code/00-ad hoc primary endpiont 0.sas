/****************************
AD HOC ANALYSIS -0 descriptive summaries
2022-4-12

Author: XInru Ren


*************************************/

%LET DIR = \\bardsar-prod\mk7264-cough\prot030-china-ext-main\adhoc\outputs;
%PUT &DIR;


data _eff_cough; set lptda.adeff;
where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
	  paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (0, 4,8,12,16,20,24);
keep usubjid subjid trtp trtpn paramcd aval sex  avisitn base ablfl chg pchg ;
run;

	
	proc sql;
	create table _eff_cough1 as select a.*, b.chro, b.chron, b.ccdurg1, b.ccdurg1n, b.age 
		from _eff_cough as a left join lptda.adbase as b on a.usubjid = b.usubjid;
	quit;

	data _eff_cough1; set _eff_cough1;
	length agegr1 $20;
	if age<60 then agegr1n = 1; else agegr1n=2;
	if agegr1n=1 then agegr1 = '<60 yrs'; else agegr1 = '>=60 yrs';
	run;


/**********************************************
COUGH FREQ BY TRT and 
gender 
chro(UCC vs RCC)
ccdurg1 (<10, >=10)
agegr1n (<60, >=60)
******************************************************/

** by gender;
proc sql;
create table _aa0 as select distinct trtpn, avisitn, sex, 
		count(distinct usubjid) as N, 
		round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
	from _eff_cough
		where ablfl = 'Y'
			group by trtpn, sex, avisitn
				order by trtpn, sex, avisitn;

create table _aa00 as select distinct trtpn, avisitn, sex, 
		count(distinct usubjid) as N, round(exp(mean(chg)), .01) as gmr, 
		round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
	from _eff_cough
		WHERE AVISITN IN (4,8,12,16,20,24) AND CHG NE .
			group by trtpn, sex, avisitn
				order by trtpn, sex, avisitn;
quit;

DATA _AA; SET _AA0 _AA00; if avisitn=0 then gmr = 1; RUN;
PROC SORT DATA=_AA; BY TRTPN SEX AVISITN; RUN;

ods rtf file="&DIR.\adhoc0plot0ot0trt0demo.rtf" ;
		proc sgplot data=_AA;
		   title1 'Summary of cough freq ratio(post/baseline) by treatment in Male';
		   where sex = 'M';
		   scatter x=avisitn y=gmr / group=trtpn groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=trtpn groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;


		proc sgplot data=_AA;
		   title1 'Summary of cough freq ratio(post/baseline) by treatment in Female';
		   where sex = 'F';
		   scatter x=avisitn y=gmr / group=trtpn groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=trtpn groupdisplay=cluster clusterwidth=0.5;

			xaxistable n / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) 
								title='N' ;
			xaxistable gm_base / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
   


		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;


		proc sgplot data=_AA;
		   title1 'Summary of cough freq ratio(post/baseline) in placebo by GENDER';
		   where trtpn=1;
		   scatter x=avisitn y=gmr / group=sex groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=sex groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=SEX colorgroup=SEX valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=sex colorgroup=sex valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;


		proc sgplot data=_AA;
		   title1 'Summary of cough freq ratio(post/baseline) in 45 mg by GENDER';
		   where trtpn=3;
		   scatter x=avisitn y=gmr / group=sex groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=sex groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=SEX colorgroup=SEX valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=sex colorgroup=sex valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;

/*
data _aa1m _aa1f _aa3m _aa3f; set _aa;
if trtpn = 1 and sex = 'M' then output _aa1m;
else if trtpn=1 and sex='F' then output _aa1f;
else if trtpn=3 and sex = 'M' then output _aa3m;
else output _aa3f;
*gmr_=round(gmr, .01);
run;

data _aa1; merge 
_aa1m(rename=(gmr = male_gmr   n=male_n   gm_base = male_base   gm_aval = male_aval) 	keep= trtpn avisitn gmr n gm_base gm_aval) 
_aa1f(rename=(gmr = female_gmr n=female_n gm_base = female_base gm_aval = female_aval) keep=avisitn gmr n gm_base gm_aval);
by avisitn;

run;

data _aa3; merge 
_aa3m(rename=(gmr = male_gmr n=male_n   gm_base = male_base   gm_aval = male_aval) 	keep= trtpn avisitn gmr n gm_base gm_aval) 
_aa3f(rename=(gmr=female_gmr n=female_n gm_base = female_base gm_aval = female_aval) keep=avisitn gmr n gm_base gm_aval);
by avisitn;
run;

** treatment side-by-side;
data _aa4; merge 
	_aa1(in=a rename=(male_gmr = plb_m_r female_gmr = plb_f_r 
					  male_n = plb_m_n  female_n = plb_f_n
					  male_base = plb_m_b female_base =plb_f_b
  					  male_aval = plb_m_p female_aval = plb_f_p)
 		) 
	_aa3(in=c  rename=(male_gmr = mk45_m_r female_gmr = mk45_f_r 
						male_n = mk45_m_n  female_n = mk45_f_n
						male_base = mk45_m_b female_base =mk45_f_b
  					  	male_aval = mk45_m_p female_aval = mk45_f_p)
		);

by avisitn;
drop trtpn;
run;
*/


** by CHRO (ucc vs ucc);
	proc sql;
	create table _bb0 as select distinct trtpn, avisitn, chro, 
			count(distinct usubjid) as N, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			WHERE ABLFL = 'Y'
				group by trtpn, chro, avisitn
					order by trtpn, chro, avisitn;

	create table _bb00 as select distinct trtpn, avisitn, chro, 
			count(distinct usubjid) as N, round(exp(mean(chg)), .01) as gmr, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			where avisitn in (4,8,12,16,20,24) and chg ne .
				group by trtpn, chro, avisitn
					order by trtpn, chro, avisitn;
	quit;

	data _bb; set _bb0 _bb00; IF AVISITN=0 THEN GMR=1;run;
	proc sort data=_bb; by trtpn chro avisitn;run;

proc sgplot data=_BB;
		   title1 'Summary of cough freq ratio(post/baseline) in placebo by primary UCC/RCC';
		   where trtpn=1;
		   scatter x=avisitn y=gmr / group=chro groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=chro groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=chro colorgroup=chro valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=chro colorgroup=chro valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;


		proc sgplot data=_BB;
		   title1 'Summary of cough freq ratio(post/baseline) in 45 mg by UCC/RCC';
		   where trtpn=3;
		   scatter x=avisitn y=gmr / group=chro groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=chro groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=chro colorgroup=chro valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=chro colorgroup=chro valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;

	**ccdurg1;
	proc sql;
	create table _cc0 as select distinct trtpn, avisitn, ccdurg1, 
			count(distinct usubjid) as N, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			WHERE ABLFL = 'Y'
				group by trtpn, ccdurg1, avisitn
					order by trtpn, ccdurg1, avisitn;

	create table _cc00 as select distinct trtpn, avisitn, ccdurg1, 
			count(distinct usubjid) as N, round(exp(mean(chg)), .01) as gmr, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			where avisitn in (4,8,12,16,20,24) and chg ne .
				group by trtpn, ccdurg1, avisitn
					order by trtpn, ccdurg1, avisitn;
	quit;

	data _cc; set _cc0 _cc00; IF AVISITN=0 THEN GMR=1;run;
	proc sort data=_cc; by trtpn ccdurg1 avisitn;run;
proc sgplot data=_CC;
		   title1 'Summary of cough freq ratio(post/baseline) in placebo by duration of chronic cough';
		   where trtpn=1;
		   scatter x=avisitn y=gmr / group=ccdurg1 groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=ccdurg1 groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=ccdurg1 colorgroup=ccdurg1 valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=ccdurg1 colorgroup=ccdurg1 valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;


		proc sgplot data=_CC;
		   title1 'Summary of cough freq ratio(post/baseline) in 45 mg by duration of chronic cough';
		   where trtpn=3;
		   scatter x=avisitn y=gmr / group=ccdurg1 groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=ccdurg1 groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=ccdurg1 colorgroup=ccdurg1 valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=ccdurg1 colorgroup=ccdurg1 valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;

	**agegr1n;
	proc sql;
	create table _dd0 as select distinct trtpn, avisitn, AGEGR1 ,
			count(distinct usubjid) as N, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			WHERE ABLFL = 'Y'
				group by trtpn, AGEGR1, avisitn
					order by trtpn, AGEGR1, avisitn;

	create table _dd00 as select distinct trtpn, avisitn, AGEGR1, 
			count(distinct usubjid) as N, round(exp(mean(chg)), .01) as gmr, 
			round(exp(mean(base)), .01) as gm_base, round(exp(mean(aval)), .01) as gm_aval
		from _eff_cough1 
			where avisitn in (4,8,12,16,20,24) and chg ne .
				group by trtpn, AGEGR1, avisitn
					order by trtpn, AGEGR1, avisitn;
	quit;

	data _dd; set _dd0 _dd00; IF AVISITN=0 THEN GMR=1;run;
	proc sort data=_dd; by trtpn AGEGR1 avisitn;run;

proc sgplot data=_dd;
		   title1 'Summary of cough freq ratio(post/baseline) in placebo by age group';
		   where trtpn = 1;
		   scatter x=avisitn y=gmr / group=agegr1 groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=agegr1 groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=agegr1 colorgroup=agegr1 valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=agegr1 colorgroup=agegr1 valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;


		proc sgplot data=_dd;
		   title1 'Summary of cough freq ratio(post/baseline) in 45 mg by age group';
		   where trtpn = 3;
		   scatter x=avisitn y=gmr / group=agegr1 groupdisplay=cluster clusterwidth=0.5
								      markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=gmr / group=agegr1 groupdisplay=cluster clusterwidth=0.5;
		   xaxistable n / x=avisitn class=agegr1 colorgroup=agegr1 valueattrs=(weight=bold) 
								title='N' ;
		   xaxistable gm_base / x=avisitn class=agegr1 colorgroup=agegr1 valueattrs=(weight=bold) 
								title='baseline GM of subjects at each post-baseline timepoint' ;
		   yaxis  type=linear label = 'Geometric Mean (SD) '  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;


		run;
		title1;

ods rtf close;


/*DESCRIPTIVE SUMMARY
0. 24-hour coughs frequency (GM) over time by treatment 
0. GM folder rise (ratio of 24-hour coughs frequency at post-baseline/baseline) over time by treatment 
*/

* FREQUENCY OVER TIME (GM);

data _eff_totfeq; set lptda.adeff;
where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
	  paramcd = 'TOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (0, 4,8,12,16,20,24);
keep usubjid subjid trtpn trtp sex  paramcd avisitn base ablfl aval chg pchg ;
run;



proc sql;

create tale _zz0 as select distinct trtpn, trtp, avisitn,
		exp(mean(aval)) as freq_gm
	from _eff_cough where ablfl = 'Y'
		group by trtpn, trtp, avisitn
			order by avisitn,trtpn, trtp;

create table _zz1 as select distinct trtpn, trtp, avisitn,  
		exp(mean(aval)) as freq_gm
	from _eff_cough where avisitn in (4,8,12,16,20,24)
		group by trtpn, trtp, avisitn
			order by avisitn, trtpn, trtp;
quit;

data _zz2; set _zz0 _zz1;run;
proc sort data=_zz2; by avisitn trtpn trtp; run;

proc transpose data=_zz2 out=_zz2t prefix = freq_gm;
by avisitn;
id trtpn;
var freq_gm;
run;

data _zz3; set _zz2t;
freq_gm1_ = round(freq_gm1, .01);
freq_gm3_ = round(freq_gm3, .01);
label freq_gm1_ = 'Placebo cough frequency GM'
freq_gm3_ = 'MK 45mg cough frequency GM'
avisitn = 'Week';
keep avisitn freq_gm1_ freq_gm3_;
run;

** GM FOLD RISE;
PROC SQL;
create table _yy1 as select distinct trtpn, trtp, avisitn,  
		count(distinct usubjid) as N, exp(mean(base)) as freqbase_gm, exp(mean(CHG)) as ratio_gm
	from _eff_cough where avisitn in (4,8,12,16,20,24) and chg ne .
		group by trtpn, trtp, avisitn
			order by avisitn, trtpn, trtp;
quit;

data _yy2; set _yy1;
freqbase_gm = round(freqbase_gm, .01);
ratio_gm = round(ratio_gm, .01);
label freqbase_gm = 'baseline cough frequency GM'
ratio_gm = 'post-baseline/baseline ratio GM'
avisitn = 'Week';
keep trtp avisitn N freqbase_gm ratio_gm;
run;

