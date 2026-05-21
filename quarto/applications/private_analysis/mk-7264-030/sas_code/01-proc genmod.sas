data _eff_cough; set lptda.adeff;
where country = 'CHN' and fascnfl = 'Y' and trt01pn in (1,3) and 
	  paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (4,8,12,16,20,24);
keep usubjid subjid trt01pn trtpn paramcd aval sex basetype dtype avisit avisitn base ablfl chg pchg ;
run;
proc sql;
*  records available at baseline and each post baseline evaluated;
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

proc mixed data=_eff_cough;
where avisitn in (4 8 12 16 20 24);
	class usubjid trtpn avisitn sex /*reg*/;
	model chg=base trtpn avisitn trtpn*avisitn  sex /*reg*//ddfm=kr solution;
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
	ods output LSMEstimates=__ot_lsm;
run;

DATA __ot_outm; set __ot_lsm;
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

DATA __ot_outdiff; set __ot_lsm;
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

data _lda_24hcough_ot; set __part1;run;



data six;
input case city$ @@;
do i=1 to 4;
input age smoke wheeze @@;
output;
end;
datalines;
1 portage 9 0 1 10 0 1 11 0 1 12 0 0
2 kingston 9 1 1 10 2 1 11 2 0 12 2 0
3 kingston 9 0 1 10 0 0 11 1 0 12 1 0
4 portage 9 0 0 10 0 1 11 0 1 12 1 0
5 kingston 9 0 0 10 1 0 11 1 0 12 1 0
6 portage 9 0 0 10 1 0 11 1 0 12 1 0
7 kingston 9 1 0 10 1 0 11 0 0 12 0 0
8 portage 9 1 0 10 1 0 11 1 0 12 2 0
9 portage 9 2 1 10 2 0 11 1 0 12 1 0
10 kingston 9 0 0 10 0 0 11 0 0 12 1 0
11 kingston 9 1 1 10 0 0 11 0 1 12 0 1
12 portage 9 1 0 10 0 0 11 0 0 12 0 0
13 kingston 9 1 0 10 0 1 11 1 1 12 1 1
14 portage 9 1 0 10 2 0 11 1 0 12 2 1
15 kingston 9 1 0 10 1 0 11 1 0 12 2 1
16 portage 9 1 1 10 1 1 11 2 0 12 1 0
;


/*** missing observation;*/
data six_missobs; set six;
if case = 1  and i=1 then delete;
if case = 4  and i=2 then delete;
if case = 7  and i=3 then delete;
if case = 13 and i=4 then delete;
run;
/*** missing outcome;*/
data six_missout; set six;
if case = 1  and i=1 then wheeze=.;
if case = 4  and i=2 then wheeze=.;
if case = 7  and i=3 then wheeze=.;
if case = 13 and i=4 then wheeze=.;
run;


proc genmod data=six_missout;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  type=ar(1) covb corrw;
run;

proc genmod data=six_missout;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ar(1)covb corrw;
run;
*same thing;




* to test whether order of timepoints matters ;
proc genmod data=six_missobs;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  type=ar(1) covb corrw;
run;

proc genmod data=six_missobs;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ar(1)covb corrw;
run;
* different;
* withinsubject matters;


proc genmod data=six_missobs;   /*missing obs*/
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ar(1) covb corrw;
run;

proc genmod data=six_missout;  /*missing outcome*/
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ar(1) covb corrw;
run;
* different;
**consulted with Xuechan and decided to add WITHINSUBJECT= argument;


** TEST IF ANY IMPUTATION WAS DONE;
proc genmod data=six_missobs;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ind covb corrw;
run;

proc genmod data=six_missout;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age type=ind covb corrw;
run;
* different;



* missing;
data _missout4; set six;
if  case=16 and i=4 then wheeze=.;
run;
proc genmod data=_missout4;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case /  withinsubject=age  type=ar(1) covb corrw;
run;

data _missobs4; set six;
where case<=15  or  (case = 16 and  i ne 4);
run;
proc genmod data=_missobs4;
class case city age;
model wheeze(event='1') = city age smoke / dist=bin;
repeated subject=case / withinsubject=age   type=ar(1) covb corrw;
run;

data sex; set _eff_cough;
where avisitn =24;
keep usubjid trtpn sex;
run;
proc freq data=sex;
table sex*trtpn;
run;
