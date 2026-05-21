
/*2022-4-27 global completer2 with 15 mg included in the model*/
libname ad030 "\\bardsar-prod\mk7264-cough\prot030\dataanalysis" access=readonly;


	data _ad030_fas; set ad030.adsl;
		where fasfl = 'Y';
		keep usubjid subjid trt01pn sex reg;
	run;

	proc sql;
	create table _ad030_fas_c2 as select a.*, b.eotstt1 
		from _ad030_fas as a left join ad030.adbase as b
			on a.usubjid=b.usubjid;
			
	quit;

	data _ad030_fas_c2; set _ad030_fas_c2;
	if upcase(eotstt1) = 'COMPLETED' then completer2 = 'Y';  else completer2 = 'N';
	run;
	proc freq data=_ad030_fas_c2;
	table completer2*trt01pn/norow nocol nopercent;
	run;

	data  _ad030_eff_cough; set ad030.adeff;
	where  fasfl = 'Y' and 
		  paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
		  avisitn in (0, 4,8,12,16,20,24);
	keep usubjid subjid trtp trtpn ADT fasfl paramcd aval sex  reg avisitn base ablfl chg;
	run;

	proc sql;
	create table _ad030_eff_coughc2 as select  a.usubjid, a.completer2, b.*
		from _ad030_fas_c2  as a left join _ad030_eff_cough as b
			on a.usubjid=b.usubjid;
	quit;


/*subjects in the model*/
proc sql;
create table _ad030_bign0 as select distinct trt01pn, count(distinct usubjid) as bign
	from _ad030_fas_c2
		where completer2 = 'Y' 
			group by trT01pn
				order by trt01pn;

create table _ad030_bign1 as select distinct trtpn, count(distinct usubjid) as bign
	from _ad030_eff_coughc2
		where completer2 = 'Y' and  chg ne .
			group by trtpn
				order by trtpn;
quit;

/*descriptive over time*/
proc sql;
create table _ad030_coughc2_base as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase
	from _ad030_eff_coughc2
		where completer2='Y' and ablfl = 'Y'
			group by avisitn, trtpn;

create table _ad030_coughc2_post as select distinct avisitn, trtpn, count(distinct usubjid) as N, round(exp(mean(base)),.01) as geobase, round(exp(mean(aval)),.01) as geopost, round(exp(mean(chg)),.01) as geochg
	from _ad030_eff_coughc2
		where completer2='Y' and avisitn in (4,8,12,16,20,24) and chg ne .
			group by avisitn, trtpn;
quit;
data _ad030_coughc2_all; set _ad030_coughc2_base _ad030_coughc2_post;run;

/*mixed model*/
proc mixed data=_ad030_eff_coughc2;
	where avisitn in (4 8 12 16 20 24) and completer2 = 'Y';
		class usubjid trtpn avisitn sex reg;
		model chg=base trtpn avisitn trtpn*avisitn base*avisitn sex reg/ddfm=kr solution ;
		repeated avisitn/subject=usubjid type=un;

		lsmestimate trtpn*avisitn 'Placebo, at Week 4'              [1, 1 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 4'        [1, 3 1]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 4'        [1, 3 1] [-1, 1 1]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 8'              [1, 1 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 8'        [1, 3 2]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 8'        [1, 3 2] [-1, 1 2]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 12'              [1, 1 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 12'        [1, 3 3]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 12'        [1, 3 3] [-1, 1 3]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 16'              [1, 1 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 16'        [1, 3 4]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 16'        [1, 3 4] [-1, 1 4]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 20'              [1, 1 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 20'        [1, 3 5]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 20'        [1, 3 5] [-1, 1 5]/exp cl;

		lsmestimate trtpn*avisitn 'Placebo, at Week 24'              [1, 1 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg, at Week 24'        [1, 3 6]/exp cl;
		lsmestimate trtpn*avisitn 'MK-7264 45 mg vs. Placebo, at Week 24'        [1, 3 6] [-1, 1 6]/exp cl;
		ods output LSMEstimates=_ad030_lsmc2;
run;


data  _ad030_lsmc2_1; set _ad030_lsmc2;
where stmtno not in (3,6,9,12,15,18);
model_gmr = round(expestimate, .01);
model_l = round(lowerexp, .01);
model_u = round(upperexp, .01);
keep stmtno label model_gmr model_l model_u;
run;

data  _ad030_lsmc2_2; set _ad030_lsmc2;
where stmtno in (3,6,9,12,15,18) ;
rr = round((expestimate-1)*100, .01);
rrl = round((lowerexp-1)*100, .01);
rru = round((upperexp-1)*100, .01);
keep  stmtno label rr rrl rru;
run;



















proc sgplot data=_eff_n2_cough;
where completer2 = 'N'
		   title1 'Summary of cough freq ratio(post/baseline) by treatment COMPLETER 1';
		   scatter x=avisitn y=chg / group=usubjid groupdisplay=cluster clusterwidth=0.5markerattrs=(size=7 symbol=circlefilled) ;
		   series x=avisitn y=chg / group=usubjid groupdisplay=cluster clusterwidth=0.5;

			*xaxistable n / x=avisitn class=trtpn colorgroup=trtpn valueattrs=(weight=bold) title='N' ;
  
		   yaxis  type=linear label = 'logchg'  /*Values=(0.3 to 1 by 0.1)*/ grid;
		   xaxis label = 'Week' values=(0 to 24 by 4) grid;

		run;

