libname cough 'C:\Users\leolarr2\OneDrive - Merck Sharp & Dohme, Corp\documents\Myfiles\Projects\Subgroup-Identification\MRCTs_Consistency\mk7264-030\sas_dataset';

* Merge with baseline adbase;

data baseline; set cough.adbase; run;
proc contents; run;

data basecovs; set baseline;
keep usubjid subjid AGE BCO24HR BCO24SCR BCSEVAS BMIBL BMIGR1 BVASSCR CCDUR CCDURG1 CHRO ETHNIC HARQTOT HEIGHTBL 
RACE REG SEX SMKSTAT SMOKE WEIGHTBL;
run;

proc contents data=cough.adeff; run;


data eff_cough; set cough.adeff;
where trt01pn in (1,3) and 
	  FASFL = 'Y' and paramcd = 'LTOTFEQ' and basetype = 'Pretreatment' and ANL01FL = 'Y' and dtype ='' and 
	  avisitn in (0,4,8,12,16,20,24);
keep usubjid subjid trtp trtpn paramcd aval sex age avisitn base ablfl chg pchg ADT ADTM ADY ontrtfl TR01EDT TR01EDTM TR01SDT TR01SDTM TR02EDT TR02EDTM TR02SDT TR02SDTM AWTARGET TRTEDT TRTEDTM TRTSDT TRTSDTM;
run;

proc freq data = eff_cough(where=(avisitn=24)); table trtp; run;
proc univariate data=eff_cough(where=(avisitn=24)); var aval base chg pchg; run;

*proc contents data=eff_cough; run;


data eff_cough_offtreat; set eff_cough(where=(ontrtfl='N')); run;

*proc print data=eff_cough_offtreat; run;


data eff_cough; set eff_cough;
if adt>=trtsdt then do;
analysis_day=(adt-trtsdt+1);
end;
else if adt<trtsdt then do;
analysis_day=(adt-trtsdt);
end;
diff_day=ady-analysis_day;
run;

proc univariate data=eff_cough; var diff_day; run;

*proc print data=eff_cough(where=(diff_day ne 0)); run;

data eff_cough; set eff_cough;
if trtedt>=trtsdt then do;
offtrt_day=(trtedt-trtsdt+1);
end;
run;

*proc print data=eff_cough(where=(ady>offtrt_day)); run;
*proc contents; run;

* Set outcomes with analysis days corresponding to ontrtfl='N' to missing;

proc sort data=eff_cough; by usubjid avisitn; run;

data eff_cough; set eff_cough; by usubjid avisitn;
if ontrtfl='N' then do; 
chg_ontrt=.; pchg_ontrt=.;
end;
else if ontrtfl='Y' then do;
chg_ontrt=chg; pchg_ontrt=pchg;
end;
run;

*proc print data=eff_cough(where=(ady>offtrt_day)); run;

data eff_week24; set eff_cough; by usubjid avisitn;
if last.usubjid then output;
run;

proc freq data=eff_week24; table trtp avisitn trtp*avisitn; run;


/* From China cohort explorations.doc
Planned Treatment (N)	WK4	WK8	WK12 WK16	WK20	WK24
1	Placebo	            65	61	59	 49	    43	41	51
3	MK-7264 45 mg BID	64	56	55	 45	49	39	46
*/

data eff_cough_out; set eff_cough(keep=usubjid subjid age sex ady avisitn trtp trtpn 
aval base chg pchg ontrtfl offtrt_day chg_ontrt pchg_ontrt); run;

* Merge with baseline confounders;

data eff_baseline;
merge eff_cough_out(in=a) basecovs(in=b); by usubjid;
if a=1;
run;


proc export data=eff_baseline
outfile="C:\Users\leolarr2\OneDrive - Merck Sharp & Dohme, Corp\documents\Myfiles\Projects\Subgroup-Identification\MRCTs_Consistency\mk7264-030\outdata\global_eff_cough.csv"
dbms=csv
replace;
run;



















