libname bpd "C:\Users\brandie\Documents\Childrens Hospital\BPD\Gerber microbiome\Data\Raw Data";
%let dataset= bpd.merged2015_07_07; 

Proc format;
	value b  2='Mild'
			3='Moderate'
			4='Severe';
run;
/*select samples closest to day 7 but only within 5-9 days*/

data rf1;
	set bpd.merged2015_07_07;
	select=abs(days_from_birth-7);
run; 
proc means data=rf1;
	where inclin=1 and seq < 2;
	class study_id;
	var select;
	output out=select min=min_day7;
run;
proc sort data=rf1;
	by study_id;
run;
data rf2;
	merge rf1 select (where=(_TYPE_=1));
	by study_id;
	if select=min_day7 and select <= 2 then cross_sec_day7=1;
	/*some day 7 samples had missing day 7 dates, still want to select them*/
	if study_id = '0124S' and day= 7 then cross_sec_day7=1;
	/*tied*/
	if study_id ='0206U' and day=0 then cross_sec_day7=.;
run;
data bpd.merged2015_07_07_2;
	set rf2;
	where inclin=1 and cross_sec_day7=1 and bpd1 in (2 3 4);
run;

/**********************/
/*  Table 1           */
/*demographics for ALL*/
/**********************/

proc sort data=&dataset.;
	by study_id;
run;
data dem1;
	set &dataset.;
	where seq=2 and bpd1 in (2 3 4);
	by study_id;
	if first.study_id then table=1;
run;
proc freq data=dem1;
	where table=1;
	table (gender preeclampsia prematrupmembran chorioamnionitis corticosteroids 
	mom_ethncty_1 mom_ethncty_2 cesarean treatmentsurfact multgest pneumonia sga)*bpd1/chisq;
run;
proc means data=dem1 n mean std median min max;
	where table=1;
	class bpd1;
	var weight zscore gestational_age mom_age days_from_birth daysmv;
run;
proc npar1way data=dem1 anova wilcoxon;
	where table=1;
	class bpd1;
	var weight zscore gestational_age mom_age days_from_birth daysmv;
run;

/****************************************************/
/****************************************************/
/****************************************************/
/****************************************************/
/****************************************************/
/****************************************************/
/****************************************************/
/****************************************************/

/*cross-sectional analysis*/
/*select samples closest to day 7 but only within 5-9 days*/

data rf1;
	set &dataset.;
	select=abs(days_from_birth-7);
run; 
proc means data=rf1;
	where seq < 2;
	class study_id;
	var select;
	output out=select min=min_day7;
run;
proc sort data=rf1;
	by study_id;
run;
data rf2;
	merge rf1 select (where=(_TYPE_=1));
	by study_id;
	if select=min_day7 and select <= 2 then cross_sec_day7=1;
	/*some day 7 samples had missing day 7 dates, still want to select them*/
	if study_id = '0124S' and day= 7 then cross_sec_day7=1;
	/*tied*/
	if study_id ='0206U' and day=0 then cross_sec_day7=.;
run;
data cross;
	set rf2;
	where inclin=1 and cross_sec_day7=1 and bpd1 in (2 3 4);
run;
/*demographics for cross-sectional*/
proc freq data=cross;
	where seq=2;
	table (gender preeclampsia prematrupmembran chorioamnionitis corticosteroids 
	mom_ethncty_1 mom_ethncty_2 cesarean treatmentsurfact multgest pneumonia sga)*bpd1/chisq;
run;
proc means data=cross n mean std median min max;
	where seq=2;
	class bpd1;
	var weight zscore gestational_age days_from_birth daysmv;
run;
proc npar1way data=cross anova wilcoxon;
	where seq=2;
	class bpd1;
	var weight zscore gestational_age mom_age days_from_birth daysmv;
run;
/*Total bacterial load*/
proc means data=cross n nmiss median min max;
	where seq < 2;
	var TBL;
run;
proc npar1way data=cross anova wilcoxon;
	where seq < 2;
	class bpd1;
	var TBL;
run;
/*diversity*/
proc npar1way data=cross anova wilcoxon;
	where seq < 2;
	class bpd1;
	var shannonH_mean shannonE_mean;
run;
/***********/
/* Figure 2*/
/***********/
ods html style=journal;
proc sgplot data=cross noautolegend ;
	where seq < 2;
	vbox shannonH_mean/category=bpd1 group=bpd1;
	format bpd1 b.;
	yaxis label="Shannon Diversity Index" labelattrs=(size=16) valueattrs=(size=14) max=2.5;
	xaxis label= "BPD Status" labelattrs=(size=16) valueattrs=(size=14);
run;

proc sgplot data=cross noautolegend;
	where seq < 2;
	vbox ShannonE_mean/category=bpd1 group=bpd1;
	format bpd1 b.;
	yaxis label="Evenness" labelattrs=(size=16) valueattrs=(size=14) max=.6;
	xaxis label= "BPD Status" labelattrs=(size=16) valueattrs=(size=14);
run;


/***********/
/* Figure 3*/
/***********/
proc means data=cross n max;
	class taxa;
	var perc_total;
	output out=maxpt max = maxpt;
run;
proc print data=maxpt;run;
proc sort data=cross;
	by taxa;
run;
data Gerber;
	merge cross maxpt (where=(_TYPE_=1));
	by taxa;
run;
title1 'BPD Status';
proc sgpanel data=Gerber;
	where maxpt > 10 ;
	panelby bpd1 /novarname  columns=3 rows=1 uniscale=row;
	vbar library/response=perc_total group=taxa groupdisplay=stack fill;
	rowaxis max=100 label="Relative Abundance";
	colaxis display=(NOVALUES);
	attrib library label=' ' ;
	format bpd1 b.;
run;


/*compare staph and ureaplasma across groups*/

proc npar1way data=cross wilcoxon;
	where taxa="Staphylococcus";
	class bpd1;
	var perc_total;
run;
proc npar1way data=cross wilcoxon;
	where taxa="Ureaplasma";
	class bpd1;
	var perc_total;
run;

/*******************************************/
/*******************************************/
/*										   */
/*										   */
/*										   */
/*change in bacterial communities over time*/
/*										   */
/*										   */
/*										   */
/*******************************************/
/*******************************************/

proc sgplot data=&dataset.;
	where seq=2 and bpd1 in ( 1 2 3 4 );
	vbox tbl/ category=day group=bpd1 grouporder=ascending;
run;
proc freq data=&dataset.;
	where seq=2 and bpd1  in (2 3 4);
	table study_id*bpd1/out=n;
run;
data longitudinal;
	merge &dataset. (where=( bpd1  in (2 3 4))) n;
	by study_id;
	if count < 2 then delete; 
	ltbl=log10(tbl);
run;
/*select one obs per subject to generate table 1*/
proc sort data=longitudinal;
	by study_id day;
run;
data dem;
	set longitudinal;
	p=lag1(study_id);
	do i=0 to 1;
		if p=' ' then visit=1;

		else if study_id = p then visit=v+i;
		else if study_id ne p then visit=1;
		v=lag1(visit);
	end;
	drop v i p;
run;

proc freq data=dem;
	where visit=1;
	table (gender preeclampsia prematrupmembran chorioamnionitis corticosteroids 
	mom_ethncty_1 mom_ethncty_2 cesarean treatmentsurfact multgest pneumonia sga)*bpd1/chisq;
run;
proc means data=dem n mean std median min max;
	where visit=1;
	class bpd1;
	var weight zscore gestational_age days_from_birth daysmv;
run;
proc npar1way data=dem anova wilcoxon;
	where visit=1;
	class bpd1;
	var weight zscore gestational_age days_from_birth daysmv;
run;

proc genmod data=longitudinal;
	class study_id bpd1/param=ref ref=first;
	model tbl=days_from_birth days_from_birth*days_from_birth days_from_birth*days_from_birth*days_from_birth/link=log;
	repeated subject=study_id;
	  assess var=(days_from_birth) / resample
                       seed=603708000;
run;	
proc genmod data=longitudinal;
	class study_id bpd1/param=ref ref=first;
	model tbl=bpd1 bpd1*days_from_birth /link=log;
	repeated subject=study_id;
run;
ods graphics on;
proc mixed data=longitudinal plots=all;
	where seq=2 and bpd1 in (2 3 4 );
	class study_id bpd1;
	model lTBL= bpd1 bpd1*days_from_birth bpd1*days_from_birth*days_from_birth/noint solution;
	random intercept days_from_birth/subject=study_id group=bpd1;
run;
/*just look at variance components*/
proc mixed data=longitudinal plots=all;
	where seq=2 and bpd1 in (2 3 4 );
	class study_id bpd1 ;
	model lTBL= bpd1/noint solution;
	random intercept /subject=study_id group=bpd1;
	repeated / subject=study_id group=bpd1;
run;
proc sgplot data=longitudinal;
	where seq=2 and bpd1 in (2 3 4 );
	series y=tbl x=days_from_birth/ group=study_id;
run;
proc sgplot data=longitudinal;
	where seq=2 and bpd1 in (2 3 4 );
	vbox tbl/ category=day group=bpd1 grouporder=ascending;
run;
proc genmod data=longitudinal;
	where seq=2 and bpd1 in (2 3 4 );
	class study_id bpd1;
	model TBL= bpd1 bpd1*days_from_birth/noint dist=normal link=log;
	repeated subject=study_id;
	estimate 'slope sev vs mild' bpd1*days_from_birth -1 0 1;
run;
proc genmod data=longitudinal;
	where seq=2 and bpd1 in (2 3 4 );
	class study_id ;
	model TBL= days_from_birth/ dist=normal link=log;
	repeated subject=study_id;
run;
/*look at all samples*/
	proc sgpanel data=longitudinal;
			panelby bpd1 ccenter /novarname uniscale=row;
			vbar library/response=perc_total group=taxa groupdisplay=stack fill;
			rowaxis max=100 label="Relative Abundance";
			attrib days_from_birth label=' ';
		run;
	proc sgpanel data=longitudinal;
			panelby bpd1 study_id /novarname uniscale=row;
			vbar days_from_birth/response=perc_total group=taxa groupdisplay=stack fill;
			rowaxis max=100 label="Relative Abundance";
			attrib days_from_birth label=' ';
		run;
/*controls*/
proc sgpanel data=&dataset.;
	where bpd1=1 and taxa in ('Corynebacterium' 'Enterobacteriaceae' 'Enterococcus' 'Escherichia-Shigella' 'Staphylococcus');
	panelby study_id /novarname noheader columns=3 rows=1 uniscale=row;
	vbar days_from_birth/response=perc_total group=taxa groupdisplay=stack fill;
	rowaxis max=100 label="Relative Abundance";
	attrib days_from_birth label='Day';
run;
proc means data=&dataset.;
	where bpd1=1;
	class taxa;
	var perc_total;
run;
		title "Decomposition of Diversity";
		/*create multivariate dataset of seq counts*/
		proc means data=longitudinal noprint;
			var seq;
			output out=numseq max=numseq;
		run;
		data _NULL_;
			set numseq;
			call symput('norg', numseq);
		run;
		%put &norg;
		%macro data ;
			%do i=1 %to &norg;
			data out&i.;
				set longitudinal;
				where seq=&i;
				seq&i.=seq_count;
				keep library study_id days_from_birth bpd1 seq&i.;
			run;
			proc sort data=out&i.;
				by library;
			run;
			%end;
		%mend;

		%data();
		data mult2;
			merge out1 - out438;
			by library;
		run;

		/*IML program that calculates alpha, beta and gamma, from Marcon 2012 paper*/
		%macro beta_decomp (data, id, out);
				proc iml ;
					start data;
						use &data.;
						read all var _all_ into seq;
						read all var _CHAR_ into type;
						read all var {&id.} into &id.;
						Nlibrary=nrow(seq);
						dim=ncol(seq);
						/*fixed error in code below from 3 to 4*/
						/*only changed the 4th decimal place */
						count0=seq [ ,4:dim];
						check=count0[+, ];
						sub=loc(check>0);
						count=count0[ , sub];
						Norg=ncol(count);
						c= Count[+,+];
						c_i=Count[ ,+];
						c_j=Count[+, ];	
						*print Nlibrary Norg;
					finish;

					Start gamma;
						/*set 0 counts to 1 prior to log transform*/
						cg=c_j/c;
						g=-cg#log( choose(cg>0, cg, .) );
						H_gamma=g[,+];
						gamma_hill=exp(H_gamma);
						*print H_gamma gamma_hill;
					Finish;
					Start alpha;
						ra=Count/c_i;
						a=-ra#log(choose(ra>0, ra, .));
						h_alpha_i=a[, +];
						alpha_i_hill=exp(h_alpha_i);
						h_alpha=t(h_alpha_i)*(c_i/c);
						alpha_hill =exp(h_alpha);
					Finish;
					start Beta;
						b0=(Count/c_i)#(c/c_j);
						ra=Count/c_i;
						b=ra#log(choose(b0>0, b0, . ));
						h_beta_i=b[ ,+];
						beta_i_hill=exp(h_beta_i);
						h_beta=t(h_beta_i)*(c_i/c);
						beta_hill =exp(h_beta);
						*print h_beta beta_hill;
						create &out. var{&id. Nlibrary c h_alpha alpha_hill h_beta beta_hill H_gamma gamma_hill};
						append; close &out.;
					finish;
				run data; run gamma; run alpha; run beta; 
		%mend;
		
		data mult_type;
			set mult2;
			type_id=compress(study_id);
		run;
		proc freq data=mult_type noprint;
			table type_id/out=id1;
		run;
		proc freq data=id1 noprint;
			table type_id/out=id2 outcum;
		run;
		proc means data=id2 noprint;
			var CUM_FREQ;
			output out=npp max=g;
		run;
		data _NULL_;
			set npp;
			call symput('g', g);
		run;
		%put &g;
		proc sort data=mult_type;
			by type_id;
		run;
		data mult_type2;
			merge  id2(keep=type_id CUM_FREQ  rename=(CUM_FREQ=group)) mult_type;
			by type_id;
		run;
		%macro data_decomp ;
			%do i=1 %to &g.;
				data g_&i.;
					set mult_type2;
					where group=&i;
					/*needed to drop 1 numeric var to make code run*/
					drop type_id ;
				run;
				%beta_decomp (g_&i., group, out_&i.);
			%end;
		%mend;

		options spool;
		%data_decomp();
		data mult_type3;
			merge  mult_type2 (keep=group library bpd1 study_id days_from_birth type_id) out_1 - out_94;
			by group;
		run;
		title "Beta diversity for divergence within subject";
		proc print data=mult_type3;
			where  H_BETA ne .;
			var type_id bpd1 study_id NLIBRARY H_Beta beta_hill H_gamma gamma_hill;
		run;
		proc means data=mult_type3 n median min max maxdec=2;
			where  H_BETA ne .;
			var beta_hill H_gamma gamma_hill;
		run;
		data bpd.mult_subj3;
			set mult_type3 (keep= bpd1 study_id NLIBRARY H_Beta beta_hill H_gamma gamma_hill);
			where  H_Beta ne .;
			NT = (BETA_HILL-1)/(NLIBRARY -1);
		run;
		proc means data=bpd.mult_subj3 n median min max maxdec=2;
			var NLIBRARY NT beta_hill H_gamma gamma_hill;
		run;
		proc print data=bpd.mult_subj3;
			var study_id bpd1  NLIBRARY NT beta_hill H_gamma gamma_hill;
		run;
		title ' ';
		proc sgplot data=bpd.mult_subj3 noautolegend;
			vbox NT/category=bpd1 group=bpd1;
			format bpd1 b.;
			xaxis label= ' ';
		run;
	
		proc npar1way data=bpd.mult_subj3 wilcoxon;
			class bpd1;
			var NT;
			format NT b.;
		run;
		proc npar1way data=bpd.mult_subj3 wilcoxon;
			where bpd1 in (2 4 );
			class bpd1;
			var NT;
			format NT b.;
		run;
		data nt;
			set bpd.mult_subj3;
			if bpd1 in (3 4) then BPD2=3;
			if bpd1 = 2 then BPD2=2;
		run;
		proc npar1way data=nt wilcoxon;
			class BPD2;
			var NT;
		run;
		proc means data=bpd.mult_subj3 n median p25 p75 maxdec=2;
			class bpd1;
			var NT;
		run;		

		proc sort data=bpd.mult_subj3;
			by study_id;
		run;
		proc sort data=bpd.merged2015_07_07;
			by study_id;
		run;
		data bpd.mergedwbeta2015_07_07;
			merge bpd.mult_subj3(in=subj) longitudinal (where=(bpd1 in (2 3 4)));
			by study_id;
			if subj=1 then subj1=1;
		run;

		proc print data=bpd.mergedwbeta2015_07_07;
			where NLIBRARY >1 and seq=2;
			*var bpd1 study_id nt NLIBRARY H_Beta beta_hill  perc_total days_from_birth day;
		run;
		
proc means data=bpd.mergedwbeta2015_07_07 n max;
	where NLIBRARY >1 ;
	class taxa;
	var perc_total;
	output out=maxpt max = maxpt;
run;

proc sort data=bpd.mergedwbeta2015_07_07;
	by taxa;
run;
data Gerber;
	merge bpd.mergedwbeta2015_07_07 maxpt (where=(_TYPE_=1));
	by taxa;
	nt2=round(nt,0.01);
	if nt2 > 0 then do;
		id=compress(nt2)||"_"||compress(study_id);
	end;
	else if nt2 < 0.01 then do;
			id="0.00_"||compress(study_id);
	end;
	if nt >= 0.25 then highbeta=1;
	if nt ne . and nt < 0.25 then highbeta=0;
run;
proc means data=gerber n median min max;
	where NLIBRARY > 1 and seq=2;
	class bpd1;
	var nt;
run;

proc freq data=Gerber noprint;
	where NLIBRARY > 1 and seq=2 and subj1=1;
	table highbeta*bpd1*study_id/out=id;
run;
proc freq data=id;
	table highbeta*bpd1/chisq trend;
run;

proc sort data=gerber;
	by bpd1;
run;

options nodate nonumber;
data attrs;
input id $ value $25. fillcolor $;
datalines;
taxa	Acinetobacter 			Blue
taxa	Bacillus 				Black 
taxa	Corynebacterium			Green 
taxa	Enterobacteriaceae		Olive 
taxa	Enterococcus  			tan
taxa	Escherichia-Shigella  	Purple
taxa	Fusobacterium  			rose
taxa	Gammaproteobacteria 	Grey
taxa	Haemophilus 			stgb 
taxa	Klebsiella 				Red 
taxa	Neisseria 				Maroon 
taxa	Prevotella 				violet 
taxa	Proteus 				lime
taxa	Pseudomonadales 		Aqua 
taxa	Pseudomonas 			Pink
taxa	Serratia 				Brown 
taxa	Sneathia 				salmon
taxa	Staphylococcus 			Steel
taxa	Stenotrophomonas  		Yellow
taxa	Streptococcus 			Orange
taxa	Ureaplasma 				pap
;;;;
run;
proc print data=attrs;
run;
/*************/
/* Figure 5  */
/*************/
		proc sgpanel data=Gerber dattrmap=attrs;
			by bpd1;
			where maxpt > 50 and NLIBRARY >1;
			panelby id /novarname rows=5 columns=8 skipemptycells;
			vbar day/response=perc_total group=taxa groupdisplay=stack fill attrid=taxa;
			rowaxis max=100 label="Relative Abundance";
			attrib days_from_birth label=' ';
		run;

proc sgplot data=bpd.mergedwbeta2015_07_07;
	where taxa='Staphylococcus';
	pbspline y=perc_total x=days_from_birth/group=bpd1 smooth=50000000;
run;
proc sgplot data=bpd.mergedwbeta2015_07_07;
	where taxa='Staphylococcus';
	reg y=perc_total x=days_from_birth/group=bpd1 degree=2;
run;
/*indicates that the Beta-Binomial might be a beter fit - its U-shaped*/
proc capability data=bpd.mergedwbeta2015_07_07;
		where taxa='Staphylococcus';
		var perc_total;
		comphistogram perc_total/class=bpd1 nrows=3;
run;
data model;
	set bpd.mergedwbeta2015_07_07;
		where taxa='Staphylococcus' and subj1 = 1;
	y=seq_count;
	days=days_from_birth;
	ltotal=log(total);
	n= total;
	r=seq_count;
	day10=max(0, days -10);
	sev=(BPD1=4);
	mod=(bpd1=3);
	lRA=log(perc_total+1);
	steroids=(corticosteroids=1);
	keep study_id days y r n ltotal day10 shannonH_mean bpd1 sev mod lra zscore daysmv steroids;
run;

proc sort data=model;
	by study_id days;
run;
proc genmod data=model;
	class study_id bpd1;
	model y= bpd1 bpd1*days bpd1*day10 daysmv steroids/ dist=nb offset=ltotal type3;
	repeated subject=study_id;
run;
proc genmod data=model;
	class study_id bpd1;
	model y= bpd1 zscore/ dist=nb offset=ltotal type3;
	repeated subject=study_id;
run;
proc mixed data=model;
	class study_id bpd1;
	model lra=bpd1 bpd1*days/solution;
	random intercept/ subject=study_id group=bpd1;
	repeated/ subject=study_id group=bpd1;
run;
	/*Beta-binomial - random intercept*/
	proc nlmixed data=model tech=quanew;*tech=newrap;
		parms b0=-1 b1=.6 b2=1 b3=0.2 b4=-0.15 b5=-.2 rho=.6 s2u=.2;* c0=1 c1=.8 c2=1.5;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + ui;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		*rho=exp(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random ui ~ NORMAL(0, s2u*s2u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 21days' b0 + 21*b3;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5;
		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5, b1 + 21*b4;

	run;
	/*join point model*/
		proc nlmixed data=model tech=quanew;*tech=newrap;
		parms b0=-1 b1=.9 b2=1.3 b3=0.2 b4=-0.15 b5=-.2 b6=-.3 b7=.2 b8=.2  s1u=.2 c0=.6 c1=0 c2=0 ;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev + u1;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		rho=(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random u1 ~ NORMAL(0, s1u*s1u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		contrast 'slope10 by group' b7, b8;
		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 10 days' b0 + 10*b3;
		estimate 'mod 10 days' b0 + b1 + 10*b3 + 10*b4;
		estimate 'sev 10 days' b0 + b2 + 10*b3 + 10*b5;
		estimate 'mild 14days' b0 + 14*b3 +4*b6;
		estimate 'mod 14 days' b0 + b1 + 14*b3 + 14*b4 + 4*b6 + 4*b7;
		estimate 'sev 14days' b0 + b2 + 14*b3 + 14*b5 + 4*b6 + 4*b8;
		estimate 'mild 21days' b0 + 21*b3 +11*b6;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4 + 11*b6 + 11*b7;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5 + 11*b6 + 11*b8;

		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 10 days sev vs mild' b2 + 10*b5;

		estimate 'diff at 14 days sev vs mild' b2 + 14*b5 + 4*b8;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5 + 11*b8;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5 + 11*b8, b1 + 21*b4 + 11*b7;
		contrast 'diff at 14 days sev, mod vs mild' b2 + 14*b5 + 4*b8, b1 + 14*b4 + 4*b7;
		contrast 'diff at 10 days sev, mod vs mild' b2 + 10*b5, b1 + 10*b4;
		contrast 'diff at 7 days sev, mod vs mild' b2 + 7*b5, b1 + 7*b4;

	*predict p out=pred;
	predict exp(b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev)/ (1+ exp(b0 + b1*mod+ b2*sev 
			+ b3*days + b4*days*mod + b5*days*sev + b6*day10 + b7*day10*mod + b8*day10*sev)) out=pred;
	run;
	/*sensitivity analysis with clinical covariates*/
	proc nlmixed data=model tech=quanew;*tech=newrap;
		parms b0=-1 b1=.9 b2=1.3 b3=0.2 b4=-0.15 b5=-.2 b6=-.3 b7=.2 b8=.2  s1u=.2 c0=.6 c1=0 c2=0 b9=.01 b10=0.1;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev + b9*daysmv + b10*steroids + u1;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		rho=(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random u1 ~ NORMAL(0, s1u*s1u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		contrast 'slope10 by group' b7, b8;
		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 10 days' b0 + 10*b3;
		estimate 'mod 10 days' b0 + b1 + 10*b3 + 10*b4;
		estimate 'sev 10 days' b0 + b2 + 10*b3 + 10*b5;
		estimate 'mild 14days' b0 + 14*b3 +4*b6;
		estimate 'mod 14 days' b0 + b1 + 14*b3 + 14*b4 + 4*b6 + 4*b7;
		estimate 'sev 14days' b0 + b2 + 14*b3 + 14*b5 + 4*b6 + 4*b8;
		estimate 'mild 21days' b0 + 21*b3 +11*b6;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4 + 11*b6 + 11*b7;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5 + 11*b6 + 11*b8;

		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 10 days sev vs mild' b2 + 10*b5;

		estimate 'diff at 14 days sev vs mild' b2 + 14*b5 + 4*b8;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5 + 11*b8;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5 + 11*b8, b1 + 21*b4 + 11*b7;
		contrast 'diff at 14 days sev, mod vs mild' b2 + 14*b5 + 4*b8, b1 + 14*b4 + 4*b7;
		contrast 'diff at 10 days sev, mod vs mild' b2 + 10*b5, b1 + 10*b4;
		contrast 'diff at 7 days sev, mod vs mild' b2 + 7*b5, b1 + 7*b4;

	*predict p out=pred;
	predict exp(b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev)/ (1+ exp(b0 + b1*mod+ b2*sev 
			+ b3*days + b4*days*mod + b5*days*sev + b6*day10 + b7*day10*mod + b8*day10*sev)) out=pred;
	run;
	proc sort data=pred; by days;
	run;
	data plot;
		set pred;
		ra=r/n;
	run;

	ods html style=journal;

	title1 h=2 "BPD status";
	proc sgpanel data=plot;
		panelby bpd1/rows=1 columns=3 novarname headerattrs=(size=16) spacing=4;
		format bpd1 b.;
		series y=ra x=days/group=study_id;
		pbspline y=pred x=days/nomarkers lineattrs=(thickness=4);
	rowaxis label="RA Staphylococcus" labelattrs=(size=16) valueattrs=(size=14);
	colaxis label= "Days" labelattrs=(size=16) valueattrs=(size=14) values=(0 to 28 by 7);
run;





data model2;
	set bpd.mergedwbeta2015_07_07;
		where taxa='Ureaplasma' and subj1 = 1;
	y=seq_count;
	days=days_from_birth;
	ltotal=log(total);
	n= total;
	r=seq_count;
	day10=max(0, days -10);
	sev=(BPD1=4);
	mod=(bpd1=3);
	lra=log(perc_total);
		steroids=(corticosteroids=1);
	keep study_id days y r n ltotal day10 shannonH_mean bpd1 sev mod lra daysmv steroids;
run;
proc sort data=model2;
	by study_id days;
run;

proc genmod data=model2;
	class study_id bpd1;
	model y= bpd1 bpd1*days bpd1*day10/ dist=nb offset=ltotal type3;
	repeated subject=study_id;
run;

proc mixed data=model2;
	class study_id bpd1;
	model lra=bpd1 bpd1*days/solution;
	random intercept/ subject=study_id ;
	repeated/ subject=study_id group=bpd1;
run;
	/*Beta-binomial - random intercept*/
	proc nlmixed data=model2 tech=quanew;*tech=newrap;
		parms b0=-1 b1=.6 b2=1 b3=0.2 b4=-0.15 b5=-.2 rho=.6 s2u=.2;* c0=1 c1=.8 c2=1.5;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + ui;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		*rho=exp(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random ui ~ NORMAL(0, s2u*s2u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 21days' b0 + 21*b3;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5;
		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5, b1 + 21*b4;

	run;
	/*join point model*/
		proc nlmixed data=model2 tech=quanew;*tech=newrap;
		parms b0=-1 b1=.9 b2=1.3 b3=0.2 b4=-0.15 b5=-.2 b6=-.3 b7=.2 b8=.2  s1u=.2 c0=.6 c1=0 c2=0;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev + u1;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		rho=(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random u1 ~ NORMAL(0, s1u*s1u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		contrast 'slope10 by group' b7, b8;

		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 10days' b0 + 10*b3;
		estimate 'mod 10 days' b0 + b1 + 10*b3 + 10*b4;
		estimate 'sev 10 days' b0 + b2 + 10*b3 + 10*b5;
		estimate 'mild 14days' b0 + 14*b3 +4*b6;
		estimate 'mod 14 days' b0 + b1 + 14*b3 + 14*b4 + 4*b6 + 4*b7;
		estimate 'sev 14days' b0 + b2 + 14*b3 + 14*b5 + 4*b6 + 4*b8;
		estimate 'mild 21days' b0 + 21*b3 +11*b6;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4 + 11*b6 + 11*b7;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5 + 11*b6 + 11*b8;

		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 10 days sev vs mild' b2 + 10*b5;
		estimate 'diff at 14 days sev vs mild' b2 + 14*b5 + 4*b8;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5 + 11*b8;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5 + 11*b8, b1 + 21*b4 + 11*b7;
		contrast 'diff at 14 days sev, mod vs mild' b2 + 14*b5 + 4*b8, b1 + 14*b4 + 4*b7;
		contrast 'diff at 10 days sev, mod vs mild' b2 + 10*b5, b1 + 10*b4;
		contrast 'diff at 7 days sev, mod vs mild' b2 + 7*b5, b1 + 7*b4;

	*predict p out=pred;
	predict exp(b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev)/ (1+ exp(b0 + b1*mod+ b2*sev 
			+ b3*days + b4*days*mod + b5*days*sev + b6*day10 + b7*day10*mod + b8*day10*sev)) out=pred2;
	run;
	/*sensitivity analysis with clinical covariates*/
	proc nlmixed data=model2 tech=quanew;*tech=newrap;
		parms b0=-1 b1=.9 b2=1.3 b3=0.2 b4=-0.15 b5=-.2 b6=-.3 b7=.2 b8=.2  s1u=.2 c0=.6 c1=0 c2=0 b9=.01 b10=0.1;
		lp =b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev + b9*daysmv + b10*steroids + u1;
		*s2u=exp(a0 + a1*mod +a2*sev);
	    p = exp( lp ) / (1+exp( lp ));
		rho=(c0 + c1*mod + c2*sev);
		M=(1-rho)/rho;
		loglike=(lgamma(n + 1)-lgamma(r + 1)-lgamma(n - r + 1))
      + lgamma( M ) - lgamma( M*p ) - lgamma( M*(1-p) )
      + lgamma( r + M*p ) + lgamma( n - r + M*(1-p) ) - lgamma( n + M );
  		model r ~ general (loglike);
		random u1 ~ NORMAL(0, s1u*s1u) subject=study_id;
		contrast 'int by group' b1, b2;
		contrast 'slope by group' b4, b5;
		contrast 'slope10 by group' b7, b8;
		estimate 'mild 7days' b0 + 7*b3;
		estimate 'mod 7 days' b0 + b1 + 7*b3 + 7*b4;
		estimate 'sev 7days' b0 + b2 + 7*b3 + 7*b5;
		estimate 'mild 10 days' b0 + 10*b3;
		estimate 'mod 10 days' b0 + b1 + 10*b3 + 10*b4;
		estimate 'sev 10 days' b0 + b2 + 10*b3 + 10*b5;
		estimate 'mild 14days' b0 + 14*b3 +4*b6;
		estimate 'mod 14 days' b0 + b1 + 14*b3 + 14*b4 + 4*b6 + 4*b7;
		estimate 'sev 14days' b0 + b2 + 14*b3 + 14*b5 + 4*b6 + 4*b8;
		estimate 'mild 21days' b0 + 21*b3 +11*b6;
		estimate 'mod 21 days' b0 + b1 + 21*b3 + 21*b4 + 11*b6 + 11*b7;
		estimate 'sev 21days' b0 + b2 + 21*b3 + 21*b5 + 11*b6 + 11*b8;

		estimate 'diff at 7 days sev vs mild' b2 + 7*b5;
		estimate 'diff at 10 days sev vs mild' b2 + 10*b5;

		estimate 'diff at 14 days sev vs mild' b2 + 14*b5 + 4*b8;
		estimate 'diff at 21 days sev vs mild' b2 + 21*b5 + 11*b8;
		contrast 'diff at 21 days sev, mod vs mild' b2 + 21*b5 + 11*b8, b1 + 21*b4 + 11*b7;
		contrast 'diff at 14 days sev, mod vs mild' b2 + 14*b5 + 4*b8, b1 + 14*b4 + 4*b7;
		contrast 'diff at 10 days sev, mod vs mild' b2 + 10*b5, b1 + 10*b4;
		contrast 'diff at 7 days sev, mod vs mild' b2 + 7*b5, b1 + 7*b4;

	*predict p out=pred;
	predict exp(b0 + b1*mod+ b2*sev + b3*days + b4*days*mod + b5*days*sev + 
			b6*day10 + b7*day10*mod + b8*day10*sev)/ (1+ exp(b0 + b1*mod+ b2*sev 
			+ b3*days + b4*days*mod + b5*days*sev + b6*day10 + b7*day10*mod + b8*day10*sev)) out=pred2;
	run;
	proc sort data=pred2; by days;
	run;
	data plot2;
		set pred2;
		ra=r/n;
	run;

		ods html style=journal;
		title1 h=2 "BPD status";
	proc sgpanel data=plot2;
		panelby bpd1/rows=1 columns=3 novarname headerattrs=(size=16) spacing=4 ;
		format bpd1 b.;
		series y=ra x=days/group=study_id;
		pbspline y=pred x=days/nomarkers lineattrs=(thickness=4);
	rowaxis label="RA Ureaplasma" labelattrs=(size=16) valueattrs=(size=14);
	colaxis label= "Days" labelattrs=(size=16) valueattrs=(size=14) values=(0 to 28 by 7);
run;





	ods html style=journal;

	title1 h=2 "BPD status";
	proc sgpanel data=plot;
		panelby bpd1/rows=1 columns=3 novarname headerattrs=(size=16) spacing=4;
		format bpd1 b.;
		series y=ra x=days/group=study_id;
		pbspline y=pred x=days/nomarkers lineattrs=(thickness=4);
	rowaxis label="RA Staphylococcus" labelattrs=(size=16) valueattrs=(size=14);
	colaxis label= "Days" labelattrs=(size=16) valueattrs=(size=14) values=(0 to 28 by 7);
run;
	proc sgpanel data=plot2;
		panelby bpd1/rows=1 columns=3 novarname headerattrs=(size=16) spacing=4 ;
		format bpd1 b.;
		series y=ra x=days/group=study_id;
		pbspline y=pred x=days/nomarkers lineattrs=(thickness=4);
	rowaxis label="RA Ureaplasma" labelattrs=(size=16) valueattrs=(size=14);
	colaxis label= "Days" labelattrs=(size=16) valueattrs=(size=14) values=(0 to 28 by 7);
run;
















/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/
/******************************/

/* time lag analysis*/
PROC IMPORT OUT= mh1  
            DATAFILE= "C:\Users\..\Data\Raw Data\BPD_Master_Morisita-Horn_Morisita-Horn_in_Workspace_Workspace_1.txt" 
            DBMS=TAB REPLACE;
RUN;
proc transpose data=mh1 out=plot1a;
proc transpose data=plot1a out=plot1b;
	id _NAME_;
	copy _NAME_;
run;
data plot_names;
	set plot1b;
	i=_N_;
	var2=_NAME_;
	keep var2 i;
run;
data plot1c;
	set plot1b;
	j=_n_;
	array x{*} _NUMERIC_;
	do i=1 to dim(x);
		mh=x{i};
		output;
    end;
	keep _NAME_ i j mh;
run;
proc sort data=plot1c;
	by i;
run;
proc sort data=plot_names;
	by i;
run;
data bpd_mhdata;
	merge plot_names plot1c;
	by i;
run;
proc print data=bpd_mhdata (obs=10);
run;
proc freq data=bpd_mhdata;
	table var2 _NAME_;
run;
data meta2a;
	set longitudinal;
	where seq < 2;
	var2=library; 
	keep library study_id days_from_birth bpd1 var2 ;
run;
proc sort data=bpd_mhdata;
	by var2;
run;
proc sort data=meta2a;
	by var2;
run;
data meta3;
	merge bpd_mhdata(where=(var2 ne ' ')) meta2a (rename = (study_id = id1 days_from_birth = days1));
	by var2;
run;
data meta2b;
	set longitudinal;
	where seq < 2;
	_NAME_=library; 
	keep library study_id days_from_birth bpd1 _NAME_ ;
run;
proc sort data=meta3;
	by _NAME_;
run;
proc sort data=meta2b;
	by _NAME_;
run;
data bpd.MH_timeLag;
	merge meta3 meta2b (rename = (study_id = id2 days_from_birth = days2));
	by _NAME_;
	if days1 =. then delete;
	if days2 = . then delete;
	if id1 ne id2 then delete;
	tlag=days2-days1;
	if tlag < 0 then delete;
	if var2=_NAME_ then delete;
	/*there is one repeat that needs to be deleted*/
	if i > j and tlag=0 then delete;
run;
proc sort data=bpd.MH_timeLag;
	by var2 days1 days2;
run;
proc sort data=bpd.MH_timeLag;
	by tlag;
run;
proc sgplot data=bpd.mh_timelag;
	vbox mh/ category=bpd1 group=id1;
run;
proc capability data=bpd.mh_timelag;
	var mh;
	comphistogram mh/class=bpd1;
	inset n median;
run;
/***********/
/* Figure 6*/
/***********/
ods html reset=all;
	ods pdf file='C:\Users\brandie\Documents\Childrens Hospital\BPD\Gerber microbiome\paper 1 fig3.pdf' dpi=500;

proc sgplot data=bpd.mh_timelag noautolegend;
	vbox mh/ category=bpd1 group=bpd1;
	format bpd1 b.;
	yaxis label="Morisita-Horn Beta Diversity" labelattrs=(size=16) valueattrs=(size=14);
	xaxis label= "BPD Status" labelattrs=(size=16) valueattrs=(size=14);
run;
ods pdf close;


ods html style=journal;
title1 'BPD status';
proc sgpanel data=bpd.mh_timelag noautolegend;
	panelby bpd1/ rows=3 columns=1 novarname;
	histogram mh;
	format bpd1 b.;
	colaxis label="Morisita-Horn Beta Diversity" labelattrs=(size=16) valueattrs=(size=14);
	rowaxis label= "Percent" labelattrs=(size=16) valueattrs=(size=14);
run;
ods graphics on;
proc genmod data=bpd.mh_timelag plots=(all);
	class id1 bpd1;
	model mh=bpd1/ noint link=log;
	repeated subject=id1;
	estimate 'mild' bpd1 1 0 0;
	estimate 'mod' bpd1 0 1 0;
	estimate 'sev' bpd1 0 0 1;
	estimate 'mod vs mild' bpd1 -1 1 0;
	estimate 'sev vs mild' bpd1 -1 0 1; 
run;
proc glimmix data=bpd.mh_timelag plots=(all);
	class id1 bpd1;
	model mh=bpd1/dist=beta noint;
	random intercept/ subject=id1;
	estimate 'mild' bpd1 1 0 0;
	estimate 'mod' bpd1 0 1 0;
	estimate 'sev' bpd1 0 0 1;
	estimate 'mod vs mild' bpd1 -1 1 0;
	estimate 'sev vs mild' bpd1 -1 0 1; 
run;
proc glimmix data=bpd.mh_timelag plots=(all);
	class id1 bpd1;
	model mh=bpd1/dist=normal link=log noint;
	random intercept/ subject=id1;
	estimate 'mild' bpd1 1 0 0;
	estimate 'mod' bpd1 0 1 0;
	estimate 'sev' bpd1 0 0 1;
	estimate 'mod vs mild' bpd1 -1 1 0;
	estimate 'sev vs mild' bpd1 -1 0 1; 
run;
proc sgpanel data=bpd.mh_timeLag;
	panelby id1;
	pbspline y=mh x=tlag;
run;
proc sgplot data=bpd.MH_timeLag;
	pbspline y=mh x=tlag/group=id1 smooth=9000000;
run;
proc sgplot data=bpd.MH_timeLag;
	pbspline y=mh x=tlag/group=bpd1 smooth=9000000;
run;
/***********/
/*figure E3*/
/***********/
proc sort data=bpd.MH_timeLag;
	by tlag;
run;
title 'BPD Status';
proc sgpanel data=bpd.MH_timeLag noautolegend;
	panelby bpd1 / columns=3 novarname;
	pbspline y=mh x=tlag/ smooth=9000000;
	series y=mh x=tlag/group=id1;
	format bpd1 b.;
	colaxis label='Days Between Specimen Collection';
	rowaxis label='Morisita-Horn Beta Diversity';
run;

ods graphics on;
proc genmod data=bpd.MH_timeLag ;
	class id1 bpd1 ;
	model mh=bpd1 tlag*bpd1/link=log noint type3;
	repeated subject=id1;
	contrast "across all 3" bpd1 -1 0 1 tlag*bpd1 -1 0 1, bpd1 -1 1 0 tlag*bpd1 -1 1 0;
	contrast "mild vs sev" bpd1 -1 0 1 tlag*bpd1 -1 0 1;
	contrast "mild vs mod" bpd1 -1 1 0 tlag*bpd1 -1 1 0;
run;


proc sgpanel data=bpd.MH_timeLag;
	panelby bpd1/ columns=3 rows=1;
	pbspline y=mh x=days2/smooth=9000000;
	*rowaxis type=log logstyle=logexpand logbase=2;
run;

