/*******************************************************/
/* Macro: InfValue
/*******************************************************/
%macro InfValue(DSin, XVar, YVarBin, M_IV);

/* Extract the frequency table using proc freq, 
   and the categories of the X variable */

proc freq data=&DSin noprint;
 table &XVar*&YvarBin /out=Temp_freqs;
 table &XVar /out=Temp_Xcats;
 run;

proc sql noprint;
  /* Count the number of obs and categories of X */
   %local R C; /* rows and columns of freq table */
   select count(*) into : R from temp_Xcats;
   select count(*) into : N from &DSin; 
quit;

  /* extract the categories of X into CatX_i */
data _Null_;
  set temp_XCats;
   call symput("CatX_"||compress(_N_), &Xvar);
run;

proc sql noprint; 
	/* extract n_i_j*/
%local i j;
   %do i=1 %to &R; 
    %do j=1 %to 2;/* we know that YVar is 1/0 - numeric */
      %local N_&i._&j;
   Select Count into :N_&i._&j from temp_freqs where &Xvar ="&&CatX_&i" and &YVarBin = %eval(&j-1);
    %end;
   %end;
quit;
  
  /* calculate N*1,N*2 */
     %local N_1s N_2s;
      %let N_1s=0;
	  %let N_2s=0;
  %do i=1 %to &r; 
	  %let N_1s=%sysevalf(&N_1s + &&N_&i._1);
	  %let N_2s=%sysevalf(&N_2s + &&N_&i._2);
   %end;

/* substitute in the equation for IV */
     %local IV;
     %let IV=0;
       %do i=1 %to &r;
          %let IV = %sysevalf(&IV + (&&N_&i._1/&N_1s - &&N_&i._2/&N_2s)*%sysfunc(log(%sysevalf(&&N_&i._1*&N_2s/(&&N_&i._2*&N_1s)))) );
       %end;

%let &M_IV=&IV; 

/* clean the workspace */
proc datasets library=work;
delete temp_freqs temp_Xcats;
quit;
%mend;




/*******************************************************/
/* Macro : PowerIV */
/*******************************************************/
%macro PowerIV(DSin, DV, IVList, DSout);
/* Decompose the input IVList into tokens and store variable
   names into macro variables */

%local i N condition VarX; 
%let i=1;
%let N=0;
%let condition = 0; 
%do %until (&condition =1);
   %let VarX=%scan(&IVList,&i);
   %if "&VarX" =""  %then %let condition =1;
  	        %else %do;
				%local Var&i;
                %let Var&i =&VarX; 
                %let N=&i;
                %let i=%eval(&i+1); 
                  %end;  
%end;

/* now we have a total of N variables
   Loop on their  names and calculate the Information value
   between the DV and each of the variables */

proc sql noprint;
 create table &DSout (VariableName char(200), 
                      InformationValue  num);
quit;

%do i=1 %to &N;
   %local IV&i;
   %let IV&i=;
	%InfValue(&DSin, &&Var&i, &DV, IV&i);
	proc sql noprint; 
     insert into &DSout  values("&&Var&i",&&IV&i);
    quit; 	 
%end;


proc sort data=&dsout;
 by descending InformationValue; 
 run;

%mend; 



/*******************************************************/
/* Macro : GrFBinDV */
/*******************************************************/


%macro  GrFBinDV(DSin, Xvar, BinDV, M_Gr, M_Fstar, M_Pvalue);
/* Calculation of the Gr and F* values for a continuous 
   variable Xvar and a binary DV
	DSin = input dataset
	BinDV = binary dependent variable (1/0 only)
	XVar = X Variable - continuous
	M_GR = returned Gr Gini ratio
	M_Fstar= returned F* 
	M_Pvalue= returned p-value of F*
*/


/* Get the categories of the BinDV*/
	proc freq data=&DSin noprint ;
	tables &BinDV /missing out=Temp_Cats;
	run;

	/* Convert the categories (Y_i) and their frequencies 
	  n_i to macro variables */

	%local N_1 N_2 N;

	Data _null_;
	retain N 0;
	  set Temp_Cats;
	   N=N+count;
       	       call symput ("n_" || left(_N_), left(count));
		   call symput ("N", left(N));
        Run;

	/* Calculate the quantities needed to substitute in 
	   SSTO, SSR, SSE, MSR, MSE, F*, Gr */

   proc sql noprint;
   /* xbar */
   %local xbar i; 
    select avg(&xVar) into :xbar from &DSin;
 
   %do i=1 %to 2;
	/* Ybar_i */
	select avg(&XVar) into :Xbar_&i 
	          from &DSin where &BinDV = %eval(&i-1);  

   %end;	

	/* SSTO, SSR, SSE */
   %local SSTO SSR SSE;
    select var(&XVar) into: SSTO from &DSin;
	%let SSTO=%sysevalf(&SSTO *(&N-1));
	%let SSR=0;
	%let SSE=0;

	%local SSEi;
    %do i=1 %to 2;
	  select var(&XVar) into: SSEi 
	            from &DSin where &BinDV=%eval(&i-1);
      %let SSE=%sysevalf(&SSE + &ssei * (&&n_&i - 1)) ; 

	  %let SSR=%sysevalf(&SSR+ &&n_&i * (&&Xbar_&i - &Xbar)*
	                                    (&&Xbar_&i - &Xbar));
    %end;

  quit; /* end of Proc SQL */

	/* MSR, MSE , F*, Gr, Pvalue */
    %local MSR MSE Fstar;
	%let MSR=&SSR;
	%let MSE=%sysevalf(&SSE/(&N-2));
	%let &M_Gr=%Sysevalf(1-(&SSE/&SSTO));
	%let Fstar=%sysevalf(&MSR/&MSE);
	%let &M_Fstar=&Fstar;
	%let &M_PValue=%sysevalf(%sysfunc(probf(&Fstar,1,&N-2)));

/* clean workspace */
	proc datasets library=work nolist;
	 delete temp_cats;
	run; quit;
%mend;



/*******************************************************/
/* Macro : PowerFG */
/*******************************************************/


%macro PowerFG(DSin, DV, IVList, DSout);

/* Decompose the input IVList into tokens and store variable
   names into macro variables */

%local i N condition VarX; 
%let i=1;
%let N=0;
%let condition = 0; 
%do %until (&condition =1);
   %let VarX=%scan(&IVList,&i);
   %if "&VarX" =""  %then %let condition =1;
  	        %else %do;
				%local Var&i;
                %let Var&i =&VarX; 
                %let N=&i;
                %let i=%eval(&i+1); 
                  %end;  
%end;


proc sql noprint;
 create table &DSout (VariableName char(200),
					  GiniRatio         num ,
                      FStar             num ,
                      Fstar_Pvalue      num);
quit;


%do i=1 %to &N;
   %local Gr&i Fstar&i pvalue&i;
   %let Gr&i=;
   %let Fstar&i=;
   %let Pvalue&i=;
    %GrFBinDV(&DSin, &&Var&i, &DV, Gr&i, Fstar&i, Pvalue&i);

proc sql noprint; 
insert into &DSout  values("&&Var&i", &&Gr&i, &&Fstar&i, &&Pvalue&i);
    quit; 	 
%end;


%mend; 



/*******************************************************/
/* Macro : ExtractTop */
/*******************************************************/

%macro ExtrctTop(DSin, VarCol, SelVar, Method, NTop, CutOff, M_VarList);

/* sort the dataset using the selection criteria in descending order */
proc sort data=&Dsin;
by descending &SelVar;
run;


%local Nact; 

 data _null_;
  set &DSin;
   by descending &SelVar;
%if (&Method=1) %then %do;
        if (_N_ le &NTop) then do;
		call symput("Var"||compress(_N_), &VarCol); 
		call symput("Nact", compress(_N_));
							  end;
%end;

%else %do;
	if (&SelVar ge &CutOff) then  do;
      call symput("Var"||compress(_N_), &VarCol); ; 
	  call symput("Nact", compress(_N_));
	                              end;
%end;

run;

 /* initialize the list and compose it using the extracted names*/
%local List;
%let List=;
 %local i;
  %do i=1 %to &Nact;
    %let List=&List &&Var&i;
  %end;

  %let &M_VarList=&List; 


%mend;


/*******************************************************/
/* Load the dataset */
/*******************************************************/
%let dir=E:\SAS\;
data credit;
set "&dir credit_mini.sas7bdat";
run;
proc contents data=credit;
run;
proc sql noprint;
select name into :varname_char separated by ' ' 
                from dictionary.columns
                where libname=upcase("work") and memname=upcase("credit") and type="char";
        %put &varname_char;
quit;
proc sql noprint;
select name into :varname_num separated by ' ' 
                from dictionary.columns
                where libname=upcase("work") and memname=upcase("credit") and type="num";
        %put &varname_num;
quit;


/*******************************************************/
/* Call the macro */
/*******************************************************/
%let DSin=Credit;
%let DV=bad_good;
%let IVList=&varname_char;
%let DSOut=Credit_IVs;

%PowerIV(&DSin, &DV, &IVList, &DSout);

/*******************************************************/
/* Extract the top variables with IV>=0.015*/
/*******************************************************/

%let DSin=Credit_IVs;
%let VarCol=VariableName;
%let SelVar=InformationValue; 
%let Method=2;
%let Ntop=0;
%let Cutoff=0.015;
%let VarList=;
%ExtrctTop(&DSin, &VarCol, &SelVar, &Method, &NTop, &CutOff, VarList);
%put Selected Variables: &VarList;

/*******************************************************/
/* Call the macro */
/*******************************************************/
%let DSin=Credit;
%let DV=bad_good;
%let IVList=&varname_num;
%let DSOut=Credit_FGs;

%PowerFG(&DSin, &DV, &IVList, &DSout);

/*******************************************************/
/* Extract the top variables fist 10*/
/*******************************************************/

%let DSin=Credit_FGs;
%let VarCol=VariableName;
%let SelVar=FStar; 
%let Method=1;
%let Ntop=11;
%let Cutoff=0;
%let VarList=;
%ExtrctTop(&DSin, &VarCol, &SelVar, &Method, &NTop, &CutOff, VarList);
%put Selected Variables: &VarList;






