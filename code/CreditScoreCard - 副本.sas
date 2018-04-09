/*********************************************/
/*********************************************/
/***** Automatically Generated Scorecard *****/
/*********************************************/
/************    SAS CODE             ********/

/* Scorecard Scale : */
/*  Odds of [ 1 : 60 ] at  [ 600 ] Points 
     with PDO of [ 20 ] */

/*********************************************/
/*********************************************/

/********** START OF SCORING DATA STEP *******/
/*********************************************/
/*********************************************/

DATA SCORING;        /********** Modify ************/
 SET ScoringDataset; /********** Modify ************/

/*********************************************/
/*********************************************/
/*********************************************/
/* Base Points   */
/*********************************************/
Points=388 ;
/*********************************************/
/* Variable : CCCC     *****/
/*********************************************/
      IF CCCC LE (0.84) THEN  Points=Points +(-3);
      IF CCCC GT (0.84) AND CCCC LE (2.52) THEN  Points=Points +(27);
      IF CCCC GT (2.52) AND CCCC LE (7.56) THEN  Points=Points +(16);
      IF CCCC GT (7.56) AND CCCC LE (9.24) THEN  Points=Points +(27);
      IF CCCC GT (9.24) THEN  Points=Points +(34);
/*********************************************/
/* Variable : DSATD     *****/
/*********************************************/
      IF DSATD LE (2145.74) THEN  Points=Points +(1);
      IF DSATD GT (2145.74) AND DSATD LE (4354.59) THEN  Points=Points +(-2);
      IF DSATD GT (4354.59) AND DSATD LE (4480.81) THEN  Points=Points +(-1);
      IF DSATD GT (4480.81) AND DSATD LE (4985.69) THEN  Points=Points +(-6);
      IF DSATD GT (4985.69) THEN  Points=Points +(5);
/*********************************************/
/* Variable : DSDMOA     *****/
/*********************************************/
      IF DSDMOA LE (909862) THEN  Points=Points +(0);
      IF DSDMOA GT (909862) AND DSDMOA LE (8188758) THEN  Points=Points +(-7);
      IF DSDMOA GT (8188758) AND DSDMOA LE (9098620) THEN  Points=Points +(-31);
      IF DSDMOA GT (9098620) AND DSDMOA LE (13647930) THEN  Points=Points +(-3);
      IF DSDMOA GT (13647930) THEN  Points=Points +(-18);
/*********************************************/
/* Variable : DSFAC     *****/
/*********************************************/
      IF DSFAC LE (4.2) THEN  Points=Points +(0);
      IF DSFAC GT (4.9) AND DSFAC LE (11.2) THEN  Points=Points +(-8);
      IF DSFAC GT (11.9) AND DSFAC LE (15.05) THEN  Points=Points +(-1);
      IF DSFAC GT (15.75) AND DSFAC LE (18.2) THEN  Points=Points +(-20);
      IF DSFAC GT (18.9) THEN  Points=Points +(-22);
/*********************************************/
/* Variable : DSLTD     *****/
/*********************************************/
      IF DSLTD LE (3241.68) THEN  Points=Points +(0);
      IF DSLTD GT (3241.68) AND DSLTD LE (3927.42) THEN  Points=Points +(-6);
      IF DSLTD GT (3927.42) AND DSLTD LE (4052.1) THEN  Points=Points +(13);
      IF DSLTD GT (4052.1) AND DSLTD LE (4737.84) THEN  Points=Points +(-4);
      IF DSLTD GT (4737.84) THEN  Points=Points +(2);
/*********************************************/
/* Variable : DSOTD     *****/
/*********************************************/
      IF DSOTD LE (2171.58) THEN  Points=Points +(15);
      IF DSOTD GT (2171.58) AND DSOTD LE (3385.11) THEN  Points=Points +(-7);
      IF DSOTD GT (3385.11) AND DSOTD LE (3832.2) THEN  Points=Points +(-26);
      IF DSOTD GT (3832.2) AND DSOTD LE (5237.34) THEN  Points=Points +(-13);
      IF DSOTD GT (5237.34) THEN  Points=Points +(10);
/*********************************************/
/* Variable : L6SMMOA     *****/
/*********************************************/
      IF L6SMMOA LE (900012) THEN  Points=Points +(0);
      IF L6SMMOA GT (900012) AND L6SMMOA LE (2700036) THEN  Points=Points +(-6);
      IF L6SMMOA GT (2700036) AND L6SMMOA LE (12600168) THEN  Points=Points +(-15);
      IF L6SMMOA GT (12600168) AND L6SMMOA LE (15300204) THEN  Points=Points +(-4);
      IF L6SMMOA GT (15300204) THEN  Points=Points +(-21);
/*********************************************/
/* Variable :       *****/
/*********************************************/
      /*IF  = "N" THEN  Points=Points +(.);
      IF  = "Y" THEN  Points=Points +(.);*/
/*********************************************/
/* Variable : CIDF     *****/
/*********************************************/
      IF CIDF = "Y" THEN  Points=Points +(0);
      IF CIDF = "N" THEN  Points=Points +(-8);
/*********************************************/
/* Variable : CSFF     *****/
/*********************************************/
      IF CSFF = "Y" THEN  Points=Points +(21);
      IF CSFF = "N" THEN  Points=Points +(-1);
/*********************************************/
/* Variable : CSPF     *****/
/*********************************************/
      IF CSPF = "Y" THEN  Points=Points +(2);
      IF CSPF = "N" THEN  Points=Points +(-29);
/*********************************************/
/* Variable : DTF     *****/
/*********************************************/
      IF DTF = "Y" THEN  Points=Points +(23);
      IF DTF = "N" THEN  Points=Points +(-3);
/*********************************************/
/* Variable : GENDER     *****/
/*********************************************/
      IF GENDER = "1" THEN  Points=Points +(-1);
      IF GENDER = "2" THEN  Points=Points +(2);
      IF GENDER = "X" THEN  Points=Points +(-4);
RUN;

/**************按照分值大小降序排列*************/
proc sort data=scoring;
by descending points;
run;

/*************END OF SCORING DATA STEP *******/
/*********************************************/
