/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include "dstamp.h"

int readcom(n,args,parms)
int n;
char **args;
struct parameters *parms;
{
	int i,j;
	int ON;
	
	FILE *IN;
	if(n<2) exit_error();
	for(i=1; i<n; ++i) {
           if(args[i][0]!='-') exit_error();
	   if((i+1)<n && (args[i+1][0]=='n' || args[i+1][0]=='N' || args[i+1][0]=='T' || args[i+1][0]=='t' || args[i+1][0]=='1' || args[i+1][0]=='O' || args[i+1][0]=='o')) ON=1;
	   else ON=0;
	   for(j=1; j<strlen(args[i]); ++j) args[i][j]=utol(args[i][j]);
           if(strcmp(&args[i][1],"f")==0) {
	     if((i+1)>=n) exit_error();
	     strcpy(&parms[0].filename[0],args[i+1]);
             i++;
           } else if(strcmp(&args[i][1],"c")==0) {
             if((i+1)>=n) exit_error();
             parms[0].TYPE=args[i+1][0];
             i++;
           } else if(strcmp(&args[i][1],"t")==0) {
             if((i+1)>=n) exit_error();
             sscanf(args[i+1],"%f",&parms[0].CUTOFF);
             i++;
           } else if(strcmp(&args[i][1],"w")==0) {
             if((i+1)>=n) exit_error();
             sscanf(args[i+1],"%d",&parms[0].WINDOW);
             i++;
	   } else if(strcmp(&args[i][1],"prefix")==0) {
             if((i+1)>=n) exit_error();
             strcpy(parms[0].prefix,args[i+1]);
	     i++;
	   } else if(strcmp(&args[i][1],"P")==0) {
	     if((IN=fopen(args[i+1],"r"))==NULL ) {
	       fprintf(stderr,"error: parameter file %s does not exist\n",args[3]);
	       exit(-1);
	     }
	     if(getpars(IN,parms)==-1) exit(-1);
	     fclose(IN);
	     i++;
	   } else if(strcmp(&args[i][1],"cutoff")==0) {
		if((i+1)>=n) exit_error();
		sscanf(args[i+1],"%f",&parms[0].CUTOFF);
		i++;
	   } else if(strcmp(&args[i][1],"window")==0) {
		if((i+1)>=n) exit_error();
		sscanf(args[i+1],"%d",&parms[0].WINDOW);
		i++;
	   } else if(strcmp(&args[i][1],"type")==0) {
		if((i+1)>=n) exit_error();
		sscanf(args[i+1],"%s",&parms[0].TYPE);
		i++;
	   } else if(strcmp(&args[i][1],"pointsize")==0) {
		if((i+1)>=n) exit_error();
		sscanf(args[i+1],"%f",&parms[0].POINTSIZE);
		i++;
	   } else if(strcmp(&args[i][1],"boxseq")==0) {
		if((i+1)>=n) exit_error();
		parms[0].BOXSEQ=ON;
		i++;
	   } else if(strcmp(&args[i][1],"boxsec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].BOXSEC=ON;
		i++;
	   } else if(strcmp(&args[i][1],"smallseq")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SMALLSEQ=ON;
		i++;
	   } else if(strcmp(&args[i][1],"smallsec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SMALLSEC=ON;
		i++;
	   } else if(strcmp(&args[i][1],"boldseq")==0) {
		if((i+1)>=n) exit_error();
		parms[0].BOLDSEQ=ON;
		i++;
	   } else if(strcmp(&args[i][1],"boldsec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].BOLDSEC=ON;
		i++;
	   } else if(strcmp(&args[i][1],"sec")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SEC=ON;
		i++;
	   } else if(strcmp(&args[i][1],"secsum")==0) {
		if((i+1)>=n) exit_error();
		parms[0].SECSUM=ON;
		i++;
	   } else if(strcmp(&args[i][1],"verbose")==0) {
		if((i+1)>=n) exit_error();
		parms[0].VERBOSE=ON;
		i++;
	   } else if(strcmp(&args[i][1],"caseseq")==0) {
		if((i+1)>=n) exit_error();
 		parms[0].CASESEQ=ON;
		i++;
	   } else if(strcmp(&args[i][1],"casesec")==0) { 
		if((i+1)>=n) exit_error();
		parms[0].CASESEC=ON;
		i++;
           } else {
	     printf("unrecognised command: %s \n",&args[i][1]);
             exit_error();
           }
        }
	
	return 0;
}
int exit_error()
{
	  fprintf(stderr,"format: dstamp -f <STAMP file> -prefix <output prefix> \n");
 	  fprintf(stderr,"       (-P (<optional paramter file> -<parameter> <value>)\n");
	  exit(-1);
}
