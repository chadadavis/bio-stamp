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
#include <malloc.h>
#include <math.h>

#include "include.h"

#define MAXLINE 100 

/* getfile: 
 *  opens a file containing a list of appropriate directories, suffixes and
 *  prefixes.  It then makes combinates of these with the supplied code
 *  to see if a file exists 
 *
 * Allocates memory for a string corresponding to the filename read in,
 *  then returns this string or '\0' if none is found */
char *getfile(code_safe,dirfile,code_length,OUTPUT)
char *code_safe;
char *dirfile;
int code_length;
FILE *OUTPUT;
{
	char *temp;
	char *dir;
	char *buff;
	char prefix[10],suffix[10];
	char *code;
	char *CODE;
	char a;

	int i,j;
	int found;

	FILE *DIR;

	dir=(char*)malloc(100*sizeof(char)); 
	buff=(char*)malloc(100*sizeof(char)); 
	temp=(char*)malloc(100*sizeof(char)); 
	CODE=(char*)malloc((strlen(code_safe)+1)*sizeof(char));
	code=(char*)malloc((strlen(code_safe)+1)*sizeof(char));

	/* copy code_safe into code to prevent the changing case from
	 *  carrying on through the program */

	a = code[code_length];
	code[code_length]='\0'; /* truncates the code if necessary */
	for(i=0; i<strlen(code_safe); ++i) code[i]=utol(code_safe[i]);
	for(i=0; i<strlen(code_safe); ++i) CODE[i]=ltou(code_safe[i]); /* upper case code */
	code[code_length]='\0';
	CODE[code_length]='\0';

	if((DIR=fopen(dirfile,"r"))==NULL) {
	   fprintf(stderr,"error: no directories file found\n");
	   temp[0]='\0';
	   return temp;
	}
	found=0;
	while((fgets(buff,MAXLINE,DIR))!=NULL) {
	  sscanf(buff,"%s %s %s",dir,&prefix[0],&suffix[0]);
	  sprintf(temp,"%s/",dir);
	  if(prefix[0]!='_') sprintf(&temp[strlen(temp)],"%s",prefix);
	  sprintf(&temp[strlen(temp)],"%s",code);
	  if(suffix[0]!='_') sprintf(&temp[strlen(temp)],"%s",suffix);
	  if(testfile(temp)) { found=1; break; } 
	  sprintf(temp,"%s/",dir);
          if(prefix[0]!='_') sprintf(&temp[strlen(temp)],"%s",prefix);
	  sprintf(&temp[strlen(temp)],"%s",CODE);
	  if(suffix[0]!='_') sprintf(&temp[strlen(temp)],"%s",suffix);
	  if(testfile(temp)) { found=1; break; }
	}

	if(found==0) {
	   temp[0]='\0';
	}
	fclose(DIR);
	free(dir); 
	free(buff);
	free(code);
	free(CODE);

	return temp;
}
