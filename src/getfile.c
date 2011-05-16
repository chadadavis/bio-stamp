/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stamp.h"

#define MAXLINE 200 

/* getfile: 
 *  opens a file containing a list of appropriate directories, suffixes and
 *  prefixes.  It then makes combinates of these with the supplied code
 *  to see if a file exists 
 *
 * Allocates memory for a string corresponding to the filename read in,
 *  then returns this string or '\0' if none is found 
 *
 * Modification March 1999 - now tries to expand into a "distr" style hierarchy */

char *getfile(char *code_safe, char *dirfile, int code_length, FILE *OUTPUT) {

	char *temp;
	char *dir;
	char *buff;
	char prefix[10],suffix[10];
	char *code;
	char *CODE;
	char *two_letter; 
	char a;

	int i,j;
	int found;

	FILE *DIR;

	dir=(char*)malloc(1000*sizeof(char)); 
	buff=(char*)malloc(1000*sizeof(char)); 
	temp=(char*)malloc(1000*sizeof(char)); 
	CODE=(char*)malloc((strlen(code_safe)+1)*sizeof(char));
	code=(char*)malloc((strlen(code_safe)+1)*sizeof(char));
	two_letter = (char*)malloc(3*sizeof(char));

	/* copy code_safe into code to prevent the changing case from
	 *  carrying on through the program */

	a = code[code_length];
	code[code_length]='\0'; /* truncates the code if necessary */
	for(i=0; i<strlen(code_safe); ++i) code[i]=utol(code_safe[i]);
	for(i=0; i<strlen(code_safe); ++i) CODE[i]=ltou(code_safe[i]); /* upper case code */
	code[code_length]='\0';
	CODE[code_length]='\0';
	two_letter[0] = code[1]; two_letter[1] = code[2]; two_letter[2] = '\0'; /* For PDB distr style */

	if((DIR=fopen(dirfile,"r"))==NULL) {
	   fprintf(stderr,"error: no directories file found\n");
	   temp[0]='\0';
	   return temp;
	}
	found=0;
	while((fgets(buff,MAXLINE,DIR))!=NULL) {

	  sscanf(buff,"%s %s %s",dir,&prefix[0],&suffix[0]);
	  if(prefix[0] == '_') { prefix[0] = '\0'; }
	  if(suffix[0] == '_') { suffix[0] = '\0'; }

	  /* Lower case, all files one directory */
	  sprintf(temp,"%s/%s%s%s",dir,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Lower case, files in 'distr' type directories */
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,code,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, all files one directory */
	  sprintf(temp,"%s/%s%s%s",dir,prefix,CODE,suffix);
	  if(testfile(temp)) { found=1; break; } 

	  /* Upper case, files in 'distr' type directories */
	  sprintf(temp,"%s/%s/%s%s%s",dir,two_letter,prefix,CODE,suffix);
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
	free(two_letter);

	return temp;
}
