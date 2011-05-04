/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char **RBR_c_split(char *str, int *n,  char delimiter) {
	/* breaks a string up into substrings using a delimiter */

	int i,j,k,len;
	char **sstr;

	sstr=(char**)malloc(sizeof(char*));
	
	(*n)=0;
	i=0;
	len=strlen(str);

 	while(i<len && str[i]!='\0' && str[i]!='\n') {
	   while(i<len && str[i]==delimiter) i++;
	   if(i>=len || str[i]=='\0' || str[i]=='\n') break;
	   /* we are at a new sub-string */
	   sstr[(*n)]=(char*)malloc(sizeof(char));
	   j=0;
	   while(i<len && str[i]!='\0' && str[i]!='\n' && str[i]!=delimiter) {
		sstr[(*n)][j]=str[i];
		j++; i++;
	        sstr[(*n)]=(char*)realloc(sstr[(*n)],(j+1)*sizeof(char));
	   }
	   sstr[(*n)][j]='\0';
	   (*n)++;
	   sstr=(char**)realloc(sstr,((*n)+1)*sizeof(char*));
	   i++;
	}
	return sstr;
}
