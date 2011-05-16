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
#include "stamp.h"

/* Get secondary structures from a file */

int getsec(struct domain_loc *domain, int ndomain, struct parameters *parms) {

	char c;
	char id[100];
	char *sec;
	char *empty;

	int i,j,k;
	int which;

	int *found;

	FILE *IN;

	if((IN=fopen(parms[0].secfile,"r"))==NULL) {
	   fprintf(stderr,"error: file %s does not exist\n",parms[0].secfile);
	   for(i=0; i<ndomain; ++i) {
	     for(j=0; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
	     domain[i].sec[j]='\0';
	   }
	   return -1;
	}
	fprintf(parms[0].LOG,"Searching for secondary structure assignments in the file %s\n",
		parms[0].secfile);
	sec=(char*)malloc(parms[0].MAX_SEQ_LEN*sizeof(char));
	found=(int*)malloc(ndomain*sizeof(int));

	for(i=0; i<ndomain; ++i) found[i]=0;
	
	while((c=getc(IN))!=(char)EOF) {
	   if(c=='>') {  
	      if(fscanf(IN,"%s",&id[0])==(int)EOF) return -1; 
	      which=-1;
	      for(i=0; i<ndomain; ++i) 
		 if(strcmp(id,domain[i].id)==0 && !found[i]) { which=i; found[i]=1; break; }
	      while((c=getc(IN))!=(char)EOF && c!='\n') if(c==(char)EOF) break;
	      while((c=getc(IN))!=(char)EOF && c!='\n') if(c==(char)EOF) break;
	      i=0;
	      while((c=getc(IN))!=(char)EOF && c!='*') { 
		if(c==(char)EOF) break;
		if(c!='\n') sec[i++]=c;
		if(i>parms[0].MAX_SEQ_LEN) break;
	      }
	      sec[i]='\0';
	      if(which!=-1) strcpy(domain[which].sec,sec);
           }
	}
	/* now see how many we missed */
	for(i=0; i<ndomain; ++i) {
	   if(found[i]) {
	      fprintf(parms[0].LOG,"\nFound secondary structure for %s\n",domain[i].id);
	      if(strlen(domain[i].sec)!=strlen(domain[i].aa)) {
		fprintf(parms[0].LOG," though the lengths differ, ");
		if(strlen(domain[i].sec)>strlen(domain[i].aa)) {
		   domain[i].sec[strlen(domain[i].aa)]='\0';
		   fprintf(parms[0].LOG,"truncating\n");
		} else {
		  for(j=strlen(domain[i].sec); j<strlen(domain[i].aa); ++j) domain[i].sec[j]='?';
		  domain[i].sec[strlen(domain[i].aa)]='\0';
		  fprintf(parms[0].LOG,"padding with ?s\n");
		}
		fprintf(parms[0].LOG," warning: assignment of secondary structure to sequence may be incorrect\n");
	     }
	   } else {
	      fprintf(parms[0].LOG,"no secondary structure found for %s\n",domain[i].id);
	      for(j=0; j<strlen(domain[i].aa); ++j) domain[i].sec[j]='?';
	   }
	   display_align(&domain[i].aa,1,&domain[i].sec,1,&domain[i].aa,&domain[i].sec,
	       empty,empty,parms[0].COLUMNS,0,0,parms[0].LOG);

	}
	return 0;
}
