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
#include "dstamp.h"

/* This routine works in conjunction with Agetbloc.
 *
 *  Given a blockfile containing a list of '>' which specifiy 
 *   character sequences, and a list of "#" which specify
 *   numeric data, this reads them in and returns them as
 *   a string of floats called 'stamp' */

int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos) {
	
	int i,j,k,n;
	float value;

	char c;
	char *tmp;

	tmp=(char*)malloc(1000*sizeof(char));

	(*nstamp)=0;
	(*nseq)=0;

	/* read in the descriptors */
	while((c=getc(IN))!=(char)EOF) {
	   (*nseq)+=(c=='>');
	   if(c=='#') {
	      fgets(tmp,100,IN);
	      stamp[(*nstamp)].what=tmp[0];
	      stamp[(*nstamp)].title=(char*)malloc(100*sizeof(char));
	      stamp[(*nstamp)].n=(float*)malloc(maxpos*sizeof(float));
	      strncpy(stamp[(*nstamp)].title,&tmp[1],99);
	      stamp[(*nstamp)].title[99]='\0';
	      (*nstamp)++;
	   }
	}
	rewind(IN);

	/* now find the first '*' and which column it is in */
	n=0;
	while((c=getc(IN))!=(char)EOF && c!='*') n=(n+(c!='\n'))*(c!='\n');
	if(c==(char)EOF) return -1;
	while((c=getc(IN))!=(char)EOF && c!='\n'); /* read to the end of the line */
	if(c==(char)EOF) return -1;

	/* now read in the file line by line */
	(*npos)=0;
	while(tmp[n]!='*') {
	  fgets(tmp,900,IN);
	  for(i=0; i<strlen(tmp); ++i) if(tmp[i]=='\n') tmp[i]='\0';
	  if(tmp[n]=='*') break;
	  /* Lets allow for missing values, and set them to zero.
	   * In other words, when the length of the string read in is
	   * equal to the number of '>' characters + n, ignore the line */
	  if(strlen(tmp)<=((*nseq)+n+1)) {
	     for(i=0; i<(*nstamp); ++i) 
		stamp[i].n[(*npos)]=-1.0;
	  } else { 
	     j=((*nseq)+n);
	     for(i=0; i<(*nstamp); ++i) {
		while(tmp[j]==' ') ++j;  /* move to the next space in the string */
		sscanf(&tmp[j],"%f",&value);
		stamp[i].n[(*npos)]=value; /* read in the next float */
		while(tmp[j]!=' ' && tmp[j]!='\0') ++j; /* move to the next space */
	     }
	  }
	  (*npos)++;
	}
	free(tmp);
	return 0;
}
