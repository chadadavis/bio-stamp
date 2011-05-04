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
