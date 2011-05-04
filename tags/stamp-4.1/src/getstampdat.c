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

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <dstamp.h>

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
