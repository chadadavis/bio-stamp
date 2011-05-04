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
#include <stamp.h>
#define MAXL 200    /* max length for input line */

/* get_dssp_sum:  return the secondary structure summary
	 from a dssp file */

int get_dssp_sum(FILE *DSSP, struct brookn begin, struct brookn end, int type,
	char *sec, int REVERSE,	int maxlen, int *count, FILE *OUT) {

    char start,datstart;
    char cid,in;
    char css;

    char number[6];
    char *line;

    int n,i;

    datstart=start=0;
    *count=0;

    line = (char *) malloc(sizeof(char)*MAXL);

    while(fgets(line,MAXL,DSSP) != NULL) {
	strncpy(&number[0],&line[5],5); number[5]='\0';
	sscanf(&number[0],"%d",&n);
	cid=line[11];
	in=line[10];
	if(datstart && (line[13]!='!') &&
	   ( (datstart && type==1) ||
	     (cid==begin.cid && type==2) ||
	     (cid==begin.cid && in==begin.in && n==begin.n && type==3) 
	    )) {
	      start=1;
	}
	if(datstart && start && type==2 && (line[13]!='!') && begin.cid!=cid) break;
	if(start && (line[13]!='!')) { /* skips ill placed chain breaks */
	    if((*count)==maxlen-1) {
	       fprintf(stderr,"error: maximum sequence length surpassed when attempting to get DSSP summary\n");
	       return -1;
	    }
	    sec[(*count)] = line[16];
	    if(sec[(*count)]==' ') sec[(*count)]='-';
	    (*count)++;
	}
	if(start && line[13]!='!' && cid==end.cid && in==end.in && n==end.n && type==3) break;
	if(line[2] == '#') datstart=1;
    }
    sec[(*count)]='\0';
    if(REVERSE) {
       for(i=0; i<(int)(strlen(sec)/2); ++i) {
	  css=sec[i];
	  sec[i]=sec[strlen(sec)-i-1];
	  sec[strlen(sec)-i-1]=css;
	}
    }
    free(line);

    return 0;
}
