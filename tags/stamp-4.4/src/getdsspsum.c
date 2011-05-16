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
