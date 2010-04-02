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
#include "stamp.h"

int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans) {

	int i,j,k,l,m;
   	fprintf(TRANS,"%s ",domain.filename);
	fprintf(TRANS,"%s {",domain.id);
	for(j=0; j<domain.nobj; ++j) {
	   if(domain.start[j].cid==' ') domain.start[j].cid='_';
	   if(domain.start[j].in==' ') domain.start[j].in='_';
	   if(domain.end[j].cid==' ') domain.end[j].cid='_';
	   if(domain.end[j].in==' ') domain.end[j].in='_';

	   if(domain.reverse[j]==1) fprintf(TRANS," REVERSE ");
	   switch(domain.type[j]) {
	      case 1: fprintf(TRANS," ALL "); break;
	      case 2: fprintf(TRANS," CHAIN %c ",domain.start[j].cid); break;
	      case 3: fprintf(TRANS," %c %d %c to %c %d %c ",
		domain.start[j].cid,domain.start[j].n,
		domain.start[j].in,
		domain.end[j].cid,domain.end[j].n,
		domain.end[j].in); break;
	   } /* end of switch... */
	   if(domain.start[j].cid=='_') domain.start[j].cid=' ';
	   if(domain.start[j].in=='_') domain.start[j].in=' ';
	   if(domain.end[j].cid=='_') domain.end[j].cid=' ';
	   if(domain.end[j].in=='_') domain.end[j].in=' ';
	}
	if(writetrans) {
	   fprintf(TRANS,"\n");
	   for(j=0; j<3; ++j) {
	     for(k=0; k<3; ++k) fprintf(TRANS,"%10.5f ",domain.R[j][k]); 
	     fprintf(TRANS,"     %10.5f ",domain.V[j]);
	     if(j!=2) fprintf(TRANS," \n");
	   }
	} /* end of if(write... */
	fprintf(TRANS," }\n");

	return 0;
}
