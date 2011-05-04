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

int newoutput(FILE *TRANS, struct domain_loc *domain, int ndomain, int writetrans) {

	int i,j,k;
        for(i=0; i<ndomain; ++i) {
	   if(printdomain(TRANS,domain[i],writetrans)==-1) return -1;
        } 
        return 0;
}

int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans) {

	int i,j,k,l,m;
   	fprintf(TRANS,"%s ",domain.filename);
	fprintf(TRANS,"%s {",domain.id);
	for(j=0; j<domain.nobj; ++j) {
	   if(domain.start[j].cid==' ') domain.start[j].cid='_';
	   if(domain.start[j].in==' ') domain.start[j].in='_';
	   if(domain.end[j].cid==' ') domain.end[j].cid='_';
	   if(domain.end[j].in==' ') domain.end[j].in='_';

	   switch(domain.type[j]) {
	      case 1: fprintf(TRANS," ALL "); break;
	      case 2: fprintf(TRANS," CHAIN %c ",domain.start[j].cid); break;
	      case 3: fprintf(TRANS," %c %d %c to %c %d %c ",
		domain.start[j].cid,domain.start[j].n,
		domain.start[j].in,
		domain.end[j].cid,domain.end[j].n,
		domain.end[j].in); break;
	   } 
	}
	if(writetrans) {
	   fprintf(TRANS,"\n");
	   for(j=0; j<3; ++j) {
	     for(k=0; k<3; ++k) fprintf(TRANS,"%10.5f ",domain.R[j][k]); 
	     fprintf(TRANS,"     %10.5f ",domain.V[j]);
	     if(j!=2) fprintf(TRANS,"\n");
	   }
	} 
	fprintf(TRANS," }\n");

	return 0;
}
