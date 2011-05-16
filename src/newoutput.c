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
