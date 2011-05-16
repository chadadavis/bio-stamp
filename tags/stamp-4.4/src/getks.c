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
#include <stdlib.h>

#include "stamp.h"

/* This routine looks for a DSSP file, if found it reads in the DSSP summary
 *  and stores it in domain[i].sec for future use */

int getks(struct domain_loc *domain, int ndomain, struct parameters *parms) {

	FILE *DSSP;
	int i,j,k;
	int include,count;
	int total,add;
	int retval;
	char c,chain;
	char *filename;
	char *label;
	char *empty = NULL; 
	
	count=0;
	retval=0;


	for(i=0; i<ndomain; ++i) {
	   fprintf(parms[0].LOG,"%s -- ",domain[i].id);
	   /* first check to see if there is a DSSP file that uses the whole ID name */
           filename=getfile(domain[i].id,parms[0].dsspfile,strlen(domain[i].id),parms[0].LOG);
	   /* if there is, use it, otherwise try to get one using the four letter code */
	   if(filename[0]=='\0') filename=getfile(domain[i].id,parms[0].dsspfile,4,parms[0].LOG);
           if(filename[0]=='\0') {
               free(filename);
               filename=getfile(domain[i].id,parms[0].dsspfile,4,parms[0].LOG);
           }
	   if(filename[0]=='\0') {
	      fprintf(parms[0].LOG," no DSSP file found for %s\n",domain[i].id);
	      for(j=0; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
	      domain[i].sec[j]='\0';
	      retval=-1; /* if any of the sequences have missing secondary structures */
	   } else {
	      DSSP=openfile(filename,"r");
	      total=0;
	      fprintf(parms[0].LOG," using file %s\n",filename);
	      for(j=0; j<domain[i].nobj; ++j) {
	         if(get_dssp_sum(DSSP,domain[i].start[j],domain[i].end[j],
		    domain[i].type[j],&domain[i].sec[total],domain[i].reverse[j],
		    (parms[0].MAX_SEQ_LEN-total),&add,parms[0].LOG)==-1) 
		    retval=-1;
	         total+=add;
		 closefile(DSSP,filename);
	         DSSP=openfile(filename,"r");
	      }
	      closefile(DSSP,filename);
	      free(filename);
	      if(total!=domain[i].ncoords) {
		 fprintf(parms[0].LOG,"warning: DSSP summary found was incomplete -- the results may have errors\n");
		 for(j=total; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
		 domain[i].sec[j]='\0';
	      }
	      display_align(&domain[i].aa,1,&domain[i].sec,1,&domain[i].aa,&domain[i].sec,empty,empty,parms[0].COLUMNS,0,0,parms[0].LOG);
	   }
	}
	return retval;
	      
} 
	  


