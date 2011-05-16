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
#include <math.h>

#include "stamp.h"

/* roughfit: this writes a bloc file which consists of the sequences aligned from
 *   their N-terminal ends, and uses this to run 'ampsfit'. It then reads in the
 *   obtained transformations to proceed with STAMP.  This avoids having to
 *   run AMPS intially.  This will not always work, especially if the sequences
 *   are of significantly different lengths. 
 *
 * modification: now doesn't run ampsfit, just aligns all sequences to the
 *  first one -- avoids wonky system calls within programs */

int roughfit(struct domain_loc *domain, int ndomain, struct parameters *parms) {

	int i,j,k,l,m;
	int counter,ndone;
	int ntofit;
	
	float rmsd;

	char sys[200]; 

	struct brookn tmps,tmpe;

	FILE *ROUGHOUT;
	

	printf("Running roughfit.\n");

	fprintf(parms[0].LOG,"\n\nROUGH FIT has been requested.\n");
	fprintf(parms[0].LOG,"  The sequences will be aligned from their N-terminal ends, and\n");
	fprintf(parms[0].LOG,"  the resulting equivalences will be used to generate an inital\n");
	fprintf(parms[0].LOG,"  superposition.\n");

	if(parms[0].roughout == 1 ) { /* output the transformations */
            fprintf(parms[0].LOG,"\nROUGH FIT transformations will be output to the file %s\n",parms[0].roughoutfile);
	    if((ROUGHOUT=fopen(parms[0].roughoutfile,"w"))==NULL) {
		fprintf(stderr,"Error opening file %s for writing\n",parms[0].roughoutfile);
	        exit(-1);
	    }
	    fprintf(ROUGHOUT,"%% Output from STAMP ROUGH FIT routine\n");
	    fprintf(ROUGHOUT,"%%  The sequences from the file %s have been aligned from their\n",parms[0].listfile);
	    fprintf(ROUGHOUT,"%%  N-terminal ends, andthe resulting equivalences were be used \n");
	    fprintf(ROUGHOUT,"%% to generate the superpositions given below.\n");
	}

        /*  We will fit all domains onto the first domain */
	for(i=1; i<ndomain; ++i) {
	   if(domain[i].ncoords>domain[0].ncoords) ntofit=domain[0].ncoords;
	   else ntofit=domain[i].ncoords;
           rmsd=matfit(domain[0].coords,domain[i].coords,domain[i].R,domain[i].V,ntofit,1,parms[0].PRECISION); 
	   fprintf(parms[0].LOG,"Domains %s onto %s RMS of %f on %d atoms\n",domain[0].id,domain[i].id,rmsd,ntofit); 
	   if(parms[0].roughout == 1) {
		fprintf(ROUGHOUT,"%% Domains %s onto %s RMS of %f on %d atoms\n",domain[0].id,domain[i].id,rmsd,ntofit); 
	   }
        }

	for(i=0; i<ndomain; ++i) {
	  fprintf(parms[0].LOG,"\nDomain %2d, %s, %d coordinates\n",i+1,domain[i].id,domain[i].ncoords);
	  fprintf(parms[0].LOG,"Applying the transformation...\n");
	  for(j=0; j<3; ++j) {
	      fprintf(parms[0].LOG,"| ");
	      for(k=0; k<3; ++k) fprintf(parms[0].LOG,"%8.5f ",domain[i].R[j][k]);
	      fprintf(parms[0].LOG," |    %8.5f\n",domain[i].V[j]);
	  }
	  fprintf(parms[0].LOG,"      ...to these coordinates.\n");
	  if(parms[0].roughout == 1) {
		fprintf(ROUGHOUT,"%s %s { ",domain[i].filename,domain[i].id);
		for(j=0; j<domain[i].nobj; ++j) {

		    if(domain[i].start[j].cid!=' ') tmps.cid=domain[i].start[j].cid;
		    else tmps.cid='_';
		    if(domain[i].end[j].cid!=' ') tmpe.cid=domain[i].start[j].cid;
                    else tmpe.cid='_';
		    if(domain[i].start[j].in!=' ') tmps.in=domain[i].start[j].in;
                    else tmps.in='_';
		    if(domain[i].end[j].in!=' ') tmpe.in=domain[i].start[j].in;
                    else tmpe.in='_';

		    if(domain[i].type[j]==1) fprintf(ROUGHOUT,"ALL");
		    else if(domain[i].type[j]==2) fprintf(ROUGHOUT,"CHAIN %c",domain[i].start[j].cid);
		    else fprintf(ROUGHOUT,"%c %d %c to %c %d %c",
			tmps.cid,domain[i].start[j].n,tmps.in,
			tmps.cid,domain[i].end[j].n,tmpe.in);
		    fprintf(ROUGHOUT," ");
		}
	        fprintf(ROUGHOUT,"\n");
		for(j=0; j<3; ++j) {
		   fprintf(ROUGHOUT,"%10.4f %10.4f %10.4f   %10.4f ",
			domain[i].R[j][0],domain[i].R[j][1],domain[i].R[j][2],domain[i].V[j]);
		   if(j==2) fprintf(ROUGHOUT," } ");
		   fprintf(ROUGHOUT,"\n");
		}
	  }
	  matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,parms[0].PRECISION);
	 
	}

	/* and we are done */
	if(parms[0].roughout == 1) fclose(ROUGHOUT);
	fprintf(parms[0].LOG,"\n");

	return 0;
}
