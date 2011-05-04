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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <math.h>

#include "include.h"

/* roughfit: this writes a bloc file which consists of the sequences aligned from
 *   their N-terminal ends, and uses this to run 'ampsfit'. It then reads in the
 *   obtained transformations to proceed with STAMP.  This avoids having to
 *   run AMPS intially.  This will not always work, especially if the sequences
 *   are of significantly different lengths. 
 *
 * modification: now doesn't run ampsfit, just aligns all sequences to the
 *  first one -- avoids wonky system calls within programs */

int roughfit(domain,ndomain,parms)
struct domain_loc *domain;
int ndomain;
struct parameters *parms;
{
	int i,j,k,l,m;
	int counter,ndone;
	int ntofit;
	
	float rmsd;

	char sys[200]; 


	FILE *TRANS;
	printf("Running roughfit.\n");

	fprintf(parms[0].LOG,"\n\nROUGH FIT has been requested.\n");
	fprintf(parms[0].LOG,"  The sequences will be aligned from their N-terminal ends, and\n");
	fprintf(parms[0].LOG,"  the resulting equivalences will be used to generate an inital\n");
	fprintf(parms[0].LOG,"  superposition.\n");


        /*  We will fit all domains onto the first domain */
	for(i=1; i<ndomain; ++i) {
	   if(domain[i].ncoords>domain[0].ncoords) ntofit=domain[0].ncoords;
	   else ntofit=domain[i].ncoords;
           rmsd=matfit(domain[0].coords,domain[i].coords,domain[i].R,domain[i].V,ntofit,1,parms[0].PRECISION);
/*	   printf("Domains %s onto %s RMS of %f on %d atoms\n",domain[0].id,domain[i].id,rmsd,ntofit); */
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
	  matmult(domain[i].R,domain[i].V,domain[i].coords,domain[i].ncoords,parms[0].PRECISION);
	}

	/* and we are done */
	fprintf(parms[0].LOG,"\n");

	return 0;
}
