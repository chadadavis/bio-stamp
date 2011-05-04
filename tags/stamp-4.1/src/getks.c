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
#include <stdlib.h>

#include <stamp.h>

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
	char *empty;
	
	count=0;
	retval=0;


	for(i=0; i<ndomain; ++i) {
	   fprintf(parms[0].LOG,"%s -- ",domain[i].id);
	   /* first check to see if there is a DSSP file that uses the whole ID name */
           filename=getfile(domain[i].id,parms[0].dsspfile,strlen(domain[i].id),parms[0].LOG);
	   /* if there is, use it, otherwise try to get one using the four letter code */
	   if(filename[0]=='\0') filename=getfile(domain[i].id,parms[0].dsspfile,4,parms[0].LOG);
	   if(filename[0]=='\0') {
	      fprintf(parms[0].LOG," no DSSP file found for %s\n",domain[i].id);
	      for(j=0; j<domain[i].ncoords; ++j) domain[i].sec[j]='?';
	      domain[i].sec[j]='\0';
	      retval=-1; /* if any of the sequences have missing secondary structures */
	   } else {
	      DSSP=fopen(filename,"r");
	      total=0;
	      fprintf(parms[0].LOG," using file %s\n",filename);
	      for(j=0; j<domain[i].nobj; ++j) {
	         if(get_dssp_sum(DSSP,domain[i].start[j],domain[i].end[j],
		    domain[i].type[j],&domain[i].sec[total],domain[i].reverse[j],
		    (parms[0].MAX_SEQ_LEN-total),&add,parms[0].LOG)==-1) 
		    retval=-1;
	         total+=add;
		 rewind(DSSP);
	      }
	      fclose(DSSP);
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
	  


