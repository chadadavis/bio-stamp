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

/* displays alignments with a maximum line width 
 *  (ie. splits them onto more than one line)
 *
 * assumes all strings are the same length 
 */

int display_align(seqa,na,seqb,nb,seca,secb,fit,value,cols,printsec,printdat,OUTPUT)
char **seqa;
int na;
char **seqb;
int nb;
char **seca,**secb;
char *fit;
char *value;
int cols;	/* maximum number of columns */
int printsec;	/* if true print secondary structure information */
int printdat; 	/* if true print fit and value, else do not */
FILE *OUTPUT;
{
	int i,j,k;
	int count;
	int max;

	count=0;
	max=strlen(seqa[0]);
	fprintf(OUTPUT,"\nThe alignment:\n");
	while(count<max) {
	   fprintf(OUTPUT,"Position ");
	   for(i=count; i<count+cols-2; ++i) {
	      if((i+3)%10==0) {
		 fprintf(OUTPUT,"%3d",i+3);
		 i+=2;
	      } else fprintf(OUTPUT," ");
	      if(i>=max-2) break;
	   }
	   fprintf(OUTPUT,"\n");
	   for(i=0; i<na; ++i) {
	     fprintf(OUTPUT,"A(aa%3d): ",i+1);
	     dispone(seqa[i],cols,count,max,OUTPUT);
	   }
	   if(printsec) {
	     for(i=0; i<na; ++i) {
		fprintf(OUTPUT,"A(ss%3d): ",i+1);
		dispone(seca[i],cols,count,max,OUTPUT);
	     }
	   }
	   if(printdat) {
	     fprintf(OUTPUT,"FIT:      ");
	     dispone(fit,cols,count,max,OUTPUT);
	     fprintf(OUTPUT,"VALUE:    ");
	     dispone(value,cols,count,max,OUTPUT);
	   }
	   for(i=0; i<nb; ++i) {
	     fprintf(OUTPUT,"B(aa%3d): ",i+1);
	     dispone(seqb[i],cols,count,max,OUTPUT);
	   }
	   if(printsec) {
	      for(i=0; i<nb; ++i) {
		fprintf(OUTPUT,"B(ss%3d): ",i+1);
		dispone(secb[i],cols,count,max,OUTPUT);
	     }
	   }
	   count+=cols;
	   fprintf(OUTPUT,"\n");
	}
	return 0;

}
int dispone(seq,cols,count,max,OUTPUT)
char *seq;
int cols,count,max;
FILE *OUTPUT;
{
	int j;
      	for(j=0; j<cols; ++j)
	   if((count+j)<max) fprintf(OUTPUT,"%c",seq[count+j]);
      	fprintf(OUTPUT,"\n");
	return 0;
}
