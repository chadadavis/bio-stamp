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
#include <string.h>
#include <stdlib.h>

/* displays alignments with a maximum line width 
 *  (ie. splits them onto more than one line)
 *
 * assumes all strings are the same length 
 */

int dispone(char *seq, int cols, int count, int max, FILE *OUTPUT);

int display_align(char **seqa, int na, char **seqb, int nb, 
	          char **seca, char **secb, char *fit, char *value, 
                  int cols, int printsec, int printdat, FILE *OUTPUT)
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
int dispone(char *seq, int cols, int count, int max, FILE *OUTPUT)
{
	int j;
      	for(j=0; j<cols; ++j)
	   if((count+j)<max) fprintf(OUTPUT,"%c",seq[count+j]);
      	fprintf(OUTPUT,"\n");
	return 0;
}
