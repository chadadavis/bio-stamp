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
