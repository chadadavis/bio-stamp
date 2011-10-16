#include <stdio.h>

void probplot(int **prob, int lena, int lenb, int n, int cutoff, FILE *OUTPUT) {

	int i,j;
	float Pij;

	for(i=0; i<lena; i++) {
            for(j=0; j<lenb; j++)  {
               Pij=(float)prob[i+1][j+1]/(float)n;
	       if(i==j)
		  fprintf(OUTPUT,"D");
               else if(Pij>(float)cutoff) 
                  fprintf(OUTPUT,"1");
               else 
                  fprintf(OUTPUT,"0");
            } /* End of for(j... */
          fprintf(OUTPUT,"\n");
          } /* End of for(i... */

} 

