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

void probplot(prob,lena,lenb,n,cutoff,OUTPUT)
int **prob;
int lena, lenb, n;
int cutoff;
FILE *OUTPUT;
{
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

} /* End of probplot. */

