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

int probcalc(atoms1,lena,atoms2,lenb,prob,parms)
/* Calculates the probability matrix, prob, after Rossmann and Argos, given
 * the two sets of coordinates, (float)atoms1 and (float)atoms2. 
 *
 * Modification January 1995, Old shite version was inefficient.  Now
 *  fewer redundent boolean tests are done, and most importantly, the
 *  matrix is only navigated twice if absolutely necessary */
 
int **atoms1, **atoms2;
int **prob;
int lena,lenb;
struct parameters *parms;
{

        int i,j,k,ii,jj;
	float sum,sumsq;
        float Dij,Cij,mean,sd,const1,const2;
	int start,end,ll;


	const1=parms[0].const1;
	const2=parms[0].const2;



	if(parms[0].BOOLEAN) {
	  for(i=0; i<lena; i++) {
           ii=i+1;
           for(j=0; j<lenb; j++)  {
              jj=j+1;
              /* The following calculates a Probability matrix after Rossmann and
               *  Argos (J.Mol.Biol., 105_, 75 (1976))...
               * The routine 'rossmann' returns both the probabilty Pij, and the
               *  pure distance parameter Di */
	       prob[ii][jj]=(rossmann(&atoms1[i],&atoms2[j], (i==0 || j==0),(i==lena-1 || j==lenb-1),
                                  const1,const2,&Dij,&Cij,parms[0].PRECISION)>=parms[0].BOOLCUT);
	       } 
          }  
	} else if(!parms[0].STATS) {
	 /* using fixed mean and sd, don't need to calculate mean or standard deviation */
         mean=parms[0].NMEAN;
         sd=parms[0].NSD;
	 for(i=0; i<lena; i++) {
           ii=i+1;
           for(j=0; j<lenb; j++)  {
              jj=j+1;
              /* The following calculates a Probability matrix after Rossmann and
               *  Argos (J.Mol.Biol., 105_, 75 (1976))...
               * The routine 'rossmann' returns both the probabilty Pij, and the
               *  pure distance parameter Dij.  */
               prob[ii][jj]=(int)
                  ((float)parms[0].PRECISION*(rossmann(&atoms1[i],&atoms2[j],
                           (i==0 || j==0),(i==lena-1 || j==lenb-1),
                           const1,const2,&Dij,&Cij,parms[0].PRECISION) - mean)/sd);
	   }      
         }
	} else {
          sum=sumsq=0.0;
          for(i=0; i<lena; i++) {
	   ii=i+1;
           for(j=0; j<lenb; j++)  {
	      jj=j+1;
	      /* The following calculates a Probability matrix after Rossmann and
	       *  Argos (J.Mol.Biol., 105_, 75 (1976))...
	       * The routine 'rossmann' returns both the probabilty Pij, and the
	       *  pure distance parameter Dij.  */
               prob[ii][jj]=(int)((float)parms[0].PRECISION*rossmann(&atoms1[i],&atoms2[j],
			   (i==0 || j==0),(i==lena-1 || j==lenb-1),
			   const1,const2,&Dij,&Cij,parms[0].PRECISION));

	       sum+=(float)prob[ii][jj]; 
               sumsq+=(float)(prob[ii][jj]*prob[ii][jj]);
           }
         }  
	 mean=((float)parms[0].PRECISION*(sum/(float)(lena*lenb)));
	 sd=(float)parms[0].PRECISION*
	      (float)sqrt( (sumsq-(sum*sum)/(float)(lena*lenb)) / (lena*lenb-1) );
	 /* Now we must find out how many SD's above the mean each value
	  *  in the probability matrix is. */
	 for(i=0; i<lena; i++) {
	     ii=i+1;
     	     for(j=0; j<lenb; ++j) {
		jj=j+1;
                prob[ii][jj]=(int)( (float)parms[0].PRECISION*((float)(prob[ii][jj]-mean)/(float)(sd)));
	     }
	   }
	 }

	return 0;

} 
