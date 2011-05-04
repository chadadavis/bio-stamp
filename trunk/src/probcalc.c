#include <stdio.h>
#include <math.h>

#include "stamp.h"
float rossmanni(int **atoms1, int **atoms2,
        const float const1, const float const2, float * const Dij, float * const Cij);

/* Calculates the probability matrix, prob, after Rossmann and Argos, given
 * the two sets of coordinates, (float)atoms1 and (float)atoms2. 
 *
 * Modification January 1995, Old shite version was inefficient.  Now
 *  fewer redundent boolean tests are done, and most importantly, the
 *  matrix is only navigated twice if absolutely necessary */

int probcalc(int **atoms1, int **atoms2, int **prob, int lena, int lenb,
	struct parameters *parms) {

        int i,j,k,ii,jj;
	float sum,sumsq;
        float Dij,Cij,mean,sd,const1,const2;
	int start,end,ll;

/* SMJS Added to hopefully speed up */
        float *DijP = &Dij;
        float *CijP = &Cij;
        int  **atom1P;
        int  ***atom2Ps;
        int    precision=parms[0].PRECISION;
        float  flprec=(float)parms[0].PRECISION;
        float  flprec2i=1.0/(float)(parms[0].PRECISION*parms[0].PRECISION);
        float  boolcut=parms[0].BOOLCUT;
        int    probij;
        float  sdi;
        int    istflag;
        int    iendflag;
        


/* SMJS Changed to use inverse */
	const1=(1.0/parms[0].const1)*flprec2i;
	const2=(1.0/parms[0].const2)*flprec2i;

/*
        if ((atom2Ps = (int ***)malloc(lenb*sizeof(int **)))==NULL)
        {
           exit(-1);
        }
        for (i=0;i<lenb;i++)
        {
           atom2Ps[i]=&(atoms2[i]);
        }
*/


	if(parms[0].BOOLEAN) {
	  for(i=0; i<lena; i++) {
           ii=i+1;
/* SMJS Added */
           atom1P = (&atoms1[i]);
/* SMJS End added */
           for(j=0; j<lenb; j++)  {
              jj=j+1;
              /* The following calculates a Probability matrix after Rossmann and
               *  Argos (J.Mol.Biol., 105_, 75 (1976))...
               * The routine 'rossmann' returns both the probabilty Pij, and the
               *  pure distance parameter Di */
/* SMJS
	       prob[ii][jj]=(rossmann(&atoms1[i],&atoms2[j], (i==0 || j==0),(i==lena-1 || j==lenb-1),
                                  const1,const2,&Dij,&Cij,parms[0].PRECISION)>=parms[0].BOOLCUT);
*/
	       prob[ii][jj]=(rossmann(atom1P,&atoms2[j], (!i || !j),(ii==lena || jj==lenb),
                                  const1,const2,DijP,CijP)>=boolcut);
	       } 
          }  
	} else if(!parms[0].STATS) {
	 /* using fixed mean and sd, don't need to calculate mean or standard deviation */
         mean=parms[0].NMEAN;
         sd=parms[0].NSD;
         sdi=1.0/parms[0].NSD;
/* SMJS Added speed optimisation by calculating edges of prob array using standard rossmann */
/*      but using new rossmanni() which does no edge checks for internal elements */
/*      Unfortunately this is a bit messy but it does give quite a good speedup */
	 for(j=0; j<lenb; j++) {
           jj=j+1;
           prob[1][jj]=(int)(flprec*(rossmann(&atoms1[0],&atoms2[j],1,(jj==lenb),
			                          const1,const2,DijP,CijP) - mean)*sdi);
           prob[lena][jj]=(int)(flprec*(rossmann(&atoms1[lena-1],&atoms2[j],(!j),1,
			                          const1,const2,DijP,CijP) - mean)*sdi);
         }
	 for(i=1; i<lena-1; i++) {
           ii=i+1;
/* SMJS Added */
           prob[ii][1]=(int)(flprec*(rossmann(&atoms1[i],&atoms2[0],1,0,
	 		                      const1,const2,DijP,CijP) - mean)*sdi);
           atom1P = (&atoms1[i]);
/* SMJS End added */
           for(j=1; j<lenb-1; j++)  {
              jj=j+1;
              /* The following calculates a Probability matrix after Rossmann and
               *  Argos (J.Mol.Biol., 105_, 75 (1976))...
               * The routine 'rossmann' returns both the probabilty Pij, and the
               *  pure distance parameter Dij.  */
/* SMJS
               prob[ii][jj]=(int)
                  ((float)parms[0].PRECISION*(rossmann(&atoms1[i],&atoms2[j],
                           (i==0 || j==0),(i==lena-1 || j==lenb-1),
                           const1,const2,&Dij,&Cij,parms[0].PRECISION) - mean)/sd);
*/
               prob[ii][jj]=(int)(flprec*(rossmanni(atom1P,&atoms2[j],
                                          const1,const2,DijP,CijP) - mean)*sdi);
	   }      
           prob[ii][lenb]=(int)(flprec*(rossmann(&atoms1[i],&atoms2[lenb-1],0,1,
			                          const1,const2,DijP,CijP) - mean)*sdi);
         }
	} else {
          sum=sumsq=0.0;
          for(i=0; i<lena; i++) {
	   ii=i+1;
/* SMJS Added */
           atom1P = (&atoms1[i]);
/* SMJS End added */
           for(j=0; j<lenb; j++)  {
	      jj=j+1;
	      /* The following calculates a Probability matrix after Rossmann and
	       *  Argos (J.Mol.Biol., 105_, 75 (1976))...
	       * The routine 'rossmann' returns both the probabilty Pij, and the
	       *  pure distance parameter Dij.  */
/* SMJS
               prob[ii][jj]=(int)((float)parms[0].PRECISION*rossmann(&atoms1[i],&atoms2[j],
			   (i==0 || j==0),(i==lena-1 || j==lenb-1),
			   const1,const2,&Dij,&Cij,parms[0].PRECISION));
*/
               probij=prob[ii][jj]=(int)(flprec*rossmann(atom1P,&atoms2[j],
                                                  (!i || !j),(ii==lena || jj==lenb),
			                          const1,const2,DijP,CijP));

	       sum+=(float)probij; 
               sumsq+=(float)(probij*probij);
           }
         }  

	 mean=(flprec*(sum/(float)(lena*lenb)));
	 sd=flprec*(float)sqrt( (sumsq-(sum*sum)/(float)(lena*lenb)) / (lena*lenb-1) );
         sdi=1.0/sd;
	 /* Now we must find out how many SD's above the mean each value
	  *  in the probability matrix is. */
	 for(i=0; i<lena; i++) {
	     ii=i+1;
     	     for(j=0; j<lenb; ++j) {
		jj=j+1;
                prob[ii][jj]=(int)( flprec*((float)(prob[ii][jj]-mean)*sdi));
	     }
	   }
	 }

/*
        free(atom2Ps);
*/

	return 0;

} 
