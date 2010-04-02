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
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <math.h>

/* returns a Murzin type P-value given a STAMP alignment of
 * See A.G. Murzin "Sweet tasting protein monellin is related to the cystatin family of 
 *  thiol proteinase inhibitors" J. Mol. Biol., 230, 689-694, 1993
 *
 * n  is the number of structurally equivalent positions
 * m  is the number of identical positions
 *
 * the routine will return -1.0 if an error is encountered */

double bico(int n,int k);
double factln(int n);
double gammln(double xx);
double bnldev(float pp, int n, long *idum);


/* SMJS This routine is not used anymore */
double Fct(int N) {
	int i;
	double F;

	if(N>170) {
		F = 1.7976931348623157E+308;
	} else {
		F = 1.0;
		for(i=1; i<=N; ++i) {
			F*=(double)i;
		}
	}
	return F;
}

double murzin_P(int n,int m,double p) {


	double Pm;
	double sigma;
	double m_o;

	double t1,t2;

/*
        printf("Entered murzin_P\n");
        fflush(stdout);
*/
	sigma = sqrt(n * p * (1-p));
	m_o   = n*p;

        if(m<=(m_o+sigma)) { /* Test for validity of P calculation */
	    Pm = 1.0;
/* SMJS Changed to use numerical recipes routines */
	} else {
           t1 = bico(n,m);
	   t2 = (double)pow((double)p,(double)m)*(double)pow((double)(1-p),(double)(n-m));
           if (t1 > -0.99)
           {
              Pm = t1*t2;
	   } 
           else 
           {
              Pm = -1.0;
           }
        }
/*
printf("\n MURZIN_P sigma = %f m_o = %f P(p=%f) = %e * %e = %e\n",sigma,m_o,p,t1,t2,Pm);    
*/

/* SMJS Added Pm>0.0 */
	if(Pm<1e-100 && Pm>0.00) { Pm=0.0; } 
	return Pm;
}

/* This routine is adapted from one in Numerical Recipes in C */
double bico(int n,int k)
{
#ifndef LN_MAXDOUBLE   
#define LN_MAXDOUBLE (M_LN2 * DMAXEXP)
#endif
   double resfact;
   
   resfact = factln(n)-factln(k)-factln(n-k);

   if (resfact > LN_MAXDOUBLE) return -1.0;
   else return floor(0.5+exp(resfact));
}

/* This routine is adapted from one in Numerical Recipes in C */
double factln(int n)
{
   static double a[101];

   if (n<=100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
   else return gammln(n+1.0);
}

/* This routine is adapted from one in Numerical Recipes in C */
/* Note the Numerical Recipes first edition code was wrong. */
/* It was corrected in the second edition */
double gammln(double xx)
{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
                         -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   int j;
  
   y = x = xx;
   tmp  = x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser  = 1.000000000190015;
   
   for (j=0;j<=5;j++)
   {
      ser += cof[j]/++y;
   }
   return -tmp+log(2.5066282746310005*ser/x);
}

