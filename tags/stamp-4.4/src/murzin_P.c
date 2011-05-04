#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* returns a Murzin type P-value given a STAMP alignment of
 * See A.G. Murzin "Sweet tasting protein monellin is related to the cystatin family of 
 *  thiol proteinase inhibitors" J. Mol. Biol., 230, 689-694, 1993
 *
 * n  is the number of structurally equivalent positions
 * m  is the number of identical positions
 *
 * the routine will return -1.0 if an error is encountered */

double Fct(int N) {
	int i;
	double F;

	if(N>69) {
		F = 1e99;
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

	sigma = sqrt(n * p * (1-p));
	m_o   = n*p;

    if(m<=(m_o+sigma)) { /* Test for validity of P calculation */
	    Pm = 1.0;
	} else {
	   t1 = Fct(n)/(Fct(m)*Fct(n-m));
	   t2 = (double)pow((double)p,(double)m)*(double)pow((double)(1-p),(double)(n-m));
	   Pm = t1 * t2;
	}
/*	printf("\n MURZIN_P sigma = %f m_o = %f P(p=%f) = %e * %e = %e\n",sigma,m_o,p,t1,t2,Pm);  */

	if(Pm<1e-100) { Pm=0.0; }
	return Pm;
}
