/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* TPW : Replaced obsolete values.h with float.h to allow compilation on OS X. This
 * also entails replacing references to DMAXEXP in the code below with DBL_MAX_EXP.
 * */
#include <float.h>
#include <math.h>

/* returns a Murzin type P-value given a STAMP alignment of
 * See A.G. Murzin "Sweet tasting protein monellin is related to the cystatin family of 
 *  thiol proteinase inhibitors" J. Mol. Biol., 230, 689-694, 1993
 *
 * n  is the number of structurally equivalent positions
 * m  is the number of identical positions
 *
 * the routine will return -1.0 if an error is encountered */

double bico(int n, int k);
double factln(int n);
double gammln(double xx);
double bnldev(float pp, int n, long *idum);


/* SMJS This routine is not used anymore */
double Fct(int N)
{
    int i;
    double F;

    if (N > 170) {
	F = DBL_MAX; /* 1.7976931348623157E+308; */
    } else {
	F = 1.0;
	for (i = 1; i <= N; ++i) {
	    F *= (double) i;
	}
    }
    return F;
}

double murzin_P(int n, int m, double p)
{


    double Pm;
    double sigma;
    double m_o;

    double t1, t2;

/*
        printf("Entered murzin_P\n");
        fflush(stdout);
*/
    sigma = sqrt(n * p * (1 - p));
    m_o = n * p;

    if (m <= (m_o + sigma)) {	/* Test for validity of P calculation */
	Pm = 1.0;
/* SMJS Changed to use numerical recipes routines */
    } else {
	t1 = bico(n, m);
	t2 = (double) pow((double) p,
			  (double) m) * (double) pow((double) (1 - p),
						     (double) (n - m));
	if (t1 > -0.99) {
	    Pm = t1 * t2;
	} else {
	    Pm = -1.0;
	}
    }
/*
printf("\n MURZIN_P sigma = %f m_o = %f P(p=%f) = %e * %e = %e\n",sigma,m_o,p,t1,t2,Pm);    
*/

/* SMJS Added Pm>0.0 */
    if (Pm < 1e-100 && Pm > 0.00) {
	Pm = 0.0;
    }
    return Pm;
}

/* This routine is adapted from one in Numerical Recipes in C */
double bico(int n, int k)
{
/* TPW : DBL_MAX_EXP replaces DMAXEXP owing to values.h being superseded by float.h */
#ifndef LN_MAXDOUBLE
#define LN_MAXDOUBLE (M_LN2 * DBL_MAX_EXP)
#endif
    double resfact;

    resfact = factln(n) - factln(k) - factln(n - k);

    if (resfact > LN_MAXDOUBLE)
	return -1.0;
    else
	return floor(0.5 + exp(resfact));
}

/* This routine is adapted from one in Numerical Recipes in C */
double factln(int n)
{
    static double a[101];

    if (n <= 100)
	return a[n] ? a[n] : (a[n] = gammln(n + 1.0));
    else
	return gammln(n + 1.0);
}

/* This routine is adapted from one in Numerical Recipes in C */
/* Note the Numerical Recipes first edition code was wrong. */
/* It was corrected in the second edition */
double gammln(double xx)
{
    double x, y, tmp, ser;
    static double cof[6] =
	{ 76.18009172947146, -86.50532032941677, 24.01409824083091,
	-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
    };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;

    for (j = 0; j <= 5; j++) {
	ser += cof[j] / ++y;
    }
    return -tmp + log(2.5066282746310005 * ser / x);
}
