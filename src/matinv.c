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
#include <math.h>
#define TINY 1.0e-20

/* Inverts a 3 x 3 matrix. This routine was written with reference to
 *   reading Numerical Recipies in C. 
 * W.H Press, B.P Flannery, S.A. Teukolsky and W.T. Vetterling, 
 * "Numerical Recipes in C: The Art of Scientific Computing",    
 * Cambridge University Press, 1988. */

void lubksb(float **A, int n, int *indx, float b[]);
int ludcmp(float **a, int n, int *indx, float *d);
float *vector(int nl, int nh);
void free_vector(float *v, int nl, int nh);


void matinv(float **a, float **y, float d, int *indx) {

int i,j,N;
float *col;



col=(float*)malloc(10*sizeof(float));
N=3;

if(ludcmp(a,N,indx,&d)==-1) { /* matrix is singular, so just copy a to y (ie. I^-1 = I) */
   for(j=1; j<=N; j++) { 
      for(i=1; i<=N; i++) {
		y[i][j]=a[i][j];
      }
   }
} else {
   for(j=1; j<=N; j++) {
      for(i=1; i<=N; i++) col[i]=0.0;
      col[j]=1.0;
      lubksb(a,N,indx,col);
      for(i=1; i<=N; i++) y[i][j]=col[i];
   }
}

free(col);


}

void lubksb(float **a, int n, int *indx, float b[])
{
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
		for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];

	}
}



int ludcmp(float **a, int n, int *indx, float *d)
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv,*vector();
	void free_vector();

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) {
			   big=temp;
			}
		if (big == 0.0)  return -1; /* matrix is singular */
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) {
			   a[i][j] *= dum;
		  	}
		}
	}
	free_vector(vv,1,n);
}
float *vector(int nl, int nh)
{
        float *v;

        v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
        if(!v) {
	    printf("allocation failure in vector()");
	    exit(-1);
	}
	
        return v-nl;
}


void free_vector(float *v, int nl, int nh)
{
        free((char*) (v+nl));
}
