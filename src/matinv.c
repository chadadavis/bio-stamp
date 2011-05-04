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
#include <malloc.h>
#include <math.h>
#define N 3
#define TINY 1.0e-20

/* Inverts a 3 x 3 matrix. This routine was written with reference to
 *   reading Numerical Recipies in C. 
 * W.H Press, B.P Flannery, S.A. Teukolsky and W.T. Vetterling, 
 * "Numerical Recipes in C: The Art of Scientific Computing",    
 * Cambridge University Press, 1988. */

int ludcmp();
void lubksb();
float *vector();
void free_vector();

void matinv(a,y,d,indx)
float **a, **y;
float d;
int *indx;

{

int i,j;
float *col;

col=(float*)malloc(10*sizeof(float));

if(ludcmp(a,N,indx,&d)==-1) { /* matrix is singular, so just copy a to y (ie. I^-1 = I) */
   for(j=1; j<=N; j++) 
      for(i=1; i<=N; i++) y[i][j]=a[i][j];
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

void lubksb(a,n,indx,b)
float **a,b[];
int n,*indx;
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



int ludcmp(a,n,indx,d)
int n,*indx;
float **a,*d;
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
float *vector(nl,nh)
int nl,nh;
{
        float *v;

        v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
        if(!v) {
	    printf("allocation failure in vector()");
	    exit(-1);
	}
	
        return v-nl;
}


void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
        free((char*) (v+nl));
}
