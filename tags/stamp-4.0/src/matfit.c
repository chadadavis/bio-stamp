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
/**********************************************************
      MATFIT

      Note that XA is ATOMS1 and XB is ATOMS2   G.J.B.

      SUBROUTINE TO FIT THE COORD SET ATOMS1(3,N) TO THE SET ATOMS2(3,N)
      IN THE SENSE OF:
             XA= R*XB +V
      R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX
      AND V IS THE OFFSET VECTOR. THIS IS AN EXACT SOLUTION

     THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S
     TECHNIQUES. SEE
      KABSCH, W. ACTA CRYST A34, 827,1978
      MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978
      WRITTEN BY S.J. REMINGTON 11/78.

      THIS SUBROUTINE USES THE IBM SSP EIGENVALUE ROUTINE 'EIGEN'

     N.B. CHANGED BY R.B.R (OCTOBER 1990) FOR USE IN THE PROGRAM
      STRUCALIGN.  SEE STRUCALIGN.F FOR DETAILS

     CHANGED BY R.B.R. (JANUARY 1991) FOR USE IN THE PROGRAM
        RIGIDBODY.
 -----------------------------------------------------------
     Translated into C by RBR in January 1991
*************************************************************/
#include <stdio.h>
#include <math.h>

/*  This function mimics the first part of the FORTRAN routine
 *    'matfit.f'.  It makes use of the FORTRAN routine 'qkfit.f'.
 *     This keeps complications involving passing great heaps
 *     of memory to FORTRAN from arising.  The resulting transformation
 *     is in the sense of:
 *	   atoms1 = R*atoms2 + V     
 *
 *
 * Special version for integer atoms */

float matfit(atoms1,atoms2,R,V,nats,entry,PRECISION)
int **atoms1,**atoms2; 
float **R, *V;
int nats, entry;
int PRECISION;
{
	/* Note that everything that is passed to the routine
	 *  qkfit must be double precision. */
	double cma[3], cmb[3];
	double umat[3][3];
	double xasq, xbsq, xni,t,rtsum,rmse;
	double r[3][3];
	int i,j,k,xn;

	xn=nats;
	xasq=0.0;
	xbsq=0.0;
	xni=1.0/((double)xn);

/* Accumulate uncorrected (for c.m.) sums and squares */

	for(i=0; i<3; ++i) {
	   cma[i]=cmb[i]=0.0;
	   for(j=0; j<=2; ++j) umat[j][i]=0.0;

           for(j=1; j<nats; ++j) {
	      for(k=0; k<3; ++k) {
	      umat[k][i]+=(double)atoms1[j][i]/(double)PRECISION*(double)atoms2[j][k]/(double)PRECISION;
	      }  /* End of for(k=... */
/*	      printf("Matfit...\n");
	      for(k=0; k<3; ++k) printf("%10.5f ",atoms1[j][k]);
	      printf("and ");
	      for(k=0; k<3; ++k) printf("%10.5f ",atoms2[j][k]); */
		
	      t=(double)atoms1[j][i]/(double)PRECISION;
	      cma[i]+=t;
	      xasq+=t*t;
	      t=(double)atoms2[j][i]/(double)PRECISION;
	      cmb[i]+=t;
	      xbsq+=t*t;

	   } /*  End of for(j=... */
	}  /* End of for(i=... */

/* Subtract cm offsets */

	for(i=0; i<=2; ++i) {
	   xasq-=cma[i]*cma[i]*xni;
	   xbsq-=cmb[i]*cmb[i]*xni;
	   for(j=0; j<=2; ++j) 
	       umat[j][i]=(umat[j][i]-cma[i]*cmb[j]*xni)*xni;
	} /* End of for(i=... */

/* Fit it */

	qkfit(&umat[0][0],&rtsum,&r[0][0],&entry);
	rmse=(xasq+xbsq)*xni - 2.0*rtsum;
	if(rmse<0.0) rmse=0.0;
	rmse=sqrt(rmse);

/* The matrix obtained from the FORTRAN routine (r) must be transfered to
 *  the 'malloc'ed C matrix (R).  */
	for(i=0; i<=2; ++i)
	   for(j=0; j<=2; ++j) R[i][j]=(float)r[j][i];


/* Calculate offset if entry=1 */

	if(!entry) return rmse;
	for(i=0; i<=2; ++i) {
	   t=0.0;
	   for(j=0; j<=2; ++j)  t+=((double)R[i][j])*cmb[j];
	   V[i]=(float)(cma[i]-t)*xni;
	} /* End of for(i=... */

	return rmse;

} /* End of function */
 
