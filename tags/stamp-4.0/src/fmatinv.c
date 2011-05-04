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
/* This was written with reference to
 *
 * W.H Press, B.P Flannery, S.A. Teukolsky and W.T. Vetterling,
 * "Numerical Recipes in C: The Art of Scientific Computing",
 * Cambridge University Press, 1988.
 */
#define N 3
/* Inverts a 3 x 3 matrix. This routine makes use of 
 *   subroutines ludcmp and lubskb. 'a' contains the
 *   original matrix, 'y' is the inverse.  'd' is + or - 1 for some
 *   reason or other, and 'indx' contains a sequence of row operations
 *   performed to obtain 'y'.   */

void matinv(a,y,d,indx)
float **a, **y;
float d;
int *indx;

{

int i,j;
float *col;

col=(float*)malloc(10*sizeof(float));

ludcmp(a,N,indx,&d);
for(j=1; j<=N; j++) {
   for(i=1; i<=N; i++) col[i]=0.0;
   col[j]=1.0;
   lubksb(a,N,indx,col);
   for(i=1; i<=N; i++) y[i][j]=col[i];
}
free(col);

}
