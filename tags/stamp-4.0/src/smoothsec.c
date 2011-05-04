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

/* takes a three state SS assignment (H,B or c) and removes those H or B runs that are less than
 *  minhelixlen or minstrandlen (replacing them by c) */
int smoothsec(sec,minhelixlen,minstrandlen)
char *sec;
int minhelixlen;
int minstrandlen;
{
	int i,j,k,l;
	int len,add;

	len=strlen(sec);
	for(j=0; j<len; ++j) {
          if(sec[j]=='H') {
	     add=0;
	     /* a helix */
	     for(k=j; k<len; ++k) {
	        if(sec[k]=='H') add++;
	        else break;
   	     }
	     if(add<minhelixlen) { /* too short -- remove */
	        for(l=j; l<k; ++l) sec[l]='c';
	     }
	     j=k-1;
	   } else  if(sec[j]=='B') {
	      add=0;
	      /* a strand */
	      for(k=j; k<len; ++k) {
                  if(sec[k]=='B') add++;
                  else break;
              }    
              if(add<minstrandlen) { /* too short -- remove */ 
                 for(l=j; l<k; ++l) sec[l]='c';
              }
	      j=k-1;
           } else if(sec[j]!='c' && sec[j]!=' ') {
		fprintf(stderr,"error: unrecognised secondary structure character found: %c\n",sec[j]);
		return -1;
	   }
        }
	return 0;
}
