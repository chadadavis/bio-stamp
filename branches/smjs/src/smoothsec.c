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

/* takes a three state SS assignment (H,B or c) and removes those H or B runs that are less than
 *  minhelixlen or minstrandlen (replacing them by c) */
int smoothsec(char *sec, int minhelixlen, int minstrandlen) {

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
