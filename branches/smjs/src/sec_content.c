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

/* return the percent helix sheet and coil for a given secondary
 *  structure string */

int sec_content(char *sec, int npos, int type, int *pera, int *perb, int *perc) {

	int i;
	int a,b,c;
	
	a=b=c=0;
	for(i=0; i<npos; ++i) {
   	   if((type==1 && (sec[i]=='H' || sec[i]=='G')) ||
	      (type==2 && (sec[i]=='H' || sec[i]=='3')) ||
	      (type==3 && (sec[i]=='H')) ) a++; 
	   else if ((type==1 && sec[i]=='E') ||
	           ((type==2 || type==3) && sec[i]=='B') ) b++; 
	   else c++;
	}
	(*pera)=(int)(100*(float)a/(float)npos);  /* Alpha */
        (*perb)=(int)(100*(float)b/(float)npos);  /* Beta */
	(*perc)=(int)(100*(float)c/(float)npos);  /* Coil */
	return 0;
}
