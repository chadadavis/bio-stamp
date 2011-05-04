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

/* return the percent helix sheet and coil for a given secondary
 *  structure string */

int sec_content(sec,npos,type,pera,perb,perc)
char *sec;
int npos;
int type;
int *pera,*perb,*perc;
{
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
