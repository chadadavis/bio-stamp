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

char *getrksum(IN)
FILE *IN;
/* returns the summary from a DEFINE file.  The DEFINE output is from a modified
 *  version of the Richards and Kundrot program.  Refer to RBR for details */
{
	char *buff;
	char *sec;
	int start,len;

	buff=(char*)malloc(1000*sizeof(char));
	sec=(char*)malloc(sizeof(char));
	start=len=0;

	while(fgets(buff,100,IN)!=NULL) {
	   if(start) {
	     if(buff[0]=='*') 
		break;
	     if(buff[20]!='-') 
		sec[len]='H';
	     else if(buff[22]!='-') 
		sec[len]='3';
	     else if(buff[26]!='-') 
		sec[len]='B';
	     else if(buff[24]!='-') 
		sec[len]='E';
	     else  
		sec[len]='-';
	     len++;
	     sec=(char*)realloc(sec,(len+1)*sizeof(char));
	    }
	    if(!start && buff[0]=='*') 
	       start=1;
	}
	sec[len]='\0';
	free(buff);
	return sec;
}
