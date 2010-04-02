/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current contact addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Bioinformatics                      EMBL-European Bioinformatics Institute
 SmithKline Beecham Pharmaceuticals  Wellcome Trust Genome Campus
 New Frontiers Science Park (North)  Hinxton, Cambridge, CB10 1SD U.K.
 Harlow, Essex, CM19 5AW, U.K.       
 Tel: +44 1279 622 884               Tel: +44 1223 494 414
 FAX: +44 1279 622 200               FAX: +44 1223 494 468
 e-mail: russelr1@mh.uk.sbphrd.com   e-mail geoff@ebi.ac.uk
                                     WWW: http://barton.ebi.ac.uk/

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>

char *getrksum(FILE *IN) { 
/* returns the summary from a DEFINE file.  The DEFINE output is from a modified
 *  version of the Richards and Kundrot program.  Refer to RBR for details */

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
