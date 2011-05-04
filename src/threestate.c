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
#include <stdlib.h>
#include <string.h>

/* threestate -- converts secondary structures into three state 
 *  you must give the routine four strings:
 *
 * sec -- the secondary structure assignment 
 * helix -- a string of helical definitions (eg. DSSP "HGI")
 * extended -- a string of extended definitions (eg. DSSP "EB")
 * coil     -- a string of coil definitions (eg. DSSP "-TS") 
 *
 * The string supplied is altered to contain only H, B and c
 *  according to the supplied definitions 
 *
 * It will leave un-recognised characters (ie. spaces) unaltered 
 *
 */

int threestate(char *sec,char *helix,char *extended,char *coil) {

	int i,j;
	int slen,hlen,elen,clen;

	slen=strlen(sec);
	hlen=strlen(helix);
	elen=strlen(extended);
	clen=strlen(coil);

	for(i=0; i<strlen(sec); ++i) {
	     /* coil is done first, since this will keep changed 'H' or 'B' from being coil */
	     for(j=0; j<clen; ++j) {
		if(sec[i]==coil[j]) {
		   sec[i]='c';
		   break;
	   	}
	     }
	     for(j=0; j<hlen; ++j) {
		if(sec[i]==helix[j]) {
		   sec[i]='H'; 
		   break;
		}
	     }
	     for(j=0; j<elen; ++j) {
		if(sec[i]==extended[j]) {
		   sec[i]='B';
		   break;
		}
	     }
	}

	return 0;
}
