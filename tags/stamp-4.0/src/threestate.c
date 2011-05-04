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

int threestate(sec,helix,extended,coil)
char *sec,*helix,*extended,*coil;
{

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
