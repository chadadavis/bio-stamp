/*
Copyright (1997,1998,1999,2010) Robert B. Russell & Geoffrey J. Barton

This file is part of STAMP.

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. A copy of the license
can be found in the LICENSE file in the STAMP installation directory.

STAMP was developed by Robert B. Russell and Geoffrey J. Barton of
current addresses:

 Prof. Robert B. Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
 Cell Networks, University of Heidelberg            College of Life Sciences
 Room 564, Bioquant                                 University of Dundee
 Im Neuenheimer Feld 267                            Dow Street
 69120 Heidelberg                                   Dundee DD1 5EH
 Germany                                            UK
                                                
 Tel: +49 6221 54 513 62                            Tel: +44 1382 385860
 Fax: +49 6221 54 514 86                            FAX: +44 1382 385764
 Email: robert.russell@bioquant.uni-heidelberg.de   E-mail g.j.barton@dundee.ac.uk
 WWW: http://www.russell.embl-heidelberg.de         WWW: http://www.compbio.dundee.ac.uk

 All use of STAMP must cite: 

 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*/
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
