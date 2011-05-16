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
