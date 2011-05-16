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
