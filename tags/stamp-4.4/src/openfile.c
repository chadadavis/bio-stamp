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

/* This checks for compression, etc and runs open/popen as appropriate */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FILE *openfile(char *filename, char *type) {

	FILE *handle;
	char command[1000];

	if(strcmp(type,"r")==0) { /* A read file */
		if(strcmp(&filename[strlen(filename)-2],".Z")==0) { /* UNIX compression */
/*		     sprintf(command,"zcat %s 2> /dev/null",filename);  */
		     sprintf(command,"zcat %s",filename);  
		     handle=popen(command,"r");
		} else if(strcmp(&filename[strlen(filename)-3],".gz")==0) { /* gzipped  */
/*		     sprintf(command,"gunzip -c %s 2> /dev/null",filename);  */
		     sprintf(command,"gunzip -c %s",filename); 
		     handle=popen(command,"r");
		} else { /* presume no compression */
		    handle=fopen(filename,"r");
		}
	} else {
	  	    handle=fopen(filename,"w");
	}
	return handle;
}

	
