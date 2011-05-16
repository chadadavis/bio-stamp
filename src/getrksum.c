/*
Copyright (1997,1998,1999,2010) Robert B Russell & Geoffrey J Barton

This file is part of STAMP

STAMP is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE  See the
GNU General Public License for more details A copy of the license
can be found in the LICENSE file in the STAMP installation directory

STAMP was developed by Robert B Russell and Geoffrey J Barton of
current addresses:

 Prof Robert B Russell (RBR)                      Prof. Geoffrey J. Barton (GJB)
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

 RB Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992)
*/
#include <stdio.h>

char *getrksum(FILE *IN) { 
/* returns the summary from a DEFINE file  The DEFINE output is from a modified
 *  version of the Richards and Kundrot program  Refer to RBR for details */

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
