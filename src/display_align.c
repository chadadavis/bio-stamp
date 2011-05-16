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
#include <string.h>
#include <stdlib.h>

/* displays alignments with a maximum line width 
 *  (ie splits them onto more than one line)
 *
 * assumes all strings are the same length 
 */

int dispone(char *seq, int cols, int count, int max, FILE *OUTPUT);

int display_align(char **seqa, int na, char **seqb, int nb, 
	          char **seca, char **secb, char *fit, char *value, 
                  int cols, int printsec, int printdat, FILE *OUTPUT)
{
	int i;
	int count;
	int max;

	count=0;
	max=strlen(seqa[0]);
	fprintf(OUTPUT,"\nThe alignment:\n");
	while(count<max) {
	   fprintf(OUTPUT,"Position ");
	   for(i=count; i<count+cols-2; ++i) {
	      if((i+3)%10==0) {
		 fprintf(OUTPUT,"%3d",i+3);
		 i+=2;
	      } else fprintf(OUTPUT," ");
	      if(i>=max-2) break;
	   }
	   fprintf(OUTPUT,"\n");
	   for(i=0; i<na; ++i) {
	     fprintf(OUTPUT,"A(aa%3d): ",i+1);
	     dispone(seqa[i],cols,count,max,OUTPUT);
	   }
	   if(printsec) {
	     for(i=0; i<na; ++i) {
		fprintf(OUTPUT,"A(ss%3d): ",i+1);
		dispone(seca[i],cols,count,max,OUTPUT);
	     }
	   }
	   if(printdat) {
	     fprintf(OUTPUT,"FIT:      ");
	     dispone(fit,cols,count,max,OUTPUT);
	     fprintf(OUTPUT,"VALUE:    ");
	     dispone(value,cols,count,max,OUTPUT);
	   }
	   for(i=0; i<nb; ++i) {
	     fprintf(OUTPUT,"B(aa%3d): ",i+1);
	     dispone(seqb[i],cols,count,max,OUTPUT);
	   }
	   if(printsec) {
	      for(i=0; i<nb; ++i) {
		fprintf(OUTPUT,"B(ss%3d): ",i+1);
		dispone(secb[i],cols,count,max,OUTPUT);
	     }
	   }
	   count+=cols;
	   fprintf(OUTPUT,"\n");
	}
	return 0;

}
int dispone(char *seq, int cols, int count, int max, FILE *OUTPUT)
{
	int j;
      	for(j=0; j<cols; ++j)
	   if((count+j)<max) fprintf(OUTPUT,"%c",seq[count+j]);
      	fprintf(OUTPUT,"\n");
	return 0;
}
