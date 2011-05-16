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
#include <math.h>

#include "stamp.h"

#define MAXLINE 100

/* Reads in a tree-orderfile and a treefile and returns a
 *  structure containing the cluster information in a 
 *  (relatively) easy to read form.  Returns NULL if an
 *  error occurs, the said structure otherwise.  
 * One important thing to know is that although the program
 *  reads in clusters/orders for elements numbered 
 *  (say) 1..N, it returns clusters numbered from
 *        0..(N-1).
 *
 * For example: 
 *  given the tree-order file:
 *          3                0.00         0
 *	    2                0.00         0
 *	    1                0.00         0
 *  and the tree file:
 * 	   1
 *     	   1
 *	   1
 * 	   2
 *         2
 *         1    2
 *	   1
 *	   3
 *  it will return the following clusters:
 *    	cl[0].a.number=1
 *	cl[0].a.member[0]=2  (ie. 3-1)
 *	cl[0].b.number=1
 *	cl[0].b.member[0]=1  (ie. 2-1)
 *	cl[1].a.number=2
 *	cl[1].a.member[0]=2, cl[1].a.member[1]=1
 *	cl[1].b.number=1
 *      cl[1].b.member[0]=0  (ie. 1-1)
 *  implying the tree:
 *      2----+
 *           |-------+
 *      1----+       |_____
 *                   |
 *      0------------+
 *      
 *	*number is the number of elements considered.
 *  Therefor a loop such as for(i=0; i<(*number-1); ++i) 
 *   will allow analysis of each cluster one at a time.
 *   (ie. at cluster i, the members of
 *      cl[i].a and cl[i].b are being brought together) 
 *
 * Change: November 20, 1991:
 *   method is a flag to specify what information is returned
 *	0  return the tree information considering the order file
 *	1  return the tree information ignoring the order file */


struct cluster *readtree(char *tordfile, char *treefile, int *number,
	int method, FILE *OUT) {

	int i,j,k,*ord;
	FILE *f;
	char *buff,*addbuff;
	struct cluster *cl;

	ord=(int*)malloc(sizeof(int));

	/* First the order must be extracted from the order file */
	if(method==0) {
	  if((f=fopen(tordfile,"r")) == NULL) {
	   fprintf(OUT,"readtree: cannot open file %s\n",tordfile);
	   return NULL;
	  } else {
	   *number=0;
	   buff=(char*)malloc((unsigned)MAXLINE*sizeof(char));
	   addbuff=buff;
	   while((buff=fgets(buff,100,f)) !=NULL) {
	     sscanf(buff,"%d ",&ord[(*number)]);
	     (*number)++;
	     ord=(int*)realloc(ord,(*number+1)*sizeof(int));
	   }
	   free(addbuff);
	  } /* end of if((f... */
	  fclose(f); 
	}

	/* Now the tree file may be opened and the tree information
	 *  read in */
	cl=(struct cluster*)malloc((*number)*sizeof(struct cluster));
	if((f=fopen(treefile,"r")) == NULL) {
	   fprintf(OUT,"readtree: cannot open file %s\n",treefile);
	   return NULL;
        } else {
	   for(i=0; i<(*number-1); ++i) {
	     /* NB: the number of clusters is ALWAYS one 
	      *   less than the number of elements */
	     fscanf(f,"%d",&cl[i].a.number);
	     cl[i].a.member=(int*)malloc(cl[i].a.number*sizeof(int));
	     for(j=0; j<cl[i].a.number; ++j) {
	       fscanf(f,"%d",&k); 
	       if(method==1) cl[i].a.member[j]=k-1;
	       else cl[i].a.member[j]=ord[k-1]-1;
	       }
	     fscanf(f,"%d",&cl[i].b.number);
	     cl[i].b.member=(int*)malloc(cl[i].b.number*sizeof(int)); 
	     for(j=0; j<cl[i].b.number; ++j) { 
	       fscanf(f,"%d",&k);  
	       if(method==1) cl[i].b.member[j]=k-1;
	       else cl[i].b.member[j]=ord[k-1]-1; 
	       }
	    } /* End of for */
	} /* End of if((f... */
	fclose(f);
	if(method==0) free(ord);
	return cl;
}
