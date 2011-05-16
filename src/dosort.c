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
#include "stamp.h"

int comp();

/*************************************************************************
dosort:  Version 2:  create the sortarr array, then copy the values
of each result into the array, freeing the memory as we go.
do qsort on the sortarr array, and return a pointer to the head of the 
array.
This  is a revised version that creates the sortarr array as we destroy the
result array.
-------------------------------------------------------------------------*/
/* SMJS Modified to remove structures  = */
struct path *dosort(struct olist *result, int *lena, int *total) {

    int i,j;
    int k=0;
    struct path *sortarr;

    sortarr = (struct path *) malloc(sizeof(struct path));

    for(i=0; i < ((*lena)-1); ++i){
	if(result[i].len > 0){	    /* if there are paths in this row */
	    for(j=0; j < result[i].len; ++j){
		sortarr = (struct path *) 
		    realloc(sortarr,sizeof(struct path) *(k+1));
/* SMJS		sortarr[k++] = result[i].res[j];  */
                CopyPath(&(sortarr[k++]),&(result[i].res[j]));
	    }
	}
	free(result[i].res); 
    }

    if(k != *total) {
       printf("k (%d) != total (%d) in dosort\n",k,*total); 
    }

    free(result);

    if(*total > 0){
	qsort((char *) sortarr, k, sizeof(struct path), comp);
    }
    return sortarr;

}

/*************************************************************************
comp:  compare two scores in the sortarr array
-------------------------------------------------------------------------*/
int comp(left,right)

struct path *left, *right;

{
    return right->score - left->score;
}
