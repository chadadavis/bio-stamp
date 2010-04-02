/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)	            Prof. Geoffrey J. Barton (GJB)
 EMBL Heidelberg                            School of Life Sciences
 Meyerhofstrasse 1                          University of Dundee
 D-69117 Heidelberg                         Dow Street
 Germany                                    Dundee, DD1 5EH
                                          
 Tel: +49 6221 387 473                      Tel: +44 1382 345860
 FAX: +44 6221 387 517                      FAX: +44 1382 345764
 E-mail: russell@embl-heidelberg.de         E-mail geoff@compbio.dundee.ac.uk
 WWW: http://www.russell.emb-heidelberg.de  WWW: http://www.compbio.dundee.ac.uk

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
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
