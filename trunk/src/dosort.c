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
