#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <stamp.h>

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
