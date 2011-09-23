/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 Appendix A of the file LICENSE.  Conditions and restrictions for use of
 this package are also in this file.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Biomolecular Modelling Laboratory   Laboratory of Molecular Biophysics
 Imperial Cancer Research Fund       The Rex Richards Building
 Lincoln's Inn Fields, P.O. Box 123  South Parks Road
 London, WC2A 3PX, U.K.              Oxford, OX1 3PG, U.K.
 Tel: +44 171 269 3583               Tel: +44 865 275368
 FAX: +44 171 269 3417               FAX: 44 865 510454
 e-mail: russell@icrf.icnet.uk       e-mail gjb@bioch.ox.ac.uk
 WWW: http://bonsai.lif.icnet.uk/    WWW: http://geoff.biop.ox.ac.uk/

 The WORK is Copyright (1992,1993,1995,1996) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <stamp.h>

/**********************************************************************
sw7: routine to do the Smith Waterman algorithm and retain history of
best paths 
Revised storage version uses a single structure for path and col.
Fastest version so far.

sw7: modification of sw2 to store the path results in an array 
rather than write each one out.
sw7: Also stores path array to allow tracing out of alignments.
see sw6 for further details.

----------------------------------------------------------------------*/
int swstruc(int  lena, int lenb, int pen, int **prob, struct olist *result,
	int *total, unsigned char **patha, int min_score) {

    register int diag, vert, horiz, rtemp,i , j, im1, done,k;

    unsigned char DIAG = 01; /* mask for diagonal move */
    unsigned char HORIZ= 02; /*          horizontal    */
    unsigned char VERT = 04; /*          vertical      */

    struct path *new;
    struct path *old;

    struct path *tempp;

    register int match = -1;

    register int minscore = min_score;

    *total = 0;

    old = (struct path *) malloc(sizeof(struct path)*lena);
    if(old == NULL) fprintf(stderr,"Cannot get space for old\n");

    new = (struct path *) malloc(sizeof(struct path)*lena);
    if(new == NULL) fprintf(stderr,"Cannot get space for new\n");

    for (i=0; i< (lena-1); ++i){
	old[i].col = 0;
	old[i].score = 0;
	old[i].start.i = i;
	old[i].start.j = 1;
	old[i].end.i = i;
	old[i].end.j = 1;
	new[i].col = 0;
	new[i].score =0;
	new[i].start.i = i;
	new[i].start.j = 1;
	new[i].end.i = i;
	new[i].end.j = 1;
	result[i].len = 0;	
	result[i].res = (struct path *) malloc(sizeof(struct path));
	patha[0][i] = 0;
    }
    for (j = 0; j  < (lenb-1); ++j)
	patha[j][0] = 0;

/*    printf("lena: %d, lendb: %d\n",lena,lenb); */
    for(j = 1; j < (lenb-1); ++j){
	for( i = 1; i < (lena-1); ++i){
	    im1 = i - 1;
	    diag = old[im1].col + prob[i][j];
	    horiz= old[i].col - pen;
	    vert = new[im1].col - pen;
	    rtemp = max4(diag,horiz,vert,0);
	    patha[j][i] = 00;
	    if(rtemp == diag)	patha[j][i] = DIAG;
	    if(rtemp == horiz)	patha[j][i] = patha[j][i] | HORIZ;
	    if(rtemp == vert)	patha[j][i] = patha[j][i] | VERT;
	    if(rtemp > 0){
		if(diag == rtemp){
		    if(old[im1].col == 0){
			new[i].start.i = i;
			new[i].start.j = j;
			new[i].score = rtemp;
			new[i].end.i = i;
			new[i].end.j = j;
		    }
		    else{
			new[i].start = old[im1].start;
			if(rtemp >= old[im1].score){
			    new[i].score = rtemp;
			    new[i].end.i = i;
			    new[i].end.j = j;
			}
			else{
			    new[i].score = old[im1].score;
			    new[i].end = old[im1].end;
			}
		    }
		}
		else if(horiz == rtemp){
		    new[i].start = old[i].start;
		    if(horiz >= old[i].score){
			new[i].score = horiz;
			new[i].end.i = i;
			new[i].end.j = j;
		    }else{
			new[i].score = old[i].score;
			new[i].end = old[i].end;
		    }
		}
		else if(vert == rtemp){
		    new[i].start = new[im1].start;
		    if(vert > new[im1].score){
			new[i].score = vert;
			new[i].end.i = i;
			new[i].end.j = j;
		    }else{
			new[i].score = new[im1].score;
			new[i].end = new[im1].end;
		    }
		}
		}
		if((i == (lena-2)) || (j == (lenb-2))){
		    if((new[i].score >= minscore) &&
			(new[i].start.i > 0) &&
			(new[i].start.i != new[i].end.i) &&
		       (new[i].start.j != new[i].end.j)){
			done = present(&new[i],&result[new[i].start.i]);
			if(!done){
			    addsco(&new[i],&result[new[i].start.i],total);
			}
			}
		}
	    else if(rtemp == 0){
		if(old[im1].score > 0){
		    if((old[im1].score >= minscore) &&
			(old[im1].start.i > 0) &&
			(old[im1].start.i != old[im1].end.i) &&
		       (old[im1].start.j != old[im1].end.j)){
			   done = present(&old[im1],&result[old[im1].start.i]);
			   if(!done){
				addsco(&old[im1],&result[old[im1].start.i],
				total);
			   }
					   
		       }
		       new[i].score = 0;
		}
	    }
	    if(rtemp > match) match = rtemp;
	    new[i].col = rtemp;
	}
        /* switch the array pointers  - old for new*/

	tempp = old;
	old = new;
	new = tempp;
	for(k=0; k<(lena-1); ++k) 
	   new[k].start=new[k].end;

    }
/*
    for(j = 1; j < (lenb-1); ++j){ 
	for( i = 1; i < (lena-1); ++i){
	   printf("%1d",(int)patha[j][i]);
	   }
	   printf("\n");
	   }
*/
    free(new);
    free(old);
    return match;  /* the biggest value found */
}

/****************************************************************************
present:  If the new path (new) has the same start point as a path in result
then check if the new score is > the old score.  If it is, then store the
new path in the result and return 1.  If it isn't then do nothing, but return
0.
----------------------------------------------------------------------------*/
present(new, result)

struct path *new;
struct olist *result;

{
    int i;
    for(i = 0; i < result->len; ++i){
	if(new->start.j == result->res[i].start.j){
	    if(new->score > result->res[i].score){
		result->res[i] = *new;
	    }
	    return 1;
	}
    }
    return 0;
}
/**************************************************************************
addsco
--------------------------------------------------------------------------*/
addsco(new,result,total)

struct path *new;
struct olist *result;
int *total;

{
    ++result->len;
    ++*total;
    if(result->len > 1){
	result->res = 
	(struct path *) realloc(result->res,sizeof(struct path)*result->len);
    }
    result->res[result->len - 1] = *new;
}
/***************************************************************************
ppath - print out the paths as stored in result array
---------------------------------------------------------------------------*/
ppath(result,lena)

struct olist *result;
int lena;

{
    int i,j;

    for(i=0; i < lena; ++i){
	if(result[i].len > 0){
	    for(j=0; j < result[i].len; ++j){
		printf("%d %d %d %d %d\n",result[i].res[j].start.i,
					  result[i].res[j].start.j,
					  result[i].res[j].end.i,
					  result[i].res[j].end.j,
					  result[i].res[j].score);
	    }
	}
    }
}

/***************************************************************************
ppath2 - print out the paths as pointed to by sortarr
---------------------------------------------------------------------------*/
ppath2(sortarr,total)

struct path **sortarr;
int total;

{


    int i;
    for(i=0; i < total; ++i){
	printf("%d %d %d %d %d %d\n",i,sortarr[i]->start.i,
				  sortarr[i]->start.j,
				  sortarr[i]->end.i,
				  sortarr[i]->end.j,
				  sortarr[i]->score);
    }

}
/*************************************************************************
ppath3 write out genplot file for the paths obtained
-------------------------------------------------------------------------*/
ppath3(sortarr,total,lena,lenb)

struct path **sortarr;
int total,lena;
{
    int i;
    FILE *fp, *fopen();
    fp = fopen("test.genplot","w");

    fprintf(fp,"!PAPER = A4H\n");
    fprintf(fp,"!XMIN = 1\n");
    fprintf(fp,"!XMAX = %d\n",lenb+1);
    fprintf(fp,"!YMIN = 1\n");
    fprintf(fp,"!YMAX = %d\n",lena+1);
    fprintf(fp,"!START_RANGE\n");

    for(i=0; i < total; ++i){
	fprintf(fp,"%d, %d, %d, %d \n",sortarr[i]->start.j,
				  sortarr[i]->start.i,
				  sortarr[i]->end.j,
				  sortarr[i]->end.i);
    }

    fprintf(fp,"!END_RANGE\n");
    fclose(fp);
}

