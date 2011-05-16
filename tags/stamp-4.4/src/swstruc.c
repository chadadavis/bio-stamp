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

#define MAX3(A,B,C) ((A)>(B)?(A)>(C)?(A):(C):(B)>(C)?(B):(C))

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
/* SMJS Modified to remove structures  = */
int swstruc(const int  lena, const int lenb, const int pen, int ** prob, struct olist * const result,
        int * total, unsigned char ** patha, const int min_score) {

    register int diag, vert, horiz, rtemp,i , j, im1, done,k;

    unsigned char DIAG = 01; /* mask for diagonal move */
    unsigned char HORIZ= 02; /*          horizontal    */
    unsigned char VERT = 04; /*          vertical      */

    struct path *new;
    struct path *old;

    struct path *tempp;

    register int match = -1;

    register int minscore = min_score;

/* SMJS Added for speedup */
    int starti,endi,startj,endj;
    struct path *pathP;
    unsigned char flags;
    unsigned char *pathaj;

    *total = 0;

/* SMJS Changed malloc to calloc for zeroing */
    old = (struct path *) calloc(lena,sizeof(struct path));
    if(old == NULL) fprintf(stderr,"Cannot get space for old\n");

    new = (struct path *) calloc(lena,sizeof(struct path));
    if(new == NULL) fprintf(stderr,"Cannot get space for new\n");
    
    for (i=0; i< (lena-1); ++i){
/* SMJS already zeroed
	old[i].col = 0;
	old[i].score = 0;
*/
	old[i].start.i = i;
	old[i].start.j = 1;
	old[i].end.i = i;
	old[i].end.j = 1;
/* SMJS already zeroed
	new[i].col = 0;
	new[i].score =0;
*/
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
        pathaj = patha[j];
	for( i = 1,im1 = 0; i < (lena-1); ++i, ++im1){
/*	    im1 = i - 1; */
	    diag = old[im1].col + prob[i][j];
	    horiz= old[i].col - pen;
	    vert = new[im1].col - pen;
/*
	    rtemp = max4(diag,horiz,vert,0);
*/
	    rtemp = MAX3(0,diag,horiz);
            rtemp = max(rtemp,vert);
/*
	    patha[j][i] = 00;
	    if(rtemp == diag)	patha[j][i] = DIAG;
	    if(rtemp == horiz)	patha[j][i] = patha[j][i] | HORIZ;
	    if(rtemp == vert)	patha[j][i] = patha[j][i] | VERT;
*/
            flags = 00;
	    if(rtemp == diag)	flags = DIAG;
	    if(rtemp == horiz)	flags |= HORIZ;
	    if(rtemp == vert)	flags |= VERT;
/* indexing from 1 therefore ++pathaj */
            *(++pathaj) = flags;
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
#ifdef ASSIGNSTRUCT
			new[i].start = old[im1].start; 
#else
			new[i].start.i = old[im1].start.i;
			new[i].start.j = old[im1].start.j;
#endif
			if(rtemp >= old[im1].score){
			    new[i].score = rtemp;
			    new[i].end.i = i;
			    new[i].end.j = j;
			}
			else{
			    new[i].score = old[im1].score;
#ifdef ASSIGNSTRUCT
			    new[i].end = old[im1].end; 
#else
                            new[i].end.i = old[im1].end.i;
                            new[i].end.j = old[im1].end.j;
#endif
			}
		    }
		}
		else if(horiz == rtemp){
#ifdef ASSIGNSTRUCT
		    new[i].start = old[i].start;
#else
		    new[i].start.i = old[i].start.i;
		    new[i].start.j = old[i].start.j;
#endif
		    if(horiz >= old[i].score){
			new[i].score = horiz;
			new[i].end.i = i;
			new[i].end.j = j;
		    }else{
			new[i].score = old[i].score;
#ifdef ASSIGNSTRUCT
			new[i].end = old[i].end;
#else
                        new[i].end.i = old[i].end.i;
                        new[i].end.j = old[i].end.j;
#endif
		    }
		}
		else if(vert == rtemp){
#ifdef ASSIGNSTRUCT
		    new[i].start = new[im1].start; 
#else
		    new[i].start.i = new[im1].start.i;
		    new[i].start.j = new[im1].start.j;
#endif
		    if(vert > new[im1].score){
			new[i].score = vert;
			new[i].end.i = i;
			new[i].end.j = j;
		    }else{
			new[i].score = new[im1].score;
#ifdef ASSIGNSTRUCT
			new[i].end = new[im1].end;
#else
			new[i].end.i = new[im1].end.i;
			new[i].end.j = new[im1].end.j;
#endif
		    }
		}
		}
		if((i == (lena-2)) || (j == (lenb-2))){

                    pathP = &(new[i]);
                    starti = pathP->start.i;
                    startj = pathP->start.j;
                    endi = pathP->end.i;
                    endj = pathP->end.j;

/*
		    if((new[i].score >= minscore) &&
			(new[i].start.i > 0) &&
			(new[i].start.i != new[i].end.i) &&
		       (new[i].start.j != new[i].end.j)){
*/
		    if((new[i].score >= minscore) &&
			(starti > 0) &&
			(starti != endi) &&
		       (startj != endj)){
/* SMJS Combined condition */
/*
			done = present(&new[i],&result[new[i].start.i]);
			if(!done){
*/
			if(!present(pathP,&result[starti])){
			    addsco(pathP,&result[starti],total);
			}
			}
		}
	    else if(rtemp == 0){
		if(old[im1].score > 0){

                    pathP = &(old[im1]);
                    starti = pathP->start.i;
                    startj = pathP->start.j;
                    endi = pathP->end.i;
                    endj = pathP->end.j;

/*
		    if((old[im1].score >= minscore) &&
			(old[im1].start.i > 0) &&
			(old[im1].start.i != old[im1].end.i) &&
		       (old[im1].start.j != old[im1].end.j)){
*/
		    if((old[im1].score >= minscore) &&
			(starti > 0) &&
			(starti != endi) &&
		       (startj != endj)){
/* SMJS Combined condition */
/*
			   done = present(&old[im1],&result[old[im1].start.i]);
			   if(!done){
*/
			   if(!present(pathP,&result[starti])){
				addsco(pathP,&result[starti],total);
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
        {
#ifdef ASSIGNSTRUCT
	   new[k].start=new[k].end;
#else
/* SMJS Changed to zero for speed. Shouldn't affect results */
/*      because this loop is basically setting the path to  */
/*      zero length */
/*
	   new[k].start.i=new[k].end.i;
	   new[k].start.j=new[k].end.j;
*/
	   new[k].start.i=0;
	   new[k].start.j=0;
#endif
        }

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
int present(struct path *new, struct olist *result)
{
    int i;
    for(i = 0; i < result->len; ++i){
	if(new->start.j == result->res[i].start.j){
	    if(new->score > result->res[i].score){
#ifdef ASSIGNSTRUCT
		result->res[i] = *new; 
#else
                CopyPath(&(result->res[i]),new);
#endif
	    }
	    return 1;
	}
    }
    return 0;
}

/* SMJS Added routine CopyPath */
void CopyPath(struct path *to,struct path *from)
{
   to->start.i = from->start.i;
   to->start.j = from->start.j;
   to->end.i = from->end.i;
   to->end.j = from->end.j;
   to->score = from->score;
   to->col   = from->col;
}

/**************************************************************************
addsco
--------------------------------------------------------------------------*/
void addsco(struct path *new,struct olist *result,int *total)
{
    ++result->len;
    ++*total;
    if(result->len > 1){
	result->res = 
	(struct path *) realloc(result->res,sizeof(struct path)*result->len);
    }
#ifdef ASSIGNSTRUCT
    result->res[result->len - 1] = *new; 
#else
    CopyPath(&(result->res[result->len - 1]),new);
#endif
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
