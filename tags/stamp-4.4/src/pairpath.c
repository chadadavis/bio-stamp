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
#include <math.h>

#include "stamp.h"


/* A path through the matrix if found using the Smith and Waterman
 * algorithm coded in Geoff Barton's routines swstruc, dosort and 
 * aliseq.  */ 

float pairpath(struct domain_loc domain1, struct domain_loc domain2, int **prob,
	long int entry, float **R2, float *V2, int *len, float *score,
	int *nfit, char *ftouse, float *fpuse, int *start1, int *end1,
	int *start2,int *end2, struct parameters *parms) {

	int i,j,k,q, match;
	int allen,total,ia,ib;
	int finala,finalb,diffab;
	int pen;
	int p1,p2;

	unsigned char **patha;
	struct olist *result;
	struct path *sortarr;

	long int  highest,second,stl;
	int best,secbest,acount,bcount;

	float ratio,letter;
	int minscore;
	int **natoms1, **natoms2;

	char *bestaseq, *bestbseq;
	char *touse;
	float *puse;
	float conf;

	float RMSE;

/* SMJS Added tmplen */
        int tmplen;
/* SMJS Added for speedup */
        int lenstr;


	/* Allocating memory */

	bestaseq = (char*)malloc((parms[0].MAX_SEQ_LEN)*sizeof(char));
	bestbseq = (char*)malloc((parms[0].MAX_SEQ_LEN)*sizeof(char));
	puse =  (float*)malloc((parms[0].MAX_SEQ_LEN)*sizeof(float));
	touse =  (char*)malloc((parms[0].MAX_SEQ_LEN)*sizeof(char));

	patha = (unsigned char**)malloc((unsigned)(domain2.ncoords+2)*sizeof(unsigned char*));
	for(i=0;i<(domain2.ncoords+2); ++i)
	    patha[i]=(unsigned char *)malloc((unsigned)(domain1.ncoords+2)*sizeof(unsigned char));
	result = (struct olist*)malloc((unsigned)sizeof(struct olist) *(domain1.ncoords+2));
	
	(*nfit)=(*len)=0;
	RMSE=100.0;
	(*score)=0.0;

	/* Calling the Smith and Waterman routine */

	minscore=0;
	highest=secbest=second=best=0;
	if(parms[0].BOOLEAN) {
	   minscore=parms[0].MINFIT;
	   pen=0;
	} else {
	   minscore=parms[0].MINFIT*(int)((float)parms[0].PRECISION*parms[0].CUTOFF);
	   pen=(int)(parms[0].PAIRPEN*(float)parms[0].PRECISION);
	}

/*
int sw7ccs(int  lena, int lenb, int **prob, int pen, struct olist *result,
        int *total, unsigned char **patha, int min_score, 
        int auto_corner, float fk);

int swstruc(int  lena, int lenb, int pen, int **prob, struct olist *result,
        int *total, unsigned char **patha, int min_score);
*/

	if(parms[0].SW==1) { 
	  match = sw7ccs (domain1.ncoords+2,domain2.ncoords+2,prob,pen,result,&total,patha,minscore,
			parms[0].CCADD,parms[0].CCFACTOR);
	} else  { 
	  match = swstruc(domain1.ncoords+2,domain2.ncoords+2,pen,prob,result,&total,patha,minscore);
        }
    (*start2)=(*end2)=0;
    if(total!=0) {
/* SMJS Added tmplen stuff to stop k!=total errors in dosort. I'm not sure if */
/*      this is correct. */
/* SMJS sortarr = dosort(result,&domain1.ncoords,&total); */
        tmplen = domain1.ncoords+1;
        sortarr = dosort(result,&tmplen,&total);
/* SMJS End */

	i=0; /* i=0 will always be the alignment we want */

/* SMJS Changed sortarr[i] to sortarr[0] for clarity */
	aliseq(domain1.aa,domain2.aa,&sortarr[0],patha,bestaseq,bestbseq,&allen,&ia,&ib);
/* SMJS Changed sortarr[i] to sortarr[0] for clarity */
	highest=stl=sortarr[0].score; best=0;
	if( (sortarr[0].start.i+allen+(strlen(domain1.aa)-sortarr[0].end.i))>parms[0].MAX_SEQ_LEN ||
                (sortarr[0].start.j+allen+(strlen(domain2.aa)-sortarr[0].end.j))>parms[0].MAX_SEQ_LEN) {
                  fprintf(stderr,"error: MAX_SEQ_LEN limit surpassed, set value > %d in parameter file or command line\n",
                    parms[0].MAX_SEQ_LEN);
                  return -1;
        }


	  

	(*len)=allen;
/* SMJS Changed sortarr[i] to sortarr[0] for clarity */
	(*start1)=sortarr[0].start.i-1;
	(*end1)=sortarr[0].end.i-1;
	(*start2)=sortarr[0].start.j-1;
	(*end2)=sortarr[0].end.j-1;

	/* Now it is necessary to make use of the results, initially it is
	 *   probably best to use the best path to arrive at a set of
	 *   coordinates for use in least squares fitting. 
	 *
	 * Locating the coordinates to use. */
#ifdef DBGPREC
        fprintf(parms[0].LOG,"start1 = %d end1 = %d start2 = %d end2 = %d len = %d highest = %d\n",
               *start1,*end1,*start2,*end2,*len,highest);
#endif

/* SMJS Changed sortarr[best] to sortarr[0] for clarity */
	acount=sortarr[0].start.i;
	bcount=sortarr[0].start.j;
	j=0;
	/* memory for natoms1 and natoms2 is allocated */
	natoms1=(int**)malloc((allen+10)*sizeof(int*));
	natoms2=(int**)malloc((allen+10)*sizeof(int*)); 

	/* i runs from start to finish of the alignment, j represents a
	 *  count of the number of atoms to be fitted, and acount & bcount
	 *  are pointers (not 'C' pointers, just pointers!) to which residue
	 *  of each sequence/structure we are presently dealing with.
	 *  In this version of pairpath, equlivalences are chosen by
	 *  1) whether or not they are aligned by the SW algorithm.
	 *  2) whether the value of prob is above a predetermined
	 *     parms[0].CUTOFF value.  */
	q=0;

#ifdef DBGSTEVE
        printf("bestaseq = %s\n",bestaseq);
        printf("bestbseq = %s\n",bestbseq);
#endif
	for(i=0; i<allen; ++i) {
	  puse[i]=0.0;
	/* Displaying the alignment:  (only if parms[0].PAIRALIGN || parms[0].PAIRALLALIGN is set to 1)
	 * touse shows which residues from the alignment are to be used
	 *   for fitting (ie. probabilities are greater than parms[0].CUTOFF. 
	 *   puse is a string of floats representing the
	 *   probabilities. */
  	if(parms[0].PAIROUTPUT || parms[0].PAIRALIGN || parms[0].PAIRALLALIGN) { touse[i]=' '; touse[i+1]='\0'; puse[i]=0.0;  }

	   if(bestaseq[i/* SMJS+1*/] !=' ' && bestbseq[i/* SMJS+1*/] !=' ') {
	      ++q;

/*	    if(parms[0].PAIROUTPUT || parms[0].PAIRALIGN || parms[0].PAIRALLALIGN) { */
	     if(!parms[0].BOOLEAN) {
	      puse[i]=(float)prob[acount][bcount]/parms[0].PRECISION;
	     } else {
	       if(prob[acount][bcount])
	          puse[i]=1.0;
	       else 
		  puse[i]=0.0;
	     }
/*	    } */
	      if( (prob[acount][bcount] >= (int)(parms[0].CUTOFF*(float)parms[0].PRECISION) && !parms[0].BOOLEAN) ||
		  (prob[acount][bcount] == 1 && parms[0].BOOLEAN) ) {
#ifdef TESTWRONG
                 if (acount<domain1.ncoords && bcount<domain2.ncoords)
                 {
	            natoms1[j]=domain1.coords[acount];
                    printf("Coords for %d %c = %d %d %d\n",acount,domain1.aa[acount],natoms1[j][0],natoms1[j][1],natoms1[j][2]);
	            natoms2[j]=domain2.coords[bcount];
	            j++;
	            touse[i]='1'; touse[i+1]='\0';
                 }
                 else
                 {
	            touse[i]='0'; touse[i+1]='\0';
                 }
#else
	         natoms1[j]=domain1.coords[acount-1];
	         natoms2[j]=domain2.coords[bcount-1];
	         j++;
	         touse[i]='1'; touse[i+1]='\0';
#endif
	      } else {
                 touse[i]='0';  /* End of if(prob... */
		 touse[i+1]='\0';
	      }
	   } else { touse[i]=' '; touse[i+1]='\0'; puse[i]=-1.0; } /* End of if(best... */
	acount += (bestaseq[i/* SMJS+1*/] != ' ');
	bcount += (bestbseq[i/* SMJS+1*/] != ' ');
	}  /* End of for loop  */
/*	if(parms[0].PAIROUTPUT || parms[0].PAIRALIGN || parms[0].PAIRALLALIGN) */
	touse[i]='\0';
	
	(*nfit)=j;

/*	if(parms[0].PAIRALIGN || parms[0].PAIROUTPUT) { */
	 /* create a copy of the current alignment in domain1.align and
	  * domain2.align -- see treepath.c for a step by step explanation */
	 finala=finalb=0;
	 for(j=0; j<(sortarr[0].start.i-sortarr[0].start.j); ++j) {
	    domain2.align[finalb]=ftouse[finalb]=' ';
	    fpuse[finalb]=-1.0;
	    finalb++;
	 }
         for(j=0; j<(sortarr[0].start.j-sortarr[0].start.i); ++j) {
	    domain1.align[finala]=ftouse[finala]=' ';
	    fpuse[finala]=-1.0;
	    finala++;
	 }
/* SMJS Added -1 because the start.i is 1 based*/
	 for(j=0; j<sortarr[0].start.i-1; ++j) {
	    domain1.align[finala]=domain1.aa[j];
	    fpuse[finala]=-1.0;
	    ftouse[finala]=' ';
	    finala++;
	 }
/* SMJS Added -1 because the start.j is 1 based*/
	 for(j=0; j<sortarr[0].start.j-1; ++j) {
	    domain2.align[finalb]=domain2.aa[j];
	    fpuse[finalb]=-1.0;
	    ftouse[finalb]=' ';
	    finalb++;
	 }
/* SMJS Changed 1 to 0 because the sequence starts at 0 */
	 sprintf(&domain1.align[finala],"%s",&bestaseq[/* SMJS1*/0]); 
	 sprintf(&domain2.align[finalb],"%s",&bestbseq[/* SMJS1*/0]);
/* SMJS Added lenstr */
         lenstr = strlen(touse);
	 for(i=0; i<lenstr; ++i) fpuse[finala+i]=puse[i];
	 sprintf(&ftouse[finalb],"%s",touse);
	 finala+=strlen(bestaseq);
	 finalb+=strlen(bestbseq);
/* SMJS Added lenstr */
         lenstr = strlen(domain1.aa);
         for(j=sortarr[0].end.i; j<lenstr; ++j) {
	   domain1.align[finala]=domain1.aa[j];
	   fpuse[finala]=-1.0;
	   ftouse[finala]=' ';
	   finala++;
	 }
/* SMJS Added lenstr */
         lenstr = strlen(domain2.aa);
	 for(j=sortarr[0].end.j; j<lenstr; ++j) {
	    domain2.align[finalb]=domain2.aa[j];
	    fpuse[finalb]=-1.0;
	    ftouse[finalb]=' ';
	    finalb++;
	 }
	 diffab=(finala-finalb);
	 for(j=0; j<(-1*diffab); ++j) {
	    domain1.align[finala]=ftouse[finalb]=' ';
	    fpuse[finala]=-1.0;
	    finala++;
	 }
	 for(j=0; j<diffab; ++j)  {
	    domain2.align[finalb]=ftouse[finalb]=' '; 
	    fpuse[finalb]=-1.0;
	    finalb++;
	 }
	 domain1.align[finala]=domain2.align[finalb]=ftouse[finala]='\0';
	 domain1.align[finalb]=domain2.align[finala]=ftouse[finalb]='\0';
#ifdef DBGSTEVE
         printf("domain1.align = %s\n",domain1.align);
         printf("domain2.align = %s\n",domain2.align);
         printf("ftouse        = %s\n",ftouse);
#endif
	 fpuse[finala]=fpuse[finalb]=-1.0;
	 if(finala!=finalb) fprintf(parms[0].LOG,"Something funny in pairpath.c\n");
/*	} */


	if((*nfit)>=3) {
	   RMSE=matfit(natoms1,natoms2,R2,V2,(*nfit),entry,parms[0].PRECISION);
	} else { 
	   fprintf(parms[0].LOG,"WARNING: matfit NOT called as there were less than three equivalences.\n");  
	   RMSE=100.0; 
	}
	if(parms[0].PAIRWISE) {
	   if((allen>0) && (domain1.ncoords>0) && (domain2.ncoords>0)) {
	     (*score)=((float)highest/(float)allen)*((float)(allen-ia)/(float)(domain1.ncoords))*((float)(allen-ib)/(float)(domain2.ncoords));
	   } else {
	     (*score) = 0.0;
	   }
	} else { /* if scanning, then we don't want to penalise for length of the database sequence */
	   if((allen>0) && (domain1.ncoords>0) && (domain2.ncoords>0)) {
	     switch(parms[0].SCANSCORE) {
	     case 1: (*score)=((float)highest/(float)allen)*((float)(allen-ia)/(float)(domain1.ncoords))*((float)(allen-ib)/(float)(domain2.ncoords)); break;
	     case 2: (*score)=((float)highest/(float)allen)*((float)(allen-ia)/(float)(allen))*((float)(allen-ib)/(float)(allen)); break;
	     case 3: (*score)=((float)highest/(float)allen)*((float)(allen-ia)/(float)(domain1.ncoords)); break;
	     case 4: (*score)=((float)highest/(float)allen)*((float)(allen-ib)/(float)(domain2.ncoords)); break;
	     case 5: (*score)=((float)highest/(float)allen); break;
	     case 6: (*score)=((float)highest/(float)domain1.ncoords)*((float)(allen-ia)/(float)(domain1.ncoords)); break;
	     default: (*score)=((float)highest/(float)allen);	
	     }
	   } else {
		(*score) = 0.0;
	   }
		 
	} 
	if(!parms[0].BOOLEAN) (*score)/=(float)parms[0].PRECISION;
	/* The above makes the score dependant upon:
	 *  1) the sum of reliability values for the found alignment,
	 *  2) the length of the alignment, and
	 *  3) the ratio of alignment length to sequence length for
  	 *     each sequence compared. 
	 */

	free(natoms1); 
	free(natoms2);

#ifdef DBGSTEVE
        if (*score>100)
        {
           for(i=0;i<domain1.ncoords+2;i++)
           {
              printf("%4d",i); 
              for (j=0;j<domain2.ncoords+2;j++)
              {
                 printf("|%4d|",prob[i][j]);
              }
              printf("\n");
           }
        }
#endif
	
	/* Modification: now we redifine the start, lengths and ends to ignore
 	 *  un-equivalent regions lying off the ends */

/* SMJS I don't see why p1 and p2 were set to sortarr[0].start.i-1 and */
/*      sortarr[0].start.j-1 respectively. The equivalent prob array element */
/*      to bestaseq[0] and bestbseq[0] is  */
/*      prob[sortarr[0].start.i][sortarr[0].start.j] */

	p1=sortarr[0].start.i/* SMJS-1*/; p2=sortarr[0].start.j/* SMJS-1*/;
	for(j=0; j<allen-1;  ++j) {
	  if(bestaseq[j]!=' ' && bestbseq[j]!=' ' && prob[p1][p2]>=parms[0].CUTOFF &&
	     bestaseq[j+1]!=' ' && bestbseq[j+1]!=' ' && 
	     p1<domain1.ncoords-1 && p2<domain2.ncoords-1 && prob[p1+1][p2+1]>=parms[0].CUTOFF) {
	     (*start1)=p1;
	     (*start2)=p2;
	     break;
	  }
	  if(bestaseq[j]!=' ') p1++;
	  if(bestbseq[j]!=' ') p2++;
#ifdef DBGSTEVE
          printf("Did a start iter\n");
#endif
	}
/* SMJS I don't see why p1 and p2 were set to sortarr[0].end.i-1 and */
/*      sortarr[0].end.j-1 respectively. The equivalent prob array element */
/*      to bestaseq[allen-1] and bestbseq[allen-1] is  */
/*      prob[sortarr[0].end.i][sortarr[0].end.j] */

	p1=sortarr[0].end.i/* SMJS-1*/; p2=sortarr[0].end.j/* SMJS-1*/;
	for(j=allen-1; j>1; --j) {
	  if(bestaseq[j]!=' ' && bestbseq[j]!=' ' && prob[p1][p2]>parms[0].CUTOFF &&
             bestaseq[j-1]!=' ' && bestbseq[j-1]!=' ' && prob[p1-1][p2-1]>parms[0].CUTOFF) {
             (*end1)=p1;
             (*end2)=p2;
             break;
          }
	  if(bestaseq[j]!=' ') p1--;
          if(bestbseq[j]!=' ') p2--;
#ifdef DBGSTEVE
          printf("Did an end iter\n");
#endif
	} 
     } else { 
       fprintf(parms[0].LOG,"No paths found\n");
     } /* end of if(total==0... */

	/* Freeing memory allocated for patha and result */
	for(i=0;i<(domain2.ncoords+2); ++i)
		free((char*)patha[i]);
	free((char*)patha);

	/* Freeing other memory */
	if(total!=0) free(sortarr);
	free(bestaseq);
	free(bestbseq);
	free(puse);
	free(touse);
	return RMSE;


}  /* pairpath */
