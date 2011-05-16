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
#include <stdlib.h>

#include "stamp.h"

int treefit(struct domain_loc *domain, int ndomain, struct cluster cl, 
	float *score, float *rms, int *length, int *nfit,
	float *Pij, float *Dij, float *dist, float *Pijp,
	int rev, int align, struct parameters *parms) {

	char **pseq1,**pseq2,**psec1,**psec2;
	char *puse,*touse;

	int i,j,k,l,m,n;
	int iter;
	int pysize, pxsize;
	int indx,indy;
	int scorerise;
	int no_comparisons,no_matrix,ave_length;
	int *xcount,*ycount;
	int **prob;

	float oldscore;
	float sum,sumsq;
	float mean,sd;
	float pairwise_mean,pairwise_sd;
	float diff;
	float D,P,C;
	float Rossmann;
	float **R2,*V2,**dR2,*dV2;
/* SMJS Added variables for inverse consts in rossmann */
        float const1,const2,prec2i;

	FILE *IN;

	R2=(float**)malloc(3*sizeof(float*));
	for(i=0; i<3; ++i) 
	  R2[i]=(float*)malloc(3*sizeof(float));
	V2=(float*)malloc(3*sizeof(float));
	puse=(char*)malloc(parms[0].MAX_SEQ_LEN*sizeof(char));
	touse=(char*)malloc(parms[0].MAX_SEQ_LEN*sizeof(char));
	xcount=(int*)malloc(ndomain*sizeof(int));
	ycount=(int*)malloc(ndomain*sizeof(int));
	pseq1=(char**)malloc(cl.a.number*sizeof(char*));
	psec1=(char**)malloc(cl.a.number*sizeof(char*));
	pseq2=(char**)malloc(cl.b.number*sizeof(char*));
	psec2=(char**)malloc(cl.b.number*sizeof(char*));

	pxsize=strlen(domain[cl.a.member[0]].align)+1;
	pysize=strlen(domain[cl.b.member[0]].align)+1;	

	prob=(int**)malloc((pxsize+2) *sizeof(int*));
	for(j=0; j<(pxsize+2); ++j)
/* SMJS Was sizeof(int *) */
	   prob[j]=(int*)malloc((pysize+2)*sizeof(int));

        iter=0; oldscore=0.0; diff=parms[0].SCORETOL+1;
	(*nfit)=100; scorerise=1;
	while(iter<parms[0].MAXTITER && diff>parms[0].SCORETOL && (*nfit)>3 && (scorerise || parms[0].SCORERISE!=1)) {
	   for(j=0; j<ndomain; ++j) xcount[j]=ycount[j]=0;
	   /* first apply the current transformation (if there is one) 
	    *  an copy the oldalignment into the current alignment */
	   if(iter>0) {
	      /* copy old alignment back into align to do another iteration */
	      for(j=0; j<cl.a.number; ++j) 
	          strcpy(domain[cl.a.member[j]].align,domain[cl.a.member[j]].oldalign);
	      for(j=0; j<cl.b.number; ++j) 
		  strcpy(domain[cl.b.member[j]].align,domain[cl.b.member[j]].oldalign);
	   }

	   sum=sumsq=0.0;
	   no_matrix=(pxsize-1)*(pysize-1);
	   no_comparisons=cl.a.number*cl.b.number;

/* SMJS Added const1 and const2 and prec2i */
          prec2i=1.0/(float)(parms[0].PRECISION*parms[0].PRECISION);
          const1=(1.0/parms[0].const1)*prec2i;
          const2=(1.0/parms[0].const2)*prec2i;

	  for(j=0; j<cl.a.number; ++j) xcount[j]=0;
          for(l=0; l<strlen(domain[cl.a.member[0]].align); ++l) {
	    for(k=0; k<cl.b.number; ++k) ycount[k]=0;
	    for(m=0; m<strlen(domain[cl.b.member[0]].align); ++m) {
	       P=0.0;
	       for(j=0; j<cl.a.number; ++j) {
		  indx=cl.a.member[j];
		  for(k=0; k<cl.b.number; ++k) {
		     indy=cl.b.member[k];
	             if(j==0 && k==0) prob[l+1][m+1]=(parms[0].BOOLEAN);
	             if(domain[indx].align[l]!=' ' && domain[indy].align[m]!=' ') {
/* SMJS Changed to use inverse constants (its faster)
	 	          Rossmann=rossmann(&domain[indx].coords[xcount[j]],&domain[indy].coords[ycount[k]],
		  	       (xcount[j]==0 || ycount[k]==0),
			       (xcount[j]==(domain[indx].ncoords-1) || ycount[k]==(domain[indy].ncoords-1)),
				parms[0].const1,parms[0].const2,&D,&C,parms[0].PRECISION);
*/
	 	          Rossmann=rossmann(&domain[indx].coords[xcount[j]],&domain[indy].coords[ycount[k]],
		  	       (xcount[j]==0 || ycount[k]==0),
			       (xcount[j]==(domain[indx].ncoords-1) || ycount[k]==(domain[indy].ncoords-1)),
				const1,const2,&D,&C);
			  if(!parms[0].BOOLEAN) P+=Rossmann;
			  else prob[l+1][m+1]*=(Rossmann>=parms[0].BOOLCUT);
			  if(P>(1.0*(float)no_comparisons)) {
			     fprintf(parms[0].LOG,"something is wrong with the coordinates\n");
			     fprintf(parms[0].LOG," at position A(%d) %d, B(%d) %d\n",j,l,k,m);
			   }
		    } else { 
		        P+=parms[0].TREEPEN;
		    }
		 }
	       }
	       if(!parms[0].BOOLEAN) {
	          sum+=P;
	          sumsq+=P*P;
	          prob[l+1][m+1]+=(int)(((float)parms[0].PRECISION)*P/(float)no_comparisons);
	       }
	       for(k=0; k<cl.b.number; ++k) ycount[k]+=(domain[cl.b.member[k]].align[m]!=' ');
	     }
	     for(j=0; j<cl.a.number; ++j) xcount[j]+=(domain[cl.a.member[j]].align[l]!=' ');
	   }

	   /* We must correct the sum and sumsq as appropriate */
	   if(!parms[0].BOOLEAN) {
	      sum/=(float)no_comparisons;
	      sumsq/=(float)(no_comparisons*no_comparisons);
	      ave_length=(int)((float)(pxsize+pysize-2)/2);
	      mean=sum/(float)no_matrix;
	      sd=sqrt( (sumsq-(sum*sum)/(float)no_matrix) / (float)(no_matrix-1) );

	      if(!parms[0].STATS) { /* correct values if required */
	       if((cl.a.number+cl.b.number)>2) {
	         pairwise_mean=exp( parms[0].NA*log((float)ave_length)+parms[0].NB );
	         pairwise_sd=exp( parms[0].NASD*log((float)ave_length)+parms[0].NBSD );
	         mean=(float)parms[0].PRECISION*parms[0].NMEAN*(mean/pairwise_mean);
	         sd=(float)parms[0].PRECISION*parms[0].NSD*(sd/pairwise_sd);
	       } else {
	         mean=parms[0].NMEAN*(float)parms[0].PRECISION;
	         sd=parms[0].NSD*(float)parms[0].PRECISION;
	       }
	      }
	      for(l=0; l<strlen(domain[indx].align); ++l) 
	        for(m=0; m<strlen(domain[indy].align); ++m) 
	    	  prob[l+1][m+1]=(int)((float)parms[0].PRECISION*((float)prob[l+1][m+1]-mean)/sd );
	   }
	   if(treepath(domain,ndomain,cl,R2,V2,prob,score,rms,length,nfit,Pij,Dij,dist,Pijp,mean,sd,puse,touse,parms)==-1) return -1;

	   diff=100.0*fabs((*score)-oldscore)/(*score);
	   fprintf(parms[0].LOG,"iteration: %3d, ",iter+1);
	   fprintf(parms[0].LOG,"RMS: %7.3f, ",(*rms));
	   fprintf(parms[0].LOG," Sc diff: %6.2f %%, Sc: %7.3f, len: %4d, nfit: %4d\n",diff,*score,*length,*nfit);
	   iter++;
	   scorerise=((*score)-oldscore>0);
	   oldscore=(*score);
	   /* transform cluster b after the current transformation */
	   for(i=0; i<cl.b.number; ++i) {
	      indy=cl.b.member[i];
	      matmult(R2,V2,domain[indy].coords,domain[indy].ncoords,parms[0].PRECISION);
	      update(R2,domain[indy].r,V2,domain[indy].v);
	   }

	}
	if(diff<parms[0].SCORETOL) fprintf(parms[0].LOG,"Convergence ");
	else fprintf(parms[0].LOG,"No convergence ");
	fprintf(parms[0].LOG," after %d iterations\n",iter);
        if(parms[0].TREEPLOT) 
           probplot(prob,pxsize-1,pysize-1,1,
 	      (int)((float)parms[0].PRECISION*parms[0].CUTOFF)*(!parms[0].BOOLEAN)+(!parms[0].BOOLEAN),
	       parms[0].LOG);
	if(parms[0].TREEALIGN && align) {
 	 /* final alignment */
	 for(i=0; i<cl.a.number; ++i) {
	   pseq1[i]=domain[cl.a.member[i]].align;
	   psec1[i]=(char*)malloc((strlen(domain[cl.a.member[i]].align)+1)*sizeof(char));
	   k=0; 
	   for(j=0; j<strlen(domain[cl.a.member[i]].align); ++j) {
	      if(domain[cl.a.member[i]].align[j]!=' ') 
		 psec1[i][j]=domain[cl.a.member[i]].sec[k++];
	      else
		 psec1[i][j]=' ';
	   }
	   psec1[i][j]='\0';
	 }
	 for(i=0; i<cl.b.number; ++i) {
	   pseq2[i]=domain[cl.b.member[i]].align;
	   psec2[i]=(char*)malloc((strlen(domain[cl.b.member[i]].align)+1)*sizeof(char));
	   k=0;  
	   for(j=0; j<strlen(domain[cl.b.member[i]].align); ++j) { 
              if(domain[cl.b.member[i]].align[j]!=' ') 
	         psec2[i][j]=domain[cl.b.member[i]].sec[k++]; 
	      else 
		 psec2[i][j]=' ';
	   } 
	   psec2[i][j]='\0'; 
	 }
	 if(strcmp(parms[0].logfile,"silent")!=0) {
	   display_align(pseq1,cl.a.number,pseq2,cl.b.number,psec1,psec2,touse,puse,parms[0].COLUMNS,1,1,parms[0].LOG);
	 } else {
	   display_align(pseq1,cl.a.number,pseq2,cl.b.number,psec1,psec2,touse,puse,parms[0].COLUMNS,1,1,stdout);
	 }
	 for(i=0; i<cl.a.number; ++i) free(psec1[i]);
	 for(i=0; i<cl.b.number; ++i) free(psec2[i]);
	}

	/* copying all the "align"s into "oldaligns" for the next cluster 
	 *  that is unless we are in the first of two sets of parameters */
	if(!rev) {
	  for(i=0; i<cl.a.number; ++i) {
	     indx=cl.a.member[i];
	     strcpy(domain[indx].oldalign,domain[indx].align);
	  }
	  for(i=0; i<cl.b.number; ++i) {
	     indy=cl.b.member[i];
	     strcpy(domain[indy].oldalign,domain[indy].align);
	  }
	}

	/* freeing memory */
        for(j=0; j<(pxsize+2); ++j) 
	   free(prob[j]); 
	free(prob);

	/* free-ing to keep purify happy */
	for(i=0; i<3; ++i) 
	  free(R2[i]);
	free(R2);
	free(V2);
	free(puse); 
	free(touse);
	free(xcount); free(ycount);
	free(pseq1); free(pseq2); free(psec1); free(psec2);
}
