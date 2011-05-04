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
#include <stdio.h>
#include "stamp.h"

/* This is igetca converted back to getca - there were some changes not introduced
 *  into wonky old getca that is only used in alignfit
 * Coordinates are multiplied by 1000 and converted to integers */

int getca(FILE *IN, float **coords, char *aa, struct brookn *numb, int *ncoord,
	struct brookn start, struct brookn end, int type, int MAXats,
	int REVERSE,  FILE *OUTPUT) {

	int i,j,k;
	int begin;
	int number,last_number;
	int seq_only;
	int *ccoord;

	char cid,in,last_cid,last_in;
	char alt;
	char tmp[10];
	char *buff,*add_buff;
	char caa;

	float x;

	struct brookn cnumb;

	buff=(char*)malloc(100*sizeof(char));
	add_buff=buff;
	ccoord=(int*)malloc(3*sizeof(int));

	begin=0;
	(*ncoord)=0;

	last_in = '?'; last_cid = '?'; last_number = 999999;
	if(coords==NULL) seq_only=1;
	else seq_only=0;
/*	printf("SEQONLY is %d\n",seq_only); */

	while((buff=fgets(buff,99,IN))!=NULL) {
	   if((strncmp(buff,"ENDMDL",6)==0 || strncmp(buff,"END   ",6)==0) && begin==1) {
                break;
           }
	   if(strncmp(buff,"ATOM  ",6)==0 && strncmp(&buff[12]," CA ",4)==0) {
	      /* get chain, number and insertion code */
	      cid=buff[21];
	      sscanf(&buff[22],"%d",&number);
	      in=buff[26];
	      alt=buff[16]; /* alternate position indicator */
	      if(!begin && 
		 ((start.cid==cid && start.n==number && start.in==in) ||
		  (start.cid==cid && type==2) ||
		  (type==1) )) begin=1;
	      if(begin && type==2 && start.cid!=cid) break;
/* SMJS Changed to be like Robs version */
	      if(begin && (alt==' ' || alt=='A' || alt=='1' || alt=='L' || alt=='O') && 
		  !(cid==last_cid && number==last_number && in==last_in)) { 
		 /* only reads in the first position if more than one are given */
		if(seq_only==0) {
		 coords[(*ncoord)]=(float*)malloc(3*sizeof(float));
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&x);
		   coords[(*ncoord)][i]=x;
		 }
	        }
		aa[(*ncoord)]=a3to1(&buff[17]);
	        if(seq_only==0) {	 
		  if(cid==' ') numb[(*ncoord)].cid='_'; 
		  else numb[(*ncoord)].cid=cid; 
		  if(in==' ') numb[(*ncoord)].in='_';
		  else numb[(*ncoord)].in=in; 
		  numb[(*ncoord)].n=number;
		}

		(*ncoord)++;

		if((*ncoord)>MAXats) {
		    fprintf(stderr,"error: number of coordinates read surpasses memory limit\n");
		    return -1;
		 }
		 aa[(*ncoord)]=' ';
		 last_cid = cid; last_in = in; last_number = number;

		 if(begin && end.cid==cid && end.n==number && end.in==in && type==3) 
		     break;
		 /* this residing after the last "if" makes the set of atoms inclusive */
	      } 
	      /* in case of missing CA in the last residue */
	      if(begin && cid==end.cid && number>end.n && type==3) break;  
	   } 
	} 
	aa[(*ncoord)]='\0';
	free(add_buff);
	if(!begin) {
	   fprintf(stderr,"error: begin of sequence not found in PDB file\n");
	   (*ncoord)=0;
	   free(ccoord);
	   return -1;
	} else {
	  /* reverse the data if necessary */
	  if(REVERSE) {
	     j=(int)((*ncoord)/2);
	     for(i=0; i<j; ++i) {
		/* save the N-terminal data */
		if(seq_only==0) {
		   ccoord[0]=coords[i][0]; ccoord[1]=coords[i][1]; ccoord[2]=coords[i][2];
	   	   cnumb.cid=numb[i].cid; cnumb.n=numb[i].n; cnumb.in=numb[i].in;
		}
		caa=aa[i];
		/* copy the C into the N-terminal */
		if(seq_only==0) {
		   coords[i][0]=coords[(*ncoord)-i-1][0]; 
		   coords[i][1]=coords[(*ncoord)-i-1][1];
		   coords[i][2]=coords[(*ncoord)-i-1][2];
		   numb[i].cid=numb[(*ncoord)-i-1].cid;
		   numb[i].n=numb[(*ncoord)-i-1].n;
		   numb[i].in=numb[(*ncoord)-i-1].in;
		}
		aa[i]=aa[(*ncoord)-i-1];
		/* and visa versa */
		if(seq_only==0) {
		   coords[(*ncoord)-i-1][0]=ccoord[0];
		   coords[(*ncoord)-i-1][1]=ccoord[1];
		   coords[(*ncoord)-i-1][2]=ccoord[2];
		   numb[(*ncoord)-i-1].cid=cnumb.cid;
		   numb[(*ncoord)-i-1].n=cnumb.n;
		   numb[(*ncoord)-i-1].in=cnumb.in;
		}
		aa[(*ncoord)-i-1]=caa;
	    }
	   }
	   free(ccoord);
	   return 0;
	}
}
