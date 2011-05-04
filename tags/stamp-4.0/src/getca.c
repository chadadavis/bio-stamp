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

 The WORK is Copyright (1995) University of Oxford
	Administrative Offices
	Wellington Square
	Oxford OX1 2JD U.K.

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include "include.h"

int getca(IN,coords,aa,numb,ncoord,start,end,type,MAXats,REVERSE,seqonly)
FILE *IN;
float **coords;
char *aa;
struct brookn *numb;
int *ncoord;
struct brookn start,end;
int type; 	/* 1 = all CA atoms in the file, 2 = single chain, 3 = specific start and end */
int MAXats;
int REVERSE;	/* if 1, then reverse the order of the data */
int seqonly;	/* 1 => just return the sequence */
{
	int i,j,k;
	int begin;
	int number;

	char cid,ins,alt;

	char *buff;
	char tmp[10];

	buff=(char*)malloc(100*sizeof(char));
	begin=0;
	(*ncoord)=0;

	while((buff=fgets(buff,99,IN))!=NULL) {
	   if((strncmp(buff,"ENDMDL",6)==0 || strcmp(buff,"END   ",6)==0) && begin==1) {
		break; 
	   }
	   if(strncmp(buff,"ATOM  ",6)==0 && strncmp(&buff[12]," CA ",4)==0) {
	      /* get chain, number and insertion code */
	      cid=buff[21];
	      sscanf(&buff[22],"%d",&number);
	      ins=buff[26];
	      alt=buff[16]; /* alternate position indicator */
	      if(!begin && 
		 ((start.cid==cid && start.n==number && start.in==ins) ||
		  (start.cid==cid && type==2) ||
		  (type==1) )) begin=1;
	      if(begin && type==2 && cid!=start.cid) 
		      break; 
	      if(begin && (alt==' ' || alt=='A' || alt=='1')) {
		if(!seqonly) {
		   coords[(*ncoord)]=(float*)malloc(3*sizeof(float));
		   for(i=0; i<3; ++i) {
		     strncpy(&tmp[0],&buff[30+i*8],8); 
		     tmp[8]='\0'; 
		     sscanf(&buff[30+i*8],"%f",&coords[(*ncoord)][i]);
		   }
		   numb[(*ncoord)].cid=cid; numb[(*ncoord)].in=ins; numb[(*ncoord)].n=number;
		 }
		 aa[(*ncoord)]=a3to1(&buff[17]);
		 (*ncoord)++;
		 if((*ncoord)>MAXats) {
		    fprintf(stderr,"error: number of coordinates read surpasses memory limit\n");
		    return -1;
		 }
		/* aa[(*ncoord)]=' ';*/
	      } /* end of if(begin... */
	     if(begin &&  end.cid==cid && end.n==number && end.in==ins && type==3) 
				break;
		/* this residing after the last "if" makes the set of atoms inclusive */
	   } /* end of if(strncmp(buff,"ATOM... */
	} /* end of while((buff=... */
	aa[(*ncoord)]='\0';
	free(buff);
	if(!begin) {
	   fprintf(stderr,"error: begin of sequence not found in PDB file\n");
	   (*ncoord)=0;
	   return -1;
	} else
	   return 0;
}
