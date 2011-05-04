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
#include <stamp.h>

/* slightly varied version of getca
 * Finds any atom type (supplied by a four character string)
 * Coordinates are multiplied by 1000 and converted to integers */

int igetgen(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
	struct brookn start, struct brookn end, int type, char *atom,
	int MAXats, int REVERSE, int PRECISION, FILE *OUTPUT) {

	int i,j,k;
	int begin;
	int number;
	int *ccoord;

	char cid,in;
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

	while((buff=fgets(buff,99,IN))!=NULL) {
	   if((strncmp(buff,"ENDMDL",6)==0 || strncmp(buff,"END   ",6)==0) && begin==1) {
                break;
           }
	   if(strncmp(buff,"ATOM  ",6)==0 && strncmp(&buff[12],atom,4)==0) {
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
	      if(begin && (alt==' ' || alt=='A' || alt=='1') ) { 
		 /* only reads in the first position if more than one are given */
		 coords[(*ncoord)]=(int*)malloc(3*sizeof(int));
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&x);
		   coords[(*ncoord)][i]=(int)(PRECISION*x);
		 }
		 aa[(*ncoord)]=a3to1(&buff[17]);
		 
		 if(cid==' ') numb[(*ncoord)].cid='_'; 
		 else numb[(*ncoord)].cid=cid; 
		 if(in==' ') numb[(*ncoord)].in='_';
		 else numb[(*ncoord)].in=in; 
		 numb[(*ncoord)].n=number;
		 (*ncoord)++;
		 if((*ncoord)>MAXats) {
		    fprintf(stderr,"error: number of coordinates read surpasses memory limit\n");
		    return -1;
		 }
		 aa[(*ncoord)]=' ';
	      } /* end of if(begin... */
	      if(begin && end.cid==cid && end.n==number && end.in==in && type==3) 
		 break;
	      /* this residing after the last "if" makes the set of atoms inclusive */
	   } /* end of if(strncmp(buff,"ATOM... */
	} /* end of while((buff=... */
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
		ccoord[0]=coords[i][0]; ccoord[1]=coords[i][1]; ccoord[2]=coords[i][2];
		cnumb.cid=numb[i].cid; cnumb.n=numb[i].n; cnumb.in=numb[i].in;
		caa=aa[i];
		/* copy the C into the N-terminal */
		coords[i][0]=coords[(*ncoord)-i-1][0]; 
		coords[i][1]=coords[(*ncoord)-i-1][1];
		coords[i][2]=coords[(*ncoord)-i-1][2];
		numb[i].cid=numb[(*ncoord)-i-1].cid;
		numb[i].n=numb[(*ncoord)-i-1].n;
		numb[i].in=numb[(*ncoord)-i-1].in;
		aa[i]=aa[(*ncoord)-i-1];
		/* and visa versa */
		coords[(*ncoord)-i-1][0]=ccoord[0];
		coords[(*ncoord)-i-1][1]=ccoord[1];
		coords[(*ncoord)-i-1][2]=ccoord[2];
		numb[(*ncoord)-i-1].cid=cnumb.cid;
		numb[(*ncoord)-i-1].n=cnumb.n;
		numb[(*ncoord)-i-1].in=cnumb.in;
		aa[(*ncoord)-i-1]=caa;
	    }
	   }
	   free(ccoord);
	   return 0;
	}
}