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

/* Takes a PDB file, a description of a domain, and a transformation.
 *  It applies the transformation to the coordinates outputs the 
 *  corresponding porition of the PDB file in PDB format */
int extract_pdb(IN,start,end,type,R,V,startats,HETERO,HOH,chainlabel,verbose,filename,OUT)
FILE *IN;
struct brookn start,end;
int type; /* 1 = all CA atoms in the file, 2 = single chain, 3 = specific start and end */
float **R;
float *V;
int startats;	/* 1 => write the transformation (ie. this is the first object in the file) */
int HETERO;	/* 1 => transform and include all HETATM regardless of their labels */
int HOH;	/* 1 => transform and include all waters */
char chainlabel;
int verbose;
char *filename;
FILE *OUT;
{
	int i,j,k;
	int begin,endnext;
	int number;
	int count;
	int found;

	char cid,in;

	char *buff;
	char tmp[10];

	float *coord;

	struct brookn old,this;

	coord=(float*)malloc(3*sizeof(float));
	buff=(char*)malloc(100*sizeof(char));
	begin=endnext=count=found=0;

	while(fgets(buff,99,IN)!=NULL) {
	   if(strncmp(buff,"ENDMDL",6)==0) break;
	   buff[strlen(buff)-2]='\0';
	   if(strncmp(buff,"ATOM  ",6)==0 || strncmp(buff,"HETATM",6)==0) {
	      /* get chain, number and insertion code */
	      this.cid=buff[21];
	      sscanf(&buff[22],"%d",&this.n);
	      this.in=buff[26];
	      if(endnext && type==3 && (this.cid!=old.cid || this.in!=old.in || this.n!=old.n)) 
		 begin=0;
	      if(!begin && 
		 ((start.cid==this.cid && start.n==this.n && start.in==this.in && type==3) ||
		  (start.cid==this.cid && type==2) ||
		  (type==1) )) {
		  begin=1; found=1;

	          if(startats) {
	            fprintf(OUT,"REMARK The following transformation has been applied to the coordinates, taken  \n");
	            fprintf(OUT,"REMARK  from the orginal PDB file: %40s     \n",filename);
	            for(i=0; i<3; ++i) {
		      fprintf(OUT,"REMARK  "); 
		      for(j=0; j<3; ++j) fprintf(OUT,"%12.5f ",R[i][j]);
		      fprintf(OUT,"    %12.5f",V[i]); 
		      for(j=0; j<17; ++j) fprintf(OUT," "); 
		      fprintf(OUT,"\n");
	            }
	            fprintf(OUT,"REMARK"); for(j=0; j<74; ++j) fprintf(OUT," ");
	            fprintf(OUT,"\n");
	          }
	      }
	      if(type==2 && start.cid!=this.cid && begin) break;
	      if((begin && strncmp(buff,"ATOM  ",6)==0) ||
                 (HETERO && strncmp(buff,"HETATM",6)==0 && 
		  (strncmp(&buff[17],"HOH",3)!=0 && strncmp(&buff[17],"WAT",3)!=0)) ||
		  (HOH==1 && (strncmp(&buff[17],"HOH",3)==0 || strncmp(&buff[17],"WAT",3)==0))) {
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&coord[i]);
		 }
		 buff[30]='\0';
		 /* transform coordinates */
/*		 if(count<10) 
		    printf("%d: %8.4f %8.4f %8.4f\n",count+1,coord[0],coord[1],coord[2]); */
		 matmult(R,V,&coord,1);
/*		 if(count<10)
		    printf("%d: %8.4f %8.4f %8.4f\n",count+1,coord[0],coord[1],coord[2]); */
		 count++;
		 if(chainlabel!='\0') buff[21]=chainlabel;
		 fprintf(OUT,"%s",buff);
 		 fprintf(OUT,"%8.3f%8.3f%8.3f",coord[0],coord[1],coord[2]); 
		 fprintf(OUT,"%s",&buff[54]);
		 for(i=0; i<(26-strlen(&buff[54])); ++i) fprintf(OUT," ");
		 fprintf(OUT,"\n");
	      } /* end of if(begin... */
	      if(begin && type==3 && end.cid==this.cid && end.n==this.n && end.in==this.in) 
	          endnext=1;
	        /* this residing after the last "if" makes the set of type==3 atoms inclusive */
	      old.cid=this.cid; old.n=this.n; old.in=this.in;
	   } else if(verbose==1 && startats==1 && strncmp(buff,"TER",3)!=0 && strncmp(buff,"SIGATM",6)!=0) {
		 fprintf(OUT,"%s",buff);
		 for(i=0; i<(80-strlen(buff)); ++i) fprintf(OUT," ");
		 fprintf(OUT,"\n");
	   } /* end of if(strncmp(buff,"ATOM... */
	} /* end of while((buff=... */
	free(buff);
	if(!found) {
	   printf("error: begin of sequence not found in PDB file\n");
	   printf("  file %s  region ",filename);	
	   if(type==1) printf(" all residues\n");
	   else if(type==2) printf(" chain %c\n",start.cid);
	   else printf(" %c%d%c to %c%d%c\n",start.cid,start.n,start.in,
		   end.cid,end.n,end.in);
	   return -1;
	} else
	   return 0;
}
