/******************************************************************************
 The computer software and associated documentation called STAMP hereinafter
 referred to as the WORK which is more particularly identified and described in 
 the LICENSE.  Conditions and restrictions for use of
 this package are also in the LICENSE.

 The WORK is only available to licensed institutions.

 The WORK was developed by: 
	Robert B. Russell and Geoffrey J. Barton

 Of current contact addresses:

 Robert B. Russell (RBR)             Geoffrey J. Barton (GJB)
 Bioinformatics                      EMBL-European Bioinformatics Institute
 SmithKline Beecham Pharmaceuticals  Wellcome Trust Genome Campus
 New Frontiers Science Park (North)  Hinxton, Cambridge, CB10 1SD U.K.
 Harlow, Essex, CM19 5AW, U.K.       
 Tel: +44 1279 622 884               Tel: +44 1223 494 414
 FAX: +44 1279 622 200               FAX: +44 1223 494 468
 e-mail: russelr1@mh.uk.sbphrd.com   e-mail geoff@ebi.ac.uk
                                     WWW: http://barton.ebi.ac.uk/

   The WORK is Copyright (1997,1998,1999) Robert B. Russell & Geoffrey J. Barton
	
	
	

 All use of the WORK must cite: 
 R.B. Russell and G.J. Barton, "Multiple Protein Sequence Alignment From Tertiary
  Structure Comparison: Assignment of Global and Residue Confidence Levels",
  PROTEINS: Structure, Function, and Genetics, 14:309--323 (1992).
*****************************************************************************/
#include <stdio.h>
#include <stamp.h>

/* Takes a PDB file, a description of a domain, and a transformation.
 *  It applies the transformation to the coordinates outputs the 
 *  corresponding porition of the PDB file in PDB format */

int extract_pdb(FILE *IN, struct brookn start, struct brookn end, int type, 
	float **R, float *V, int startats, int HETERO, int NUCLEIC, int HOH,
	char chainlabel, int verbose, char *filename, FILE *OUT) {

	int i,j,k;
	int begin,ended,endnext;
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

	ended=0; begin=0;
	while(fgets(buff,99,IN)!=NULL) {
	   if(strncmp(buff,"ENDMDL",6)==0) ended=1;
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
	      if(type==2 && start.cid!=this.cid && begin) ended=1;
	      if((!ended && begin && strncmp(buff,"ATOM  ",6)==0) ||
                 (HETERO && strncmp(buff,"HETATM",6)==0 && 
		  (strncmp(&buff[17],"HOH",3)!=0 && strncmp(&buff[17],"WAT",3)!=0 && strncmp(&buff[17],"DOD",3)!=0)) ||
		 (NUCLEIC && strncmp(buff,"ATOM  ",6)==0 &&
                        (strncmp(&buff[17],"  A",3)==0 || strncmp(&buff[17],"  G",3)==0 ||
                         strncmp(&buff[17],"  T",3)==0 || strncmp(&buff[17],"  C",3)==0 ||
                         strncmp(&buff[17],"  U",3)==0)) ||
		  (HOH==1 && (strncmp(&buff[17],"HOH",3)==0 || strncmp(&buff[17],"WAT",3)==0 || strncmp(&buff[17],"DOD",3)==0))) {
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&coord[i]);
		 }
		 buff[30]='\0';
		 /* transform coordinates */
/*		 if(count<10) 
		    printf("%d: %8.4f %8.4f %8.4f\n",count+1,coord[0],coord[1],coord[2]); */
		 fmatmult(R,V,&coord,1);
/*		 if(count<10)
		    printf("%d: %8.4f %8.4f %8.4f\n",count+1,coord[0],coord[1],coord[2]); */
		 count++;
		 if(chainlabel!='\0') buff[21]=chainlabel;
		 fprintf(OUT,"%s",buff);
 		 fprintf(OUT,"%8.3f%8.3f%8.3f",coord[0],coord[1],coord[2]); 
		 fprintf(OUT,"%s",&buff[54]);
		 for(i=0; i<(26-strlen(&buff[54])); ++i) fprintf(OUT," ");
		 fprintf(OUT,"\n");
	      } 
	      if(begin && type==3 && end.cid==this.cid && end.n==this.n && end.in==this.in) 
	          endnext=1;
	        /* this residing after the last "if" makes the set of type==3 atoms inclusive */
	      old.cid=this.cid; old.n=this.n; old.in=this.in;
	   } else if(!ended && verbose==1 && startats==1 && strncmp(buff,"TER   ",6)!=0 && strncmp(buff,"SIGATM",6)!=0) {
		 fprintf(OUT,"%s",buff);
		 for(i=0; i<(80-strlen(buff)); ++i) fprintf(OUT," ");
		 fprintf(OUT,"\n");
	   } 
	} 
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
