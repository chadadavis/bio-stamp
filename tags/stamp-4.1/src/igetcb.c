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
#include <math.h>
#include <stamp.h>
#include <igetcb.h>
#define PI 3.141592653589793

/* Reads in C beta coordinates as integers (like igetca) 
 *
 * If a glycine is encountered, a "ghost" CB is built around an
 *  ideal geometry (see RBR_build_cb()) */

float *RBR_build_cb(float *GLY_N,
		    float *GLY_CA,
		    float *GLY_CO,
		    float angle,
		    float distance,
		    FILE *OUT
);


int igetcb(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
	struct brookn start, struct brookn end, int type, int MAXats,
	int REVERSE, int PRECISION, FILE *OUTPUT) {

	int i,j,k;
	int begin;
	int number;
	int got_gly_ca,got_gly_co,got_gly_n;
	int *ccoord;


	char cid,in;
	char alt;
	char tmp[10];
	char *buff,*add_buff;
	char caa;

	float x;
	float angle,distance;

	float *GLY_N,*GLY_CA,*GLY_CO;
	float *GLY_CB;

	struct brookn cnumb;



	angle=54.0*PI/180;
	distance=1.54;
	buff=(char*)malloc(100*sizeof(char));
	add_buff=buff;
	ccoord=(int*)malloc(3*sizeof(int));
	GLY_N=(float*)malloc(3*sizeof(float));
	GLY_CA=(float*)malloc(3*sizeof(float));
	GLY_CO=(float*)malloc(3*sizeof(float));

	begin=0;
	(*ncoord)=0;
	got_gly_ca=got_gly_co=got_gly_n=0;

	while((buff=fgets(buff,99,IN))!=NULL) {
	   if((strncmp(buff,"ENDMDL",6)==0 || strncmp(buff,"END   ",6)==0) && begin==1) {
                break;
           }
	   if(strncmp(buff,"ATOM  ",6)==0 && strncmp(&buff[17],"ACE",3)!=0 && strncmp(&buff[17],"FOR",3)!=0) {
	    alt=buff[16]; /* alternate position indicator */
	    /* to make ghost CB for glycine we must read in the N, CA and C coordinates */
	    if(strncmp(&buff[17],"GLY",3)==0) {
	       if(strncmp(&buff[12]," CA ",4)==0 && (alt==' ' || alt=='A' || alt=='1')) {
		  for(i=0; i<3; ++i) {
		    strncpy(&tmp[0],&buff[30+i*8],8);
		    tmp[8]='\0';
		    sscanf(&buff[30+i*8],"%f",&x);
		    GLY_CA[i]=x;
		  }
		  got_gly_ca=1;
/*		  printf("%s",buff); */
		}
		if(strncmp(&buff[12]," C  ",4)==0 && (alt==' ' || alt=='A' || alt=='1')) {
		   for(i=0; i<3; ++i) {
		     strncpy(&tmp[0],&buff[30+i*8],8);
		     tmp[8]='\0';
		     sscanf(&buff[30+i*8],"%f",&x);
		     GLY_CO[i]=x;
		   }
		   got_gly_co=1;
/*		   printf("%s",buff); */
		}
		if(strncmp(&buff[12]," N  ",4)==0 && (alt==' ' || alt=='A' || alt=='1')) {
		   for(i=0; i<3; ++i) { 
		      strncpy(&tmp[0],&buff[30+i*8],8);
		      tmp[8]='\0';
		      sscanf(&buff[30+i*8],"%f",&x);
		      GLY_N[i]=x;
		   }
		   got_gly_n=1;
/*		   printf("%s",buff); */
		}
	    } else {
	      got_gly_ca=got_gly_co=got_gly_n=0;
	    }
	    if(strncmp(&buff[12]," CB ",4)==0 || (strncmp(&buff[17],"GLY",3)==0 && got_gly_n && got_gly_ca && got_gly_co)) {
	      /* get chain, number and insertion code */
	      cid=buff[21];
	      sscanf(&buff[22],"%d",&number);
	      in=buff[26];
	      if(!begin && 
		 ((start.cid==cid && start.n==number && start.in==in) ||
		  (start.cid==cid && type==2) ||
		  (type==1) )) begin=1;
	      if(begin && type==2 && start.cid!=cid) break;
	      if(begin && strncmp(&buff[17],"GLY",3)==0 && got_gly_n && got_gly_ca && got_gly_co) { 
	        /* if we have all the coordinates, build the CB for glycine */ 
		coords[(*ncoord)]=(int*)malloc(3*sizeof(int));
		 if((GLY_CB=RBR_build_cb(GLY_N,GLY_CA,GLY_CO,angle,distance,stdout))==NULL) {
		    return -1;
		  }
		  for(i=0; i<3; ++i) 
		     coords[(*ncoord)][i]=(int)(PRECISION*GLY_CB[i]);
		  got_gly_n = got_gly_ca = got_gly_co = 0;
		  aa[(*ncoord)]=a3to1(&buff[17]);
		  if(cid==' ') numb[(*ncoord)].cid='_';
		  else numb[(*ncoord)].cid=cid;
		  if(in==' ') numb[(*ncoord)].in='_';
	 	  else numb[(*ncoord)].in=in;
		  numb[(*ncoord)].n=number;
		  (*ncoord)++;
	       } else if(begin && alt==' ' || alt=='A' || alt=='1') {
		 coords[(*ncoord)]=(int*)malloc(3*sizeof(int));
		 /* only reads in the first position if more than one are given */
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&x);
		   coords[(*ncoord)][i]=(int)(PRECISION*x);
	          }
/*		  printf("CB %4d: %8d %8d %8d\n",
		     (*ncoord),coords[(*ncoord)][0],coords[(*ncoord)][1],coords[(*ncoord)][2]); */
		  aa[(*ncoord)]=a3to1(&buff[17]);
		  if(cid==' ') numb[(*ncoord)].cid='_'; 
		  else numb[(*ncoord)].cid=cid; 
		  if(in==' ') numb[(*ncoord)].in='_';
		  else numb[(*ncoord)].in=in; 
		  numb[(*ncoord)].n=number;
		  (*ncoord)++;
		}
		if((*ncoord)>MAXats) {
		    fprintf(stderr,"error: number of coordinates read surpasses memory limit\n");
		    return -1;
		}
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
	   free(GLY_N);
	   free(GLY_CA);
	   free(GLY_CO);
	   return 0;
	}
}

float *RBR_build_cb(
		    float *GLY_N,
		    float *GLY_CA,
		    float *GLY_CO,
		    float angle,
		    float distance,
		    FILE *OUT
)
{
	int i;
	float *V1,*V2;
	float *average;
	float *cross;
	float *cb,*GLY_CB;

	cb=(float*)malloc(3*sizeof(float));

	/* this sets the origin to the position of GLY_CA */


	V1=RBR_vector_diff(GLY_N,GLY_CA);
	V2=RBR_vector_diff(GLY_CO,GLY_CA);

	/* set lengths to unity */
	RBR_vector_unify(V1);
	RBR_vector_unify(V2);

	/* get average of V1 and V2 */
	average=RBR_vector_ave(V1,V2);
	/* reverse the sign to make it point in the opposite direction */
	for(i=0; i<3; ++i) average[i]=-1*average[i]; 
	RBR_vector_unify(average);

	/* cross product */
	cross=RBR_vector_cross(V1,V2); 
	RBR_vector_unify(cross);

	/* find coordinates based on ideal angle */
	for(i=0; i<3; ++i) {
	  cb[i]=cross[i]*(float)cos((double)angle)+average[i]*(float)cos((double)angle);
	}

	/* length of cb is unity, set to ideal */
	GLY_CB=RBR_vector_set_dist(cb,distance);

	/* this vector treats CA as the origin, we must
	 * get the "real" coordinate back */
	


	/* done, return the result */
	for(i=0; i<3; ++i) 
	  GLY_CB[i]=GLY_CB[i]+GLY_CA[i];

	free(cb);
	return GLY_CB;


}
