#include <stdio.h>
#include <math.h>
#include "stamp.h"
#define PI 3.141592653589793
#define MAX_NATS 50


int res_name_comp (struct brookn name1, struct brookn name2);

/* Reads in side chain atoms a la igetca 
 * RBR January 1996
 */
float *RBR_build_cb(float *GLY_N,
		    float *GLY_CA,
		    float *GLY_CO,
		    float angle,
		    float distance,
		    FILE *OUT
);


int igetside(FILE *IN, struct side_chain *side, char *aa, struct brookn *numb, int *nres,
	struct brookn start, struct brookn end, int type, int MAXats,
	int REVERSE, int PRECISION, FILE *OUTPUT) {

	int i,j,k;
	int begin;
	int got_gly_ca,got_gly_co,got_gly_n;
	int *ccoord;


	char alt;
	char tmp[10];
	char *buff,*add_buff;
	char caa;

	float x;
	float angle,distance;
	float *GLY_N,*GLY_CA,*GLY_CO;
        float *GLY_CB;



	struct brookn cnumb;
	struct brookn this,last;


	angle=54.0*PI/180;
	distance=1.54;
	buff=(char*)malloc(100*sizeof(char));
	add_buff=buff;
	ccoord=(int*)malloc(3*sizeof(int));

	begin=0;
	(*nres)=0;

	while((buff=fgets(buff,99,IN))!=NULL) {
	   if((strncmp(buff,"ENDMDL",6)==0 || strncmp(buff,"END   ",6)==0) && begin==1) {
                break;
           }
	   if(strncmp(buff,"ATOM  ",6)==0 && strncmp(&buff[17],"ACE",3)!=0 && strncmp(&buff[17],"FOR",3)!=0) {
	      alt=buff[16]; /* alternate position indicator */
	      /* get chain, number and insertion code */
	      this.cid=buff[21];
	      sscanf(&buff[22],"%d",&this.n);
	      this.in=buff[26];
	      if(!begin && 
		 ((start.cid==this.cid && start.n==this.n && start.in==this.in) ||
		  (start.cid==this.cid && type==2) ||
		  (type==1) )) begin=1;
	      if(begin && type==2 && start.cid!=this.cid) break;
	      /* work out if we are in a new residue */
	      if(res_name_comp(this,last)!=0) { /* New residue */
		/*
		struct side_chain {
		        int n;
		        int **names;
		        int **coords;
		};
		*/
		   last.cid=this.cid; last.n=this.n; last.in=this.in;
                   side[(*nres)].n=0;
		   /* We'll skip any attempt to be clever with memory for
		    * the moment */
                   side[(*nres)].coords=(int**)malloc(MAX_NATS*sizeof(int *));
		   side[(*nres)].names=(char**)malloc(MAX_NATS*sizeof(char*));
		   for(i=0; i<MAX_NATS; ++i) {
			side[(*nres)].coords[i]=(int*)malloc(3*sizeof(int));
			side[(*nres)].names[i]=(char*)malloc(6*sizeof(char));
		   }
	      }
	      if(begin && strncmp(&buff[17],"GLY",3)==0 && got_gly_n && got_gly_ca && got_gly_co) { 
	        /* if we have all the coordinates, build the CB for glycine */ 
		 if((GLY_CB=RBR_build_cb(GLY_N,GLY_CA,GLY_CO,angle,distance,stdout))==NULL) {
		    return -1;
		  }
		  for(i=0; i<3; ++i) 
		     side[(*nres)].coords[side[(*nres)].n][i]=(int)(PRECISION*GLY_CB[i]);
		  got_gly_n = got_gly_ca = got_gly_co = 0;
		  aa[(*nres)]=a3to1(&buff[17]);
		  if(this.cid==' ') numb[(*nres)].cid='_';
		  else numb[(*nres)].cid=this.cid;
		  if(this.in==' ') numb[(*nres)].in='_';
	 	  else numb[(*nres)].in=this.in;
		  numb[(*nres)].n=this.n;
		  (*nres)++;
/* SMJS Changed to be like Robs version */
/* SMJS Changed to add () around ||s */
	       } else if(begin && (alt==' ' || alt=='A' || alt=='1' || alt=='L' || alt=='O')) {
		 /* only reads in the first position if more than one are given */
		 for(i=0; i<3; ++i) {
		   strncpy(&tmp[0],&buff[30+i*8],8); 
		   tmp[8]='\0'; 
		   sscanf(&buff[30+i*8],"%f",&x);
		   side[(*nres)].coords[side[(*nres)].n][i]=(int)(PRECISION*x);
	          }
		  strncpy(side[(*nres)].names[side[(*nres)].n],&buff[12],4);
		  
		  printf(" %4d: %5s %8d %8d %8d\n",
		     (*nres),side[(*nres)].names[side[(*nres)].n],side[(*nres)].coords[side[(*nres)].n][0],side[(*nres)].coords[side[(*nres)].n][1],side[(*nres)].coords[side[(*nres)].n][2]);
		  aa[(*nres)]=a3to1(&buff[17]);
		  if(this.cid==' ') numb[(*nres)].cid='_'; 
		  else numb[(*nres)].cid=this.cid; 
		  if(this.in==' ') numb[(*nres)].in='_';
		  else numb[(*nres)].in=this.in; 
		  numb[(*nres)].n=this.n;
		  (*nres)++;
		}
		if((*nres)>MAXats) {
		    fprintf(stderr,"error: number of coordinates read surpasses memory limit\n");
		    return -1;
	        }
	      if(begin && end.cid==this.cid && end.n==this.n && end.in==this.in && type==3) 
		 break;
	      /* this residing after the last "if" makes the set of atoms inclusive */
	    } /* end of if(strncmp(buff,"ATOM... */
	} /* end of while((buff=... */
	aa[(*nres)]='\0';
	free(add_buff);
	if(!begin) {
	   fprintf(stderr,"error: begin of sequence not found in PDB file\n");
	   (*nres)=0;
	   free(ccoord);
	   return -1;
	} else {
	  /* reverse the data if necessary */
/*	  if(REVERSE) {
	     j=(int)((*nres)/2);
	     for(i=0; i<j; ++i) {
		ccoord[0]=coords[i][0]; ccoord[1]=coords[i][1]; ccoord[2]=coords[i][2];
		cnumb.cid=numb[i].cid; cnumb.n=numb[i].n; cnumb.in=numb[i].in;
		caa=aa[i];
		coords[i][0]=coords[(*nres)-i-1][0]; 
		coords[i][1]=coords[(*nres)-i-1][1];
		coords[i][2]=coords[(*nres)-i-1][2];
		numb[i].cid=numb[(*nres)-i-1].cid;
		numb[i].n=numb[(*nres)-i-1].n;
		numb[i].in=numb[(*nres)-i-1].in;
		aa[i]=aa[(*nres)-i-1];
		coords[(*nres)-i-1][0]=ccoord[0];
		coords[(*nres)-i-1][1]=ccoord[1];
		coords[(*nres)-i-1][2]=ccoord[2];
		numb[(*nres)-i-1].cid=cnumb.cid;
		numb[(*nres)-i-1].n=cnumb.n;
		numb[(*nres)-i-1].in=cnumb.in;
		aa[(*nres)-i-1]=caa;
	    }
	   }
*/
	   free(ccoord);
	   return 0;
	}
}


int res_name_comp (struct brookn name1, struct brookn name2) {

	if(name1.cid==name2.cid && name1.n==name2.n && name1.in==name2.in) return 0;
	else return 1;

}
