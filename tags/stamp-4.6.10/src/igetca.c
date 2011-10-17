#include <stdio.h>
#include <stamp.h>

/* slightly varied version of getca.  
 * Coordinates are multiplied by 1000 and converted to integers */

/* April 2007 - change to allow:
 *  1. slop on the N-/C- terminus */

int igetca(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
	   struct brookn start, struct brookn end, int type, int MAXats,
	   int REVERSE, int PRECISION, FILE *OUTPUT) {

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
    if(buff == NULL) {
	/* this should probably go to stderr... but sending to OUTPUT to conform to other messages */
	fprintf(OUTPUT, "error: malloc failed on buff\n");
	return -1;
    }
    add_buff=buff;
    ccoord=(int*)malloc(3*sizeof(int));
    if(ccoord == NULL) {
	/* this should probably go to stderr... but sending to OUTPUT to conform to other messages */
	fprintf(OUTPUT, "error: malloc failed on ccoord\n");
	return -1;
    }

    begin=0;
    (*ncoord)=0;

    last_in = '?'; last_cid = '?'; last_number = 999999;
    if(coords==NULL) seq_only=1;
    else seq_only=0;

    while((buff=fgets(buff,99,IN))!=NULL) {
	if((strncmp(buff,"ENDMDL",6)==0 || strncmp(buff,"END   ",6)==0) && begin==1) {
	    break;
	}
	if((strncmp(buff,"ATOM  ",6)==0) || (strncmp(buff,"HETATM",6)==0)) {
            /* look at non-CA to see if we are in/out of range */
	    /* get chain, number and insertion code */
	    cid=ltou(buff[21]);
	    sscanf(&buff[22],"%d",&number);
	    in=buff[26];
	    alt=buff[16]; /* alternate position indicator */
	    if(!begin && 
	       ((start.cid==cid && start.n==number && start.in==in) ||
		(start.cid==cid && type==2) ||
		(type==1) )) {
		//printf("Began at %c %4d %c vs %c %4d %c (type %d)\n",cid,number,in,start.cid,start.n,start.in,type);
		begin=1;
            }
	    if(begin && type==2 && start.cid!=cid) break;
            if(strncmp(&buff[12]," CA ",4)==0) {
		/* SMJS Changed to be like Robs version */
		if(begin && (alt==' ' || alt=='A' || alt=='1' || alt=='L' || alt=='O') && 
		   !(cid==last_cid && number==last_number && in==last_in)) { 
		    /* only reads in the first position if more than one are given */
		    //printf("%s",buff); 
		    if(seq_only==0) {
			coords[(*ncoord)]=(int*)malloc(3*sizeof(int));
			if(coords[(*ncoord)] == NULL) {
			    /* this should probably go to stderr... but sending to OUTPUT to conform to other messages */
			    fprintf(OUTPUT, "error: malloc failed on coords[(*ncoord)]\n");
			    return -1;
			}
			for(i=0; i<3; ++i) {
			    strncpy(&tmp[0],&buff[30+i*8],8); 
			    tmp[8]='\0'; 
			    sscanf(&buff[30+i*8],"%f",&x);
			    coords[(*ncoord)][i]=(int)(PRECISION*x);
			}
		    }
		    aa[(*ncoord)]=a3to1(&buff[17]);
		    if(seq_only==0) {	 
			if(cid==' ') {
			    numb[(*ncoord)].cid='_';
			}
			else {
			    numb[(*ncoord)].cid=cid;
			}

			if(in==' ') {
			    numb[(*ncoord)].in='_';
			}
			else {
			    numb[(*ncoord)].in=in; 
			}
			numb[(*ncoord)].n=number;
		    }

		    last_cid = cid; last_in = in; last_number = number;
		    (*ncoord)++;

		    if((*ncoord)>MAXats) {
			fprintf(OUTPUT,"error: number of coordinates read surpasses memory limit %d\n",MAXats);
			return -1;
		    }
		    aa[(*ncoord)]=' ';
		} 
		if(begin && end.cid==cid && end.n==number && end.in==in && type==3) break;
		/* this residing after the last "if" makes the set of atoms inclusive */
	    } 
            /* in case of missing CA in the last residue */
	    if(begin && cid==end.cid && number>end.n && type==3) break;  
	} 
    } 
    aa[(*ncoord)]='\0';
    free(add_buff);

    if(!begin) {
	fprintf(OUTPUT,"error: begin of sequence not found in PDB file\n");
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
