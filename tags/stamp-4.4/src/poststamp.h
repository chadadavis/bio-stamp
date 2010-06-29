#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
};

struct domain_loc{		/* This structure allows rather complex domains to be described */
   char filename[4096];
   char id[100];
   int nobj;			/* The number of objects considered within the named file */
   int *type;			/* The type that each object is:
					0 ==> an error
					1 ==> All of the residues in the file
					2 ==> A particular chain
					3 ==> A particular named region (eg. A 25 _ to B 10 _ ) */
   struct brookn *start;	/* There will be a start and end for each 'object' */
   struct brookn *end;
   int *reverse;
   int **coords;
   char *aa;
   char *align;
   char *oldalign;
   char *sec;
   struct brookn *numb;
   int ncoords;
   int *use;
   float **R;			/* Initial transformation */
   float *V;
   float **r;			/* current transformation, when STAMP is done, we must update(r,v,R,V) to get the final transformation */
   float *v;
   float value;
   };

/* structure for storing protein sequence data */
struct seqdat {	/* all lengths include char terminator and [0] */
    int ilen;	/* length of identifier*/
    char *id;	/* identifier */
    int tlen;	/* length of title */
    char *title;	/* title */
    int slen;	/* length of sequence*/
    char *seq;	/* sequence */
};
/*  mseq Structure for storing multiply aligned sequence data - 
    arranged so that
    aligned positions in different sequences are adjacent 
    This structure could easily be extended to include additional annotations
    (eg. the name of the alignment, its history, etc. etc.) can be added
    without having to modify existing code (I hope..)
    */

struct stampdat{
   char what;
   char *title;
   float *n;
};

int *getstampsecsum();

int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos);

int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window);

float rossmann(int **atoms1, int **atoms2, int start, int end,
        float const1, float const2, float *Dij, float *Cij, int PRECISION);


int count_domain(FILE *IN);
int matmult(float **r, float *v, int **coord, int n, int PRECISION);
int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans);
void exit_error();
int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, 
        int *gottrans, char *env, int DSSP, FILE *OUTPUT);
int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);
int igetca(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type, int MAXats,
        int REVERSE, int PRECISION, FILE *OUTPUT);

void closefile(FILE *handle,char *filename);
FILE *openfile(char *filename, char *type);
