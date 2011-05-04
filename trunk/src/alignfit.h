#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXslen  10000
#define MAXnbloc 10000
#define MAXtlen  200

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
   };

struct domain_loc {		/* This structure allows rather complex domains to be described */
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
   float **coords;
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

struct parameters {
   char bloc_file[4096];
   char dom_file[4096];
   char parm_file[4096];
   char prefix[4096];
   int MAX_SEQ_LEN;
   int PAIRWISE;
   int TREEWISE;
   int OLDFORMAT;
   char MATFILE[4096];
   char TREEFILE[4096];
   char ORDFILE[4096];
   char TRANSFILE[4096];
   char STAMPDIR[4096];
   };
/* Standard structure for storing protein sequence data */

struct seqdat {	/* all lengths include char terminator and [0] */
    int ilen;	/* length of identifier*/
    char *id;	/* identifier */
    int tlen;	/* length of title */
    char *title;/* title */
    int slen;	/* length of sequence*/
    char *seq;	/* sequence */
};
/* tree structures */

#if !defined(CLUST_STRUCT) 
struct indclust {
	int number;
	int *member;
	};
struct cluster {
	struct indclust a;
	struct indclust b;
	};
#define CLUST_STRUCT
#endif
#include <gjutil.h>

struct cluster *readtree(char *tordfile, char *treefile, int *number, int method, FILE *OUT);
float fmatfit(float **atoms1, float **atoms2, float **R, float *V, int nats, int entry);
struct cluster *get_clust(double **matrix, char **ids, int ndomain, char *noc_parms);
void closefile(FILE *handle,char *filename);
FILE *openfile(char *filename, char *type);
char ltou(char c);
char utol(char c);
void fmatmult(float **r, float *v, float **coord, int n);
void update(float **dR, float **R, float *dV, float *V);
int count_domain(FILE *IN);
void rmsp(char *c);
int getca(FILE *IN, float **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type,  int MAXats, int REVERSE, FILE *OUTPUT);
int newoutput(FILE *TRANS, struct domain_loc *domain, int ndomain, int writetrans);
void exit_error();
int getpars(FILE *fp, struct parameters *var);
int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);
int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, 
        int *gottrans, char *env, int DSSP, FILE *OUTPUT);

int getdomain_error(char *buff); 
