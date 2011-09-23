
#define MAXslen 2000
#define MAXnbloc 200
#define MAXtlen 200

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
   };

struct domain_loc {		/* This structure allows rather complex domains to be described */
   char filename[100];
   char id[30];
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
   char bloc_file[200];
   char dom_file[200];
   char parm_file[200];
   char prefix[200];
   int MAX_SEQ_LEN;
   int PAIRWISE;
   int TREEWISE;
   int OLDFORMAT;
   char MATFILE[100];
   char TREEFILE[100];
   char ORDFILE[100];
   char TRANSFILE[100];
   char STAMPDIR[200];
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
#include "gjutil.h"
struct cluster *readtree();
float matfit();
struct cluster *get_clust();
