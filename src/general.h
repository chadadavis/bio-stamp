/* #include <malloc.h> */
#include <stdlib.h>

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
};

struct indclust {
	int number;
	int *member;
	};
struct cluster {
	struct indclust a;
	struct indclust b;
	};
/* Matrix file structure - eg. Dayhoffs matrix */
struct pmatrix {
    char *title;    /* title of matrix */
    char *indx;	    /* Index to matrix - ie amino acid codes */
    int  inlen;	    /* number of residues per line */
    int **array;    /* the actual pair score array */
};
/* Standard structure for storing protein sequence data */
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

struct sident {
    char *id;	    /* sequence id */
    char *title;    /* sequence title */
    int uid;	    /* unique identifier number for sequence */
};

struct alseq {
    char *line;	    /* one position in the sequence alignment */
};

struct mseq {
    int nseq;		/* number of sequences in this aligned bloc */
    int blen;		/* overall aligned length of the bloc */
    struct sident *itd;	/* list of id,title,uid structures */
    struct alseq *bloc;	/* pointer to list of blen lines of the alignment
			   structure used to allow future expansion of info
			   on a per-aligned line basis (eg. presence of 
			   flexible gap, etc)*/
};

struct coord{
	int i;		/* to store i and j of seqa[i],seqb[j] */
	int j;
};

struct path{
	struct coord start;	/* path start point i and j*/
	struct coord end;	/* path end point i and j*/
	int score;		/* max score on this path */
	int col;		/* current score on this path */
};

struct olist{
	int len;
	struct path *res;
};
struct side_chain {
        int n;
	int **names;
	int **coords;
};
	
