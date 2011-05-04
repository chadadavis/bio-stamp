#include <malloc.h>

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

char **RBR_c_split(char *str, int *n,  char delimiter);
struct seqdat *RBR_getmsfalign(FILE *IN, int *nseq);
