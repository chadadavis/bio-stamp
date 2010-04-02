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

/* smseq Structure to store the details of an alignment - the alignment
itself is stored in an array of seqdat structures
*/

struct smseq {
	int nseq;		/* number of seqs in the alignment */
	int blen;		/* length of the alignment = no of residues +2 */
	char *title;		/* a name for the alignment */
	int ninfo;              /* number of lines of info */
	char **info;		/* optional information about the alignment */
	char *gapchars;		/* string of valid gap-characters */
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

/* structure to hold default values for id lengths etc. */
struct defstr {
	int maxilen;	/* max id len */
	int maxtlen;	/* max title len */
	int maxslen;	/* max sequence length */
};





