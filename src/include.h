#include <malloc.h>
#define defaultsfile "$STAMPDIR/stamp.defaults"
#define max(A,B) ((A) > (B) ? (A) : (B))
#define leq(A,B) ((A) <= (B) ? (A) : (B))
#define max4(A,B,C,D) ((A)>(B))?(((A)>(C))?(((A)>(D))?(A):(D)):(((C)>(D))?(C):(D))):(((B)>(C))?(((B)>(D))?(B):(D)):(((C)>(D))?(C):(D)))

struct brookn { /* structure to represent brookhaven residue numbers */
    int n;      /* numerical part of number */
    char cid;   /* chain identifier */
    char in;    /* insertion code */
};

struct domain_loc{		/* This structure allows rather complex domains to be described */
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
   int *reverse;			/* if 1, then reverse invert the object N to C */
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

struct parameters {
	float E1,E2;			/* Rossmann and Argos Parameters */
	float first_E1,first_E2;	/* if more than one pass is suggested, use these first */
	float second_E1,second_E2;	/* if more than one pass is suggested, use these second */
	float const1,const2;		/* 2*E1*E1, 2*E2*E2 */
	float CUTOFF;			/* Cutoff value for fitting */
	float first_CUTOFF;		/* when NPASS == 2 */
	float second_CUTOFF;		/* when NPASS == 2 */
	float ADD;			/* Value added to matrix */
	float NMEAN;			/* Statistics values from the paper version */
	float NSD,NA,NB;
	float NASD,NBSD;
	float SCORETOL;			/* Score tolerance for convergence */
	float THRESH;			/* Threshold for fitting */
        float CPUtime;			/* running total of CPU time used by the program */
	int PRECISION;			/* Precision for float to int conversion */
	int MAX_SEQ_LEN;		/* Maximum number length of protein that can be used */
	float PAIRPEN; 
	float first_PAIRPEN;
	float second_PAIRPEN; 		/* Gap penalties */
	float first_THRESH;		/* if the Sc value is not >= this number then do not do the second fit (NPASS==2 only) */
	float TREEPEN;
	float first_TREEPEN;
	float second_TREEPEN;		
	int DSSP;			/* if true, read coordinates from DSSP files, not PDB files */
	int PAIROUTPUT;			/* if true, output all pairwise alignments */
	int MAXPITER,MAXTITER;		/* Maximum number of iterations */		
	int CLUSTMETHOD;		/* What to derive tree from 1=rms, 2=Sc */
	int MAXLEN;			/* Maximum sequence length */
	int MAXATS;			/* Maximum number of atoms */
	int PAIRPLOT,TREEPLOT;		/* 1 = plot distance matrix */
	int PAIRALIGN,TREEALIGN;	/* 1 = display pairwise or treewise alignments */
	int PAIRALLALIGN, TREEALLALIGN; /* 1 = display pairwise or treewise alignments during each iteration */
        int NALIGN;			/* Number of alignments displayed */
	int DISPALL; 			/* 1 = display all alignments */
	int HORIZ;			/* 1 = display horizontal alignments */
	int STATS;			/* 1 = use statistics rather than supplied parameters */
	int PAIRWISE; 			/* 1 = do pairwise calculations */
	int TREEWISE;			/* 1 = do treewise calculations */
	int SCAN;			/* 1 = scan the named domain database */
	int SCANMODE;			/* 0 = report score only; 1 = report alignment and transformation */
	int SCANALIGN;			/* 1 = display scan alignment */
	float SCANCUT;			/* only report (1 above) for alignments > SCANCUT */
	int SLOWSCAN;			/* new method of obtaining initial superimpositions */
	int SCANSLIDE;			/* align the N-terminus of the query with every SCANSLIDEth amino acid on the database structure */
	int SECSCREEN;			/* When scanning, if SECSCREEN is set to one, then an initial comparison of 
					 *  3 state secondary structure content will be performed (if secondary
					 *  structure information is available.  If the total difference in content
					 *  (ie. diff helix + diff sheet + diff other) is greater than SECSCREENMAX,
					 *  then the comparison will be skipped */
	float SECSCREENMAX;		/* See above, the value is in percent */
	int SCANTRUNC;			/* 1 = do truncation of structures to be scanned (fraction of size of query) */
	float SCANTRUNCFACTOR;		/* Factor for the above */
	int NPASS;			/* 1 = just do one comparison with E1=first_E1,E2=first_E2
					 * 2 = do a comparison with E1=first_E1, E2=first_E2 then one with
					 *      E1=second_E1, E2=second_E2 (more sensitive but slower) */
	int BOOLEAN;			/* 1 = use new boolean method */
	float BOOLCUT;			/* cutoff for Pij if boolean method is used */
	int BOOLMETHOD;			/* 0 = use a 1 or a 0 in the matrix if all pairwise comparisons are favorable (pen = 0)
					 * 1 = each position will be the total number of equivalent pairwise comparisons (pen =
					 *     PAIRPEN * Nx(N-1)/2) */
	int BOOLFRAC;			/* BOOLMETHOD=1 only; minimum fraction of Nx(N-1)/2 of comparisons required to have
					 *  Pij>BOOLCUT allowed for a position to be considered equivalent */
	float first_BOOLCUT;
	float second_BOOLCUT;
	int SCANSEC;
	int CO;				/* 1 = Cut output (i.e. cut domain descriptors in output)
					 * 0 = Leave descriptors as is */
	int SEC;			/* 0 = do not use secondary structure
					 * 1 = use Kabsch and Sander's DSSP for secondary structures 
					 * 2 = read in assignments from a supplied summary file */
	int ROUGHFIT;			/* 1 = calculate a rough initial transformation by simply aligning
					 * the sequences from their N-terminal ends */
	int CLUST;			/* 1 = use single linkage cluster analysis (maketree) to generate
					 *  a tree and order file, else the program expects a tree 
					 *  file to be supplied (ie generated by another method, such as UPGMA */

	int COLUMNS;			/* number of columns for alignment purposes */
	int SW;				/* SW = 0 ==> normal Smith Waterman Routine; 1 ==> cut corners */
	float CCFACTOR;			/* Cutting corners factor */
	int CCADD;			/* 1 ==> add sequence length difference to CCFACTOR; 0 ==> do not */
	int MINFIT;			/* minimum number of fitted atoms (in residues -- the program calculates a clever minscore to use) */
	float MIN_FRAC;			/* For scanning, the minimum length (fraction) of database sequence to be compared (things shorter are to be ignored) */
	int SCORERISE;			/* If set TRUE, then a drop in score will result in the end of the comparison */
	int SKIPAHEAD;			/* For scans. If set TRUE, then skip ahead if a transformation has been output (ie. avoid using the same part of the structure again) */
	int SCANSCORE;			/* See doc/stamp.doc */
	int ALLPAIRS;			/* All pairs of comparisons to be performed */
	char listfile[100];		/* List of domains to be used */
	char secfile[100];		/* used when SEC==3 to supply secondary structure assignments */
	char ordfile[100];		/* Tree order file */
	char treefile[100];		/* Tree file */
	char plotfile[100];		/* Tree plot file */
	char transprefix[100];		/* Transformation file prefix */
	char matfile[100];		/* Matrix file */
	char scanfile[100];		/* Scan file */
	char database[100];		/* Data base for scanning */
	char logfile[100];		/* Stream output file  -- replaces standard output */
	char dsspfile[100];		/* file containing a list of DSSP directories */
	char roughalign[100];		/* ROUGHFIT aligment file */
	char stampdir[100];		/* $STAMPDIR environment variable */
	FILE *LOG;			/* The opened version of the above (passed everywhere) */
	int verbose;
	};
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
struct cluster *readtree();
float rossmann();
float pairpath();
int swstruc();
float matfit();
struct path *dosort();
void update();
float idist();
char *getfile();
char *getrksum();
struct cluster *get_clust();
