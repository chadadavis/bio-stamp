#define MAX_STAMP_NUM 10
#define MAX_N_SEQ 1000
struct parameters{
   char filename[200];
   char TYPE;
   float CUTOFF;
   int WINDOW;
   float POINTSIZE;
   int SEC;
   int BOXSEQ,BOXSEC;
   int SMALLSEQ,SMALLSEC;
   int BOLDSEQ,BOLDSEC;
   int SECSUM;
   int CASESEQ,CASESEC;
   int VERBOSE;
   char prefix[200];
};
struct stampdat{
   char what;
   char *title;
   float *n;
};
struct seqdat { /* all lengths include char terminator and [0] */
    int ilen;   /* length of identifier*/
    char *id;   /* identifier */
    int tlen;   /* length of title */
    char *title;        /* title */
    int slen;   /* length of sequence*/
    char *seq;  /* sequence */
};
int getstampdat();
int *getstamprel();
int *getstampsecsum();
