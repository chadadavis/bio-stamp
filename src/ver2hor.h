#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_STAMP_NUM 10
#define MAX_N_SEQ 1000
struct parameters{
   char filename[200];
   char TYPE;
   float CUTOFF;
   int WINDOW;
   int SECDISP;
   int SMALLSEQ,SMALLSEC;
   int SECSUM;
   int CASESEQ,CASESEC;
   int VERBOSE;
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
int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos);
int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window);
int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);
char utol(char c);
char ltou(char c);
void exit_error();
