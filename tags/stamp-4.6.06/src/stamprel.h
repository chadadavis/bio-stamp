struct stampdat{
   char what;
   char *title;
   float *n;
};

int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos);

int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window);

int stamp_clean_block(struct seqdat *bloc, int nbloc, int window, float min_frac,struct stampdat *stamp, int nstamp);


