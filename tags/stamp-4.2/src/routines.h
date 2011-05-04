char a3to1(char *a3);

int Agetbloc(FILE *bfile, struct seqdat *bloc, int *nbloc);

int getpars(FILE *fp, struct parameters *var);

int aliseq(char *seqa, char *seqb, struct path *apath, char **patha, 
           char *alseq, char *blseq, int *k, int *hgap, int *vgap);

int ccprobcalc(int **atoms1, int **atoms2, int **prob, int lena, int lenb, 
	struct parameters *parms);

int clean_block(struct seqdat *bloc, int nbloc, int window);

int count_domain(FILE *IN);

int disp(struct domain_loc domain, FILE *OUTPUT);

int display_align(char **seqa, int na, char **seqb, int nb, 
            char **seca, char **secb, char *fit, char *value, 
            int cols, int printsec, int printdat, FILE *OUTPUT);

struct path *dosort(struct olist *result, int *lena, int *total);

int readcom(int n, char **args, struct parameters *parms);

int extract_dssp(FILE *IN, struct brookn start, struct brookn end,int type,
        float **R, float *V, int startats, char chainlabel, FILE *OUT);

int extract_pdb(FILE *IN, struct brookn start, struct broon end, int type,
        float **R, float *V, int startats, int HETERO, int NUCLEIC, int HOH,
        char chainlabel, int verbose, char *filename, FILE *OUT);

float fmatfit(float **atoms1, float **atoms2, float **R, float *V, int nats, int entry);

void fmatmult(float **r, float *v, float **coord, int n);

struct cluster *get_clust(double **matrix, char **ids, int ndomain, char *noc_parms);

int getca(FILE *IN, float **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type,  int MAXats, int REVERSE, FILE *OUTPUT);
        
int getdomain(FILE *IN, struct domain_loc *domains, int *ndomain, int maxdomain, 
        int *gottrans, char *env, int DSSP, FILE *OUTPUT);

int get_dssp_sum(FILE *DSSP, struct brookn begin, struct brookn end, int type,
        char *sec, int REVERSE, int maxlen, int *count, FILE *OUT);

char *getfile(char *code_safe, char *dirfile, int code_length, FILE *OUTPUT);

int getks(struct domain_loc *domain, int ndomain, struct parameters *parms);

char *getrksum(FILE *IN);

int getsec(struct domain_loc *domain, int ndomain, struct parameters *parms);

int getstampdat(struct stampdat *stamp, FILE *IN, int *nstamp, int *nseq, int *npos, int maxpos);

int *getstamprel(struct stampdat *stamp, int nval, int npos, char type, float cutoff, int window);

float idist(int *atm1, int *atm2, int PRECISION);

int igetca(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct broonk end, int type, int MAXats,
        int REVERSE; int PRECISION; FILE *OUTPUT);

int igetcadssp(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type, int MAXats, 
        int REVERSE, int PRECISION, FILE *OUTPUT);

int igetcb(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type, int MAXats,
        int REVERSE, int PRECISION, FILE *OUTPUT);

int igetgen(FILE *IN, int **coords, char *aa, struct brookn *numb, int *ncoord,
        struct brookn start, struct brookn end, int type, char *atom,
        int MAXats, int REVERSE, int PRECISION, FILE *OUTPUT);

int igetside(FILE *IN, struct side_chain *side, char *aa, struct brookn *numb, int *nres,
        struct brookn start, struct brookn end, int type, int MAXats,
        int REVERSE, int PRECISION, FILE *OUTPUT);

char ltou(char c);

int makefile(struct domain_loc *domain, int ndomain, struct cluster cl, 
        int nclust, float score, floatrms, int length, int nfit,
        float *Pij, float *Dij, float *dist, float *Pijp,
        int PAIRWISE, struct parameters *parms);

float matfit(int **atoms1, int **atoms2, float **R, float *V,
        int nats, int entry, int PRECISION);


void matinv(float **a, float **y, float d, int *indx);

int matmult(float **r, float *v, int **coord, int n, int PRECISION);

int matprod(float **P, float **A, float **B, FILE *OUTPUT);

int matvecprod(float **A, float *C, float *B, FILE *OUTPUT);

int newoutput(FILE *TRANS, struct domain_loc *domain, int ndomain, int writetrans);

int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans);

int pairfit(struct domain_loc *domain1, struct domain_loc *domain2, float *score, float *rms,
        int *length, int *nfit, struct parameters *parms, int rev,
        int *start1, int *end1, int *start2, int *end2,
        float *seqid, float *secid, int *nequiv, int *nsec, int **hbcmat,
        int ALIGN, int count, int FINAL);

float pairpath(struct domain_loc domain1, struct domain_loc domain2, int **prob,
        long int entry, float **R2, float *V2, int *len, float *score,
        int *nfit, char *ftouse, float *fpuse, int *start1, int *end1,
        int *start2,int *end2, struct parameters *parms);

int pairwise(struct domain_loc *domain, int ndomain, struct parameters *parms);

int printdomain(FILE *TRANS, struct domain_loc domain, int writetrans);

int printmat(float **R, float *V, int n, FILE *OUTPUT);

int probcalc(int **atoms1, int **atoms2, int **prob, int lena, int lenb,
        struct parameters *parms);

void probplot(int **prob, int lena, int lenb, int n, int cutoff, FILE *OUTPUT);

int qkfit(doublereal *umat, doublereal *rtsum, doublereal *r, integer *entry_);

struct cluster *readtree(char *tordfile, char *treefile, int *number,
        int method, FILE *OUT);

int reval(char *seq, int start, int end);

int revmatmult(float **r, float *v, int **coord, int n, int PRECISION);

void rmsp(char *c);

float rossmann(int **atoms1, int **atoms2, int start, int end,
        float const1, float const2, float *Dij, float *Cij, int PRECISION);

int roughfit(struct domain_loc *domain, int ndomain, struct parameters *parms);

int scan(struct domain_loc domain, struct parameters *parms);

int sec_content(char *sec, int npos, int type, int *pera, int *perb, int *perc);

float seq_identity(char *seq1, char *seq2, int *nid, int *align_len, FILE *OUTPUT);

int slow_scan(struct domain_loc qdomain, struct parameters *parms);

int smoothsec(char *sec, int minhelixlen, int minstrandlen);

int stamp_clean_block(struct seqdat *bloc, int nbloc, int window,
        struct stampdat *stamp, int nstamp);

int getpars(FILE *fp, struct parameters *var);

int sw7ccs(int  lena, int lenb, int **prob, int pen, struct olist *result,
        int *total, unsigned char **patha, int min_score, 
        int auto_corner, float fk);

int swstruc(int  lena, int lenb, int pen, int **prob, struct olist *result,
        int *total, unsigned char **patha, int min_score);

int testfile(char *file);

int threestate(char *sec,char *helix,char *extended,char *coil);

int treefit(struct domain_loc *domain, int ndomain, struct cluster cl, 
        float *score, float *rms, int *length, int *nfit,
        float *Pij, float *Dij, float *dist, float *Pijp,
        int rev, int align, struct parameters *parms);

int treepath(struct domain_loc *domain, int ndomain, struct cluster cl,
        float **R2, float *V2, int **prob, float *score, float *rms,
        int *length, int *nfit, float *Pij, float *Dij, float *distance,
        float *Pijp, float mean, float sd, char *fpuse, char *ftouse,
        struct parameters *parms);

int treewise(struct domain_loc *domain, long int ndomain,       
        struct parameters *parms);

void update(float **dR, float **R, float *dV, float *V);

char utol(char c);

int RBR_print_vector(float *V);

int RBR_vector_unify(float *V);

float *RBR_vector_diff(float *V1, float *V2);

float *RBR_vector_ave(float *V1, float *V2);

float *RBR_vector_cross(float *V1, float *V2);

float *RBR_vector_set_dist(float *V1, float R);
