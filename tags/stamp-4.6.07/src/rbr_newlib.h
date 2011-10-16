#define RBR_MAX_ID_LEN 50
#define RBR_MAX_TITLE_LEN 200

char *RBR_a1to3();
char RBR_a3to1();
int  RBR_readcbin();
char *RBR_getfile();
float RBR_idist();
float RBR_dist();
float **RBR_readpm();
struct brookn *RBR_sum_get_brk();
struct range RBR_sum_get_brange();
int *RBR_sum_get_index();
int RBR_compseq();
struct seqdat RBR_readseq();
float RBR_r_of_gyration();
float RBR_r_of_gyration_pep();
int *RBR_c_of_m();
int *RBR_c_of_m_pep();
float RBR_ipowdist();
char *RBR_skiptononspace();
int **RBR_readmdm();
int **RBR_readsolv();
float RBR_imatfit(int **atoms1, int **atoms2, float **R, float *V, int nats, int entry, int PRECISION);
float RBR_fmatfit();
float RBR_dayhoff_score();
float RBR_sequence_identity();
float RBR_ran3();
int RBR_getdsspsum();
FILE *RBR_zopen(char *filename, char *type);
struct seqdat *RBR_getseq(FILE *IN, int *nseq);
struct seqdat *RBR_getseqfasta(FILE *IN, int *nseq);
char **RBR_c_split(char *str, int *n,  char delimiter);
int *RBR_seq_index(char *seq1, char *seq2);

float RBR_iwmatfit(int **atoms1, int **atoms2, float *w, float **R, float *V,  int nats, int entry, int PRECISION);
