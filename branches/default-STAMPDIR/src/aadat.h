#define RBR_AA1  "ABCDEFGHIJKLMNOPQRST_VWXYZ________c____h"
#define RBR_AA1a "A_CDEFGHI_KLMN_PQRST_VW_Y_______________"
#define RBR_AA1b "ARNDCQEGHILKMFPSTWYVBZX*"
#define RBR_AA3 "ALA ASX CYS ASP GLU PHE GLY HIS ILE CSH LYS LEU MET ASN CSS PRO GLN ARG SER THR ___ VAL TRP UNK TYR GLX ___ ___ ___ ___ ___ ___ ___ ___ CYH"
#define RBR_SS "-STGHIBE"
#define RBR_AA20 "A_CDEFGHI_KLMN_PQRST_VW_Y_"
#define RBR_AA20A  "A_CDEFGHI_KLMN_PQRST_VW_Y "
#define RBR_SSPRED "HBc_"
/* General notes:
 *
 * When wishing to refer to C_ss and C_sh seperately, I define C_sh as one letter 
 *  code "J" and  C_ss as one letter code "O".
 *
 * Given a one-letter code, to get the data for that amino acid one need merely
 *  refer to <data structure>[(int)(<code> - 'A')].  For example, to get the
 *  accessibility of Phenylalanine:
 *      AA_ACC[(int)('F'-'A')]
 *
 */

/* Values of standard state mean surface area.  Taken from:
 *  Rose, GD and Dworkin, JE `Hydrophobicity Profile', In:
 *  *Prediction of Protein Structure and the Principles of Protein
 *   Conformation*, Fasman, GD editor, Plenum press (1989). Page 629.
 * Value for Asx is the average of Asn/Asp; value for Glx is the 
 *  average for Gln/Glu; value for Unk is a complete guess. */
static int RBR_AA_ACC[40] =
 /* Ala, Asx, Cys, Asp, Glu, Phe, Gly, His, Ile, Cyh, Lys, Leu, Met, Asn, Cys */
 { 118, 162, 146, 158, 186, 222,  88, 203, 181, 146, 226, 193, 203, 166,  146,
 /* Pro, Gln, Arg, Ser, Thr, ___, Val, Trp, Unk, Tyr, Glx,---------------,Cys */ 
   147, 193, 256, 130, 153,   1, 165, 266, 400, 237, 190,1,1,1,1,1,1,1,1,146,1,1,1,1,400 };

/* Taken from Rost & Sander, Proteins, 20, 216-226, 1994 */
static int RBR_ROST_ACC[40] =
 /* Ala, Asx, Cys, Asp, Glu, Phe, Gly, His, Ile, Cyh, Lys, Leu, Met, Asn, Cys */
 { 106, 160, 135, 163, 194, 197,  84, 184, 169, 135, 205, 164, 188, 157,  135,
 /* Pro, Gln, Arg, Ser, Thr, ___, Val, Trp, Unk, Tyr, Glx,---------------,Cys */
   136, 198, 248, 130, 142,   1, 142, 227, 400, 222, 196,1,1,1,1,1,1,1,1,135,1,1,1,1,400 };


 /* Amino acid abundance as calculated on REP.SEQ during the environment study */
 static float RBR_AA_ABUN[40] =
 /* Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
 {  8.5,  0.1,  1.0,  5.8,  5.7,  3.6,  8.5,  2.6,  5.3,  1.0,  6.0,  8.1,  2.0,  
 /* Asn,  Css,  Pro, Gln, Arg, Ser, Thr, ___, Val, Trp, Unk, Tyr, Glx */
    4.7,  1.0,  4.7, 3.5, 4.2, 6.9, 6.2, 0.1, 7.7, 1.3, 0.0, 3.6, 0.1,
 /* ------,                         cys-sh */ 
    0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1.0,0.0,0.0,0.0,0.0,0.0 };

 /* Molecular weights from Rawn (1984) 
  * Note that these are the weights for isolated amino acids */
 static int RBR_AA_MASS_ALL[40] = 
 /* Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
 {   89,  133,  121,  133,  147,  165,   75,  155,  131,  121,  146,  131,  149,
/*  Asn,  Cys,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */ 
    132,  121,  115,  146,  174,  105,  119,    0,  117,  204,  110,  181,  147,
/*  ...            Cys-sh */
    0,0,0,0,0,0,0,0, 121,0,0,0,0,110 };

/* The above, less 18.0 for each one (i.e. for calculating the weight of a polypeptide */
static int RBR_AA_MASS[40] = 
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
 {   71,  115,  103,  115,  129,  147,   57,  137,  113,  103,  128,  113,  131,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
    114,  103,   97,  128,  156,   87,  101,    0,   99,  186,   92,  163,  129,
/*  ...            Cys-sh */
    0,0,0,0,0,0,0,0, 103,0,0,0,0,110 };

/* number of heavy atoms in the amino acid (including main chain) */
static int RBR_AA_NATS[40] =
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
{     5,    8,    6,    8,    9,   11,    4,   10,    8,   10,    9,    8,    8,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
      8,    6,    7,    9,   11,    6,    7,    0,    7,   14,   99,   12,    9,
/*  ...            Cys-sh */
    0,0,0,0,0,0,0,0,6,0,0,0,0,99 };

/* The next two scales were taken from Branden and Tooze, page 210 */
/* Kyte and Doolittle hydrophobicity scale -- multiplied by 10 */
static int RBR_HYDR_KD[40] =
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
{    18,  -35,   25,  -35,  -35,   28,   -4,  -32,   45,   25,  -39,   38,   19,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
    -35,   25,  -16,  -35,  -45,   -8,  -13,    0,   42,   -9,    0,  -13,  -35,
/*  ...            Cys-sh */
      0,0,0,0,0,0,0,0,25 ,0,0,0,0,0};

/* Engleman, Steits and Goldman hydrophobicity scale -- multiplied by 10 */
static int RBR_HYDR_ESG[40] =
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
{    16,  -70,   20,  -92,  -82,   37,   10,  -30,   31,   20,  -88,   28,   34,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
    -48,   20,   -2,  -41, -123,    6,   -7,    0,   26,   19,    0,   -7,  -62,
/*  ...            Cys-sh */
     0,0,0,0,0,0,0,0,20,0,0,0,0,0 };

/* At pH 7 */
static int RBR_AA_CHARGE[40] =
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
{     0,    0,    0,   -1,   -1,    0,    0,    0,    0,    0,   +1,    0,    0,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
      0,    0,    0,    0,   +1,    0,    0,    0,    0,    0,    0,    0,    0,
/*  ...            Cys-sh */
     0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

/* A value of 100 implies an unknown or very high value */
static float RBR_AA_PK[40] =
/*  Ala,  Asx,  Cys,  Asp,  Glu,  Phe,  Gly,  His,  Ile,  Csh,  Lys,  Leu,  Met */
{   100,  100,  100,  3.9,  4.3,  100,  100,  6.0,  100,  100, 10.5,  100,  100,
/*  Asn,  Css,  Pro,  Gln,  Arg,  Ser,  Thr,  ___,  Val,  Trp,  Unk,  Tyr,  Glx */
    100,  100,  100,  100, 12.5,  100,  100,  100,  100,  100,  100, 10.1,  100,
/*  ...            Cys-sh */
     0,0,0,0,0,0,0,0,100,0,0,0,0,0 };
