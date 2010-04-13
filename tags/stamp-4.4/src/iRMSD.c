#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#define TINY 1.0e-20

// A faster iRMSD calculation
// Takes as input:
//   1. transformations (in one file)
//   2. cofm data
//   3. interactions (dom identifier pairs)
// Then outputs iRMSDs


struct trans {
   float *V;
   float **R;
};

struct trans_set {
  struct trans *t;
  char **id;
  int *nid;
  char **f;
  char **des;
  int n;
};

struct cofm {
  char *id;
  int nid;
  int N;
  float Rg;
  float Rmax;
  float Ro[3];
  float coord[10][3];
};

struct trans_index {
   int n;
   int *t;
   int *p;
};


char **c_split(char *str, int *n,  char delimiter) {
        /* breaks a string up into substrings using a delimiter */

        int i,j,k,len;
        char **sstr;

        sstr=(char**)malloc(sizeof(char*));
        
        (*n)=0;
        i=0;
        len=strlen(str);

        while(i<len && str[i]!='\0' && str[i]!='\n') {
           while(i<len && str[i]==delimiter) i++;
           if(i>=len || str[i]=='\0' || str[i]=='\n') break;
           /* we are at a new sub-string */
           sstr[(*n)]=(char*)malloc(sizeof(char));
           j=0;
           while(i<len && str[i]!='\0' && str[i]!='\n' && str[i]!=delimiter) {
                sstr[(*n)][j]=str[i];
                j++; i++;
                sstr[(*n)]=(char*)realloc(sstr[(*n)],(j+1)*sizeof(char));
           }
           sstr[(*n)][j]='\0';
           (*n)++;
           sstr=(char**)realloc(sstr,((*n)+1)*sizeof(char*));
           i++;
        }
//        printf("String %s breaks into %d sub-strings delimited by '%c' characters\n",
//                str,(*n),delimiter);
//        for(i=0; i<(*n); ++i) printf("%s, ",sstr[i]);
//        printf("\n"); 

        return sstr;
}

struct trans_set *get_trans(char *transfile, int *nt) {
    int i, j, k;
    int in_trans;
    int line;
    int ntokens;
    int nd;

    char buff[1000];
    char **tokens;

    FILE *IN;

    if((IN = fopen(transfile,"r")) == NULL) {
	fprintf(stderr, "Error reading file %s\n", transfile);
	exit(-1);
    }

    struct trans_set *TR;

    in_trans = 0;
    line = 0;

    (*nt) = -1;
    nd = -1;

    TR = (struct trans_set *) malloc(sizeof(struct trans_set));
 
    while((fgets(buff, 999, IN)) != NULL) {
	if(strstr(buff,"TRANS_BEGIN") != NULL) {
	    in_trans = 1;
	    line = 0;
	    nd = -1;
	    (*nt)++;

	    /*
	    if((((*nt) + 1) % 10) == 0) {
		printf("Total number of trans entries: %10d", (*nt) + 1);
	    }
	    if((((*nt) % 10) == 0) && ((*nt) > 0)) {
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    }
	    */

	    TR = (struct trans_set *) realloc(TR, ((*nt) + 2) * sizeof(struct trans_set));
	    TR[(*nt)].n = 0;
	    TR[(*nt)].f = (char **) malloc(sizeof(char *));
	    TR[(*nt)].id = (char **) malloc(sizeof(char *));
	    TR[(*nt)].nid = (int *) malloc(sizeof(int));
	    TR[(*nt)].des = (char **) malloc(sizeof(char *));
	    TR[(*nt)].t = (struct trans *) malloc(sizeof(struct trans));
	    //printf("So far %d trans buff %s",(*nt),buff);
	}
	else if(strstr(buff, "TRANS_END") != NULL) {
	    in_trans = 0;
	}
	else if((in_trans == 1) && (strchr(buff, '{') != NULL)) {
	    line = 1;
	}
	if((in_trans == 1) && (line > 0)) {
	    //printf("TRANS %d LINE %d %s",(*nt),line,buff);

	    tokens = c_split(buff, &ntokens, ' ');

	    if(line == 1) { // Title line
		TR[(*nt)].n++;
		nd++;
		TR[(*nt)].f[nd] = (char *) malloc(100 * sizeof(char));
		TR[(*nt)].f = (char **) realloc(TR[(*nt)].f, ((nd + 2) * sizeof(char *)));
		strcpy(TR[(*nt)].f[nd], tokens[0]);
		TR[(*nt)].id[nd] = (char *) malloc(100 * sizeof(char));
		TR[(*nt)].id = (char **) realloc(TR[(*nt)].id, ((nd + 2) * sizeof(char *)));
		strcpy(TR[(*nt)].id[nd], tokens[1]);
		TR[(*nt)].nid = (int* ) realloc(TR[(*nt)].nid,((nd + 2) * sizeof(int)));

		TR[(*nt)].des[nd] = (char *) malloc(100 * sizeof(char));
		TR[(*nt)].des = (char **) realloc(TR[(*nt)].des, ((nd + 2) * sizeof(char *)));
		i = 0;
		while(buff[i] != '{') { i++; }
		while(buff[i] != ' ') { i++; }
		for(j = i; j < strlen(buff); ++j) {
		    if(buff[j] == '\n') { break; }
		    TR[(*nt)].des[nd][j - i] = buff[j];
		}
		TR[(*nt)].des[nd][j - i]='\0';

		TR[(*nt)].t = (struct trans *) realloc(TR[(*nt)].t, ((nd + 2) * sizeof(struct trans)));
		TR[(*nt)].t[nd].R = (float **) malloc(3 * sizeof(float *));
		for(i = 0; i < 3; ++i) {
		    TR[(*nt)].t[nd].R[i] = (float *) malloc(3 * sizeof(float));
		}
		TR[(*nt)].t[nd].V = (float *) malloc(3 * sizeof(float));

		//printf("Assigning domain %d in trans %d |%s|%s|%s|\n",nd,(*nt),
		//TR[(*nt)].f[nd],TR[(*nt)].id[nd],TR[(*nt)].des[nd]);
		//printf("Buff is %s",buff);

	    }
	    else if(line == 2) { // 1st row
		sscanf(
		       buff,
		       "%f %f %f %f",
		       &TR[(*nt)].t[nd].R[0][0],
		       &TR[(*nt)].t[nd].R[0][1],
		       &TR[(*nt)].t[nd].R[0][2],
		       &TR[(*nt)].t[nd].V[0]
		       );
	    }
	    else if(line == 3) { // 2nd row
		sscanf(
		       buff,
		       "%f %f %f %f",
		       &TR[(*nt)].t[nd].R[1][0],
		       &TR[(*nt)].t[nd].R[1][1],
		       &TR[(*nt)].t[nd].R[1][2],
		       &TR[(*nt)].t[nd].V[1]
		       );
	    }
	    else if(line == 4) { // 3rd row
		sscanf(
		       buff,
		       "%f %f %f %f",
		       &TR[(*nt)].t[nd].R[2][0],
		       &TR[(*nt)].t[nd].R[2][1],
		       &TR[(*nt)].t[nd].R[2][2],
		       &TR[(*nt)].t[nd].V[2]
		       );
	    }

	    for(i = 0; i < ntokens; ++i) { free(tokens[i]); }
	    free(tokens);
	}
	line++;
    }
    fclose(IN);
    (*nt)++;
    return(TR);
}

struct cofm *get_cofm(char *cofmfile, int *np) {
    int i, j, k;
    int ntokens;
    int natom;
    char buff[1000];
    char temp[10];
    char **tokens;
    struct cofm *C;
    FILE *IN;

    (*np) = -1;

    if((IN = fopen(cofmfile,"r")) == NULL) {
	fprintf(stderr,"Error reading file %s\n",cofmfile);
	exit(-1);
    }

    C = (struct cofm *) malloc(sizeof(struct cofm));
  
    natom = 0;
    while((fgets(buff, 999, IN)) != NULL) {
	if(strncmp(buff, "REMARK", 6) == 0) {
	    tokens = c_split(buff, &ntokens, ' ');

	    //REMARK Domain   1 Id 1uvyA.a.1.1.1-1 N = 116 Rg =  13.074 Rmax =  18.921 Ro   =   49.879   25.058   31.369
	    (*np)++; 

	    /*
	    if((((*np)+1)%100)==0) {
		printf("Total number of cofm entries:  %10d",(*np)+1);
	    }
	    if((((*np)%100)==0) && ((*np)>0)) {
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    }
	    printf("%s",buff);
	    */

	    C[(*np)].id = (char * ) malloc(100 * sizeof(char));
	    C = (struct cofm *) realloc(C, (((*np) + 2) * sizeof(struct cofm)));
	    strcpy(C[(*np)].id, tokens[4]);
	    sscanf(tokens[7], "%d", &C[(*np)].N);
	    sscanf(tokens[10], "%f", &C[(*np)].Rg);
	    sscanf(tokens[13], "%f", &C[(*np)].Rmax);
	    sscanf(tokens[16], "%f", &C[(*np)].Ro[0]);
	    sscanf(tokens[17], "%f", &C[(*np)].Ro[1]);
	    sscanf(tokens[18], "%f", &C[(*np)].Ro[2]);
	    natom = 0;

	    for(i = 0; i < ntokens; ++i) {
		free(tokens[i]);
	    }
	    free(tokens);
	}
	else if(strncmp(buff, "ATOM  ", 6) == 0) {
	    // ATOM lines are not space-delimited so can't use c_split
	    
	    memset(temp, '\0', 10);
	    strncpy(temp, buff + 30, 8);
	    sscanf(temp, "%f", &C[(*np)].coord[natom][0]);

	    memset(temp, '\0', 10);
	    strncpy(temp, buff + 38, 8);
	    sscanf(temp, "%f", &C[(*np)].coord[natom][1]);

	    memset(temp, '\0', 10);
	    strncpy(temp, buff + 46, 8);
	    sscanf(temp, "%f", &C[(*np)].coord[natom][2]);

	    natom++;
	}
    }
    fclose(IN);
    (*np)++;
    return C;
}

int get_id(char *id, struct cofm *C, int np) {
    int i, j;

    for(i = 0; i < np; ++i) {
	if(strcmp(id, C[i].id) == 0) {
	    return i;
	}
    }
    return -1;
}

float *vector(int nl, int nh) {
        float *v;

        v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
        if(!v) {
	    printf("allocation failure in vector()");
	    exit(-1);
	}
	
        return v-nl;
}

void free_vector(float *v, int nl, int nh) {
        free((char*) (v+nl));
}

void lubksb(float **a, int n, int *indx, float b[]) {
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
		for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];

	}
}

int ludcmp(float **a, int n, int *indx, float *d) {
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv,*vector();
	void free_vector();

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) {
			   big=temp;
			}
		if (big == 0.0)  return -1; /* matrix is singular */
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) {
			   a[i][j] *= dum;
		  	}
		}
	}
	free_vector(vv,1,n);
        return 0;
}

void matinv(float **a, float **y, float *d, int *indx) {

    /* Inverts a 3 x 3 matrix. This routine was written with reference to
     *   reading Numerical Recipies in C. 
     * W.H Press, B.P Flannery, S.A. Teukolsky and W.T. Vetterling, 
     * "Numerical Recipes in C: The Art of Scientific Computing",    
     * Cambridge University Press, 1988. */
 
    int i,j,N;
    float *col;

    col = (float *) malloc(10 * sizeof(float));
    N = 3;

    if(ludcmp(a, N, indx,d ) == -1) { /* matrix is singular, so just copy a to y (ie. I^-1 = I) */
	for(j = 1; j <= N; j++) { 
	    for(i = 1; i <= N; i++) {
		y[i][j] = a[i][j];
	    }
	}
    }
    else {
	for(j = 1; j <= N; j++) {
	    for(i = 1; i <= N; i++) col[i] = 0.0;
	    col[j] = 1.0;
	    lubksb(a, N, indx, col);
	    for(i = 1; i <= N; i++) y[i][j] = col[i];
	}
    }
    free(col);
}

int matvecprod(float **A, float *C, float *B, FILE *OUTPUT) {
    /*
     * Multiplies a vector by a matrix.
     * The sense is: C = A x B
     */

    int i, j ,k;
    float temp;
    float t[3];

    /*
    printf("matvecprod:\n");
    printf("A x V:\n");
    printmat(A,B,3,OUTPUT);
    */
	
    for(j = 0; j < 3; ++j) {
	t[j] = 0.0;
	for(k = 0; k < 3; ++k) {
	    t[j] += A[j][k] * B[k];
	}
    }
    for(i = 0; i < 3; ++i) C[i] =t[i];

    /*
    printf("product:\n");
    for(i=0; i<3; ++i) printf("%8.5f\n",C[i]);
    */

    return 0;

}

int matprod(float **P, float **A, float **B, FILE *OUTPUT) {
    /*
     * multiplies two 3x3 matrices 
     * the sense is P = A x B
     */

    int i, j, k, l;
    float **T;
    float *V;

    T = (float **) malloc(3 * sizeof(float *));
    V = (float *) malloc(3 * sizeof(float));
    for(i = 0; i < 3; ++i) {
	V[i] = 0;
	T[i] = (float *) malloc(3 * sizeof(float));
	for(j = 0; j < 3; ++j) T[i][j] = 0.0;
    }

    for(i = 0; i < 3; ++i) {
	for(j = 0; j < 3; ++j) {
            for(k = 0; k < 3; ++k) {
                T[i][j] += A[i][k] * B[k][j];
	    }
	}
    }

    /*
    printf("matprod:\n");
    printf("A:\n");
    printmat(A, V, 3, OUTPUT);
    printf("x B:\n");
    printmat(B, V, 3, OUTPUT);
    printf("is C:\n");
    printmat(T, V, 3, OUTPUT);
    */

    for(i = 0; i < 3; ++i) {
	for(j = 0; j < 3; ++j) {
	    P[i][j] = T[i][j];
	}
    }

    free(V);
    for(i = 0; i < 3; ++i) free(T[i]);
    free(T);

    return 0;
}

int printmat(float **R, float *V, int n, FILE *OUTPUT) {
        int i,j;
        for(i=0; i<n; ++i) {
           for(j=0; j<n; ++j) fprintf(OUTPUT,"%10.5f ",R[i][j]);
           fprintf(OUTPUT,"         %10.5f \n",V[i]);
        }
        return 0;
}

int back_trans(float **Rnew, float *Vnew, float **R1, float *V1, float **R2, float *V2) {
    // Apply transformation R2/V2 backwards to R1/V1 - that is center on R2/V2
    // Result is returned to Rnew/Vnew
    // Doesn't fucking work
   
    int i,j,k;

    int *indx;
    float sign;
    float *negvec;
    float **invmat;
    float **R,**RI;

    indx = (int *) malloc(100 * sizeof(int));
    invmat = (float **) malloc(3 * sizeof(float *));
    R = (float **) malloc(4 * sizeof(float *));
    RI = (float **) malloc(4 * sizeof(float *));
    negvec = (float *) malloc(3 * sizeof(float));

    for(i = 0; i < 4; ++i) { 
	R[i] = (float *) malloc(4 * sizeof(float));
	RI[i] = (float *) malloc(4 * sizeof(float));
    }
    for(i = 0; i < 3; ++i) {
	invmat[i] = (float *) malloc(3 * sizeof(float));
    }

    /* get the inverse */
    for(i = 0; i < 3; ++i) {
	for(j = 0; j < 3; ++j) {
	    R[i + 1][j + 1] = R2[i][j];
	}
    }
    matinv(R, RI, &sign, indx);

    for(i = 0; i < 3; ++i) { 
	for(j = 0; j < 3; ++j) {
	    invmat[i][j] = RI[i + 1][j + 1];
	}
    }
    /* store the negative of V2 in negvec */
    for(i = 0; i < 3; ++i) {
	negvec[i] = -1 * V2[i];
    }

    /*
    printf("\n");
    printf("The original:\n");
    printmat(R2, V2, 3,stdout);
    printf("The inverse:\n");
    printmat(invmat, negvec, 3, stdout);
    */

    /* Now apply the reverse transformation to each of the matrices 
     *
     * Think, if X is the 
     *    X' = RxX + Vx
     * then,
     *    X = Rx^-1(X'-Vx)
     *
     * where X and X' are the old and new coordinates.  Now we
     * want to apply these to all other transformations, Y, where:
     *
     *  Y' = RyY + Vy
     *
     * Therefor,
     * 
     *  Y'' = Rx^-1(Y'-Vx)
     *      = Rx^-1(RyY + Vy - Vx)
     *
     * We must thus first subtract the translation, then apply the 
     *  inverse matrix, giving the new transformations as:
     *
     *  Ry' = Rx^-1 Ry
     *  Vy' = Rx^-1(Vy - Vx)
     */

    for(j = 0; j < 3; ++j) {
	Vnew[j] = V1[j] + negvec[j];
    }
    matvecprod(invmat, Vnew, Vnew, stdout);

    /* now apply the inverse to the old matrix */
    matprod(Rnew, invmat, R1, stdout);
    //matprod(Rnew, R1, invmat, stdout);
   
    /*
     * Rnew = R1 * invmat
     * Vnew = V1 + negvec
     */

    free(indx);
    free(negvec);
    for(i = 0; i < 3; ++i) { 
	free(invmat[i]); 
    }
    for(i = 0; i < 4; ++i) { 
	free(RI[i]); 
	free(R[i]); 
    }
    free(R);
    free(RI);

    return 1;
}

int transform(float **c, int n, float **R, float *V) {
    int i, j, k;
    float *nc;

    nc = (float *) malloc(3 * sizeof(float));

    for(i = 0; i < n; ++i) {
	for(j = 0; j < 3; ++j) {
	    nc[j] = 0.0;
	}

	for(j = 0; j < 3; ++j) {
	    for(k=0; k<3; ++k) {
		nc[j] += R[j][k] * c[i][k];
	    }
	}
	for(j = 0; j < 3; ++j) {
	    c[i][j] = nc[j] + V[j];
	}
    }
    free(nc);
}

float RMSD(float **c1, float **c2, int n) {
    int i, j, k;
    float D;

    D = 0.0;

    for(i = 0; i < n; ++i) {
	/*
	printf(
	       "Comparing %8.3f %8.3f %8.3f and %8.3f %8.3f %8.3f\n",
	       c1[i][0],
	       c1[i][1],
	       c1[i][2],
	       c2[i][0],
	       c2[i][1],
	       c2[i][2]
	       );
	*/

	for(j = 0; j < 3; ++j) {
	    D += pow((c1[i][j] - c2[i][j]), 2);
	}
    }

    return sqrt(D / (float) n);
}



int main(int argc, char *argv[]) {
    int i, j, k;
    int nt, np;
    int nid1, nid2, nid3, nid4;
    int ntokens;
    int id1_trans, id2_trans, id3_trans, id4_trans;
    int trans1_3;
    int trans2_4;
    float **R1, **R2;
    float **coords1, **coords2;
    float *V1, *V2;
    float iRMSD;
    char transfile[1000];
    char cofmfile[1000];
    char intspairsfile[1000];
    char buff[1000];
    char ** tokens;
    char **doms;
    struct trans_set *TR;
    struct cofm *C;
    struct trans_index *tind;

    FILE *IN;
 
    if(argc != 3) {
	fprintf(stderr, "Usage: iRMSD [transfile] [cofmfile] < [intspairsfile]\n");
	exit(-1);
    }
    strcpy(transfile,argv[1]);
    strcpy(cofmfile,argv[2]);
    //strcpy(intspairsfile,argv[3]);

    C = get_cofm(cofmfile, &np);
    printf("%%Total number of cofm entries:  %10d\n", np);

    TR = get_trans(transfile, &nt);
    printf("%%Total number of trans entries: %10d\n", nt);

    printf("%%Indexing domains in transformations\n");
    tind = (struct trans_index *) malloc(np * sizeof(struct trans_index));
    for(i = 0; i < np; ++i) {
	C[i].nid = i;
	tind[i].n = 0;
	tind[i].t = (int *) malloc(sizeof(int));
	tind[i].p = (int *) malloc(sizeof(int));
    }
    for(i = 0; i < nt; ++i) {
	if(((i + 1) % 100) == 0) {
	    printf("%8d so far\n", i + 1);
	}
	for(j = 0; j < TR[i].n; ++j) {
	    for(k = 0; k < np; ++k) {
		if(strcmp(TR[i].id[j], C[k].id) == 0) {
		    TR[i].nid[j] = C[k].nid;
		    tind[k].t[tind[k].n] = i;
		    tind[k].p[tind[k].n] = j;

		    tind[k].n++;
		    tind[k].t = (int *) realloc(tind[k].t, (tind[k].n + 1) * sizeof(int));
		    tind[k].p = (int *) realloc(tind[k].p, (tind[k].n + 1) * sizeof(int));

		    /*
		    printf(
			   "%s found in trans %d (as %s in trans)\n",
			   C[k].id,
			   i,
			   TR[i].id[j]
			   );
		    */

		    break;
		}
	    }
	}
    }
    printf("\n");

    /*
    for(i = 0; i < nt; ++i) {
	printf("TRANS %d contains %d elements\n", i, TR[i].n);
	for(j = 0; j < TR[i].n; ++j) {
	    printf("%s %s { %s\n", TR[i].f[j], TR[i].id[j], TR[i].des[j]);
	    for(k = 0; k < 3; ++k) {
		printf(
		       "%10.5f %10.5f %10.5f   %10.5f\n",
		       TR[i].t[j].R[k][0],
		       TR[i].t[j].R[k][1],
		       TR[i].t[j].R[k][2],
		       TR[i].t[j].V[k]
		       );
	    }
	}
    }
    printf("\n");
    */

    /*
    if((IN = fopen(intspairsfile,"r")) == NULL) {
	fprintf(stderr, "Error reading file %s\n", intspairsfile);
	exit(-1);
    }
    */

    V1 = (float *) malloc(3 * sizeof(float));
    V2 = (float *) malloc(3 * sizeof(float));

    R1 = (float **) malloc(3 * sizeof(float*));
    R2 = (float **) malloc(3 * sizeof(float*));
    for(i = 0; i < 3; ++i) {
	R1[i] = (float *) malloc(3 * sizeof(float));
	R2[i] = (float *) malloc(3 * sizeof(float));
    }

    coords1 = (float **) malloc(20 * sizeof(float *));
    coords2 = (float **) malloc(20 * sizeof(float *));
    for(i = 0; i < 20; ++i) {
	coords1[i] = (float *) malloc(3 * sizeof(float));
	coords2[i] = (float *) malloc(3 * sizeof(float));
    }

    while((fgets(buff, 999, stdin)) != NULL) {
	tokens = c_split(buff, &ntokens, ' ');
	if(ntokens < 4) { 
	    fprintf(stderr,"Error in interaction files - should be space not tab separated\n");
	    fprintf(stderr, "Line was %s ", buff);
	    exit(-1);
	}
	nid1 = get_id(tokens[0], C, np);
	nid2 = get_id(tokens[1], C, np);
	nid3 = get_id(tokens[2], C, np);
	nid4 = get_id(tokens[3], C, np);

	// printf("\n#################\n");
	printf("%s::%s %s::%s ", tokens[0], tokens[1], tokens[2], tokens[3]);

	if((nid1 != -1) && (nid2 != -1) && (nid3 != -1) && (nid4 != -1)) {
	    /*
	    printf("\n");
	    printf(
		   "Pair of pairs %s %s vs %s %s (%d %d vs %d %d)\n",
		   tokens[0],
		   tokens[1],
		   tokens[2],
		   tokens[3],
		   nid1,
		   nid2,
		   nid3,
		   nid4
		   );
	    printf(
		   "Translated as %s %s vs %s %s (%d %d vs %d %d)\n",
		   C[nid1].id,
		   C[nid2].id,
		   C[nid3].id,
		   C[nid4].id,
		   nid1,
		   nid2,
		   nid3,
		   nid4
		   );

	    printf("Total trans for %s %d\n", tokens[0], tind[nid1].n);
	    printf("Total trans for %s %d\n", tokens[1], tind[nid2].n);
	    printf("Total trans for %s %d\n", tokens[2], tind[nid3].n);
	    printf("Total trans for %s %d\n", tokens[3], tind[nid4].n);

	    for(i = 0; i < tind[nid1].n; ++i) {
		printf(
		       "Id1 %s (%d) found in trans %d pos %d\n",
		       tokens[0],
		       nid1,
		       tind[nid1].t[i],
		       tind[nid1].p[i]
		       );
	    }
	    for(j = 0; j < tind[nid3].n; ++j) {
		printf(
		       "Id3 %s (%d) found in trans %d pos %d\n",
		       tokens[2],
		       nid3,
		       tind[nid3].t[j],
		       tind[nid3].p[j]
		       );
	    }
	    */
	    
	    trans1_3 = -1;
	    for(i = 0; i < tind[nid1].n; ++i) {
		for(j = 0; j < tind[nid3].n; ++j) {
		    if(tind[nid1].t[i] == tind[nid3].t[j]) {
			trans1_3 = tind[nid1].t[i];
			id1_trans = tind[nid1].p[i];
			id3_trans = tind[nid3].p[j];
			break;
		    }
		}
		if(trans1_3 != -1) {
		    break;
		}
	    }
	    if(trans1_3 != -1) {
		//printf("1/3 transformation (%s/%s) is %d\n",tokens[0],tokens[2],trans1_3);
		if(strcmp(C[nid1].id, TR[trans1_3].id[id1_trans]) != 0) {
		    fprintf(
			    stderr,
			    "Error in assignment first id %s given %s in transfile\n",
			    C[nid1].id,
			    TR[trans1_3].id[id1_trans]
			    );
		    exit(-1);
		}

		/*
		printf("%s\n", TR[trans1_3].id[id1_trans]);
		for(i = 0; i < 3; ++i) {
		    printf(
			   "%10.5f %10.5f %10.5f  %10.5f\n",
			   TR[trans1_3].t[id1_trans].R[i][0],
			   TR[trans1_3].t[id1_trans].R[i][1],
			   TR[trans1_3].t[id1_trans].R[i][2],
			   TR[trans1_3].t[id1_trans].V[i]
			   );
		}
		printf("ID real3 %s\n",C[nid3].id);
		*/

		if(strcmp(C[nid3].id, TR[trans1_3].id[id3_trans]) != 0) {
		    fprintf(
			    stderr,
			    "Error in assignment third id %s given %s in transfile\n",
			    C[nid3].id,
			    TR[trans1_3].id[id3_trans]
			    );
		    exit(-1);
		}

		/*
		printf("%s\n", TR[trans1_3].id[id3_trans]);
		for(i = 0; i < 3; ++i) {
		    printf(
			   "%10.5f %10.5f %10.5f  %10.5f\n",
			   TR[trans1_3].t[id3_trans].R[i][0],
			   TR[trans1_3].t[id3_trans].R[i][1],
			   TR[trans1_3].t[id3_trans].R[i][2],
			   TR[trans1_3].t[id3_trans].V[i]
			   );
		}
		*/
	    }
	    else {
		printf("ERR: No 1/3 trans ");
	    }

	    /*
	    for(i = 0; i < tind[nid2].n; ++i) {
		printf(
		       "Id2 %s (%d) found in trans %d pos %d\n",
		       tokens[1],
		       nid2,
		       tind[nid2].t[i],
		       tind[nid2].p[i]
		       );
	    }
	    for(j = 0; j < tind[nid4].n; ++j) {
		printf(
		       "Id4 %s (%d) found in trans %d pos %d\n",
		       tokens[3],
		       nid4,
		       tind[nid4].t[j],
		       tind[nid4].p[j]
		       );
	    }
	    */

	    trans2_4 = -1;
	    for(i = 0; i < tind[nid2].n; ++i) {
		for(j = 0; j < tind[nid4].n; ++j) {
		    if(tind[nid2].t[i] == tind[nid4].t[j]) {
			//printf("Found equivalence trans %d\n", tind[nid2].t[i]);
			trans2_4 = tind[nid2].t[i];
			id2_trans = tind[nid2].p[i];
			id4_trans = tind[nid4].p[j];
			break;
		    }
		}
		if(trans2_4 != -1) {
		    break;
		}
	    }
	    if(trans2_4 != -1) {
		//printf("2/4 transformation is %d\n", trans2_4);
		if(strcmp(C[nid2].id, TR[trans2_4].id[id2_trans]) != 0) {
		    fprintf(
			    stderr,
			    "Error in assignment first id %s given %s in transfile\n",
			    C[nid2].id,
			    TR[trans2_4].id[id2_trans]
			    );
		    exit(-1);
		}

		/*
		printf("%s\n", TR[trans2_4].id[id2_trans]);
		for(i = 0; i < 3; ++i) {
		    printf(
			   "%10.5f %10.5f %10.5f  %10.5f\n",
			   TR[trans2_4].t[id2_trans].R[i][0],
			   TR[trans2_4].t[id2_trans].R[i][1],
			   TR[trans2_4].t[id2_trans].R[i][2],
			   TR[trans2_4].t[id2_trans].V[i]
			   );
		}
		printf("%s\n",TR[trans2_4].id[id4_trans]);
		for(i = 0; i < 3; ++i) {
		    printf(
			   "%10.5f %10.5f %10.5f  %10.5f\n",
			   TR[trans2_4].t[id4_trans].R[i][0],
			   TR[trans2_4].t[id4_trans].R[i][1],
			   TR[trans2_4].t[id4_trans].R[i][2],
			   TR[trans2_4].t[id4_trans].V[i]
			   );
		}
		*/
	    }
	    else {
		printf("ERR: No 2/4 trans");
	    }
	    
	    if((trans1_3 != -1) && (trans2_4 != -1)) {
		// Make two transformations & coordinates

		// Transformation of 1 on to 3
		back_trans(
			   R1,
			   V1,
			   TR[trans1_3].t[id1_trans].R,
			   TR[trans1_3].t[id1_trans].V,
			   TR[trans1_3].t[id3_trans].R,
			   TR[trans1_3].t[id3_trans].V
			   );

		/*
		printf("\n\n");
		printf("%s on to %s\n", TR[trans1_3].id[id1_trans], TR[trans1_3].id[id3_trans]);
		printf("The old:\n");
		printmat(TR[trans1_3].t[id1_trans].R, TR[trans1_3].t[id1_trans].V, 3, stdout);
		printf("The new:\n");
		printmat(R1, V1, 3, stdout);

		printf("\nDom file:\n");
		printf(
		       "%s %s {%s\n",
		       TR[trans1_3].f[id1_trans],
		       TR[trans1_3].id[id1_trans],
		       TR[trans1_3].des[id1_trans]
		       );
		printmat(R1, V1, 3, stdout);
		printf("}\n");
		printf(
		       "%s %s {%s\n",
		       TR[trans2_4].f[id2_trans],
		       TR[trans2_4].id[id2_trans],
		       TR[trans2_4].des[id2_trans]
		       );
		printmat(R1, V1, 3, stdout);
		printf("}\n");
		printf(
		       "%s %s {%s}\n",
		       TR[trans1_3].f[id3_trans],
		       TR[trans1_3].id[id3_trans],
		       TR[trans1_3].des[id3_trans]
		       );
		printf(
		       "%s %s {%s}\n",
		       TR[trans2_4].f[id4_trans],
		       TR[trans2_4].id[id4_trans],
		       TR[trans2_4].des[id4_trans]
		       );
		*/

		// Transformation of 2 on to 4
		back_trans(
			   R2,
			   V2,
			   TR[trans2_4].t[id2_trans].R,
			   TR[trans2_4].t[id2_trans].V,
			   TR[trans2_4].t[id4_trans].R,
			   TR[trans2_4].t[id4_trans].V
			   );

		/*
		printf("\n");
		printf("%s on to %s\n", TR[trans2_4].id[id2_trans], TR[trans2_4].id[id4_trans]);
		printf("The old:\n");
		printmat(TR[trans2_4].t[id2_trans].R, TR[trans2_4].t[id2_trans].V, 3, stdout);
		printf("The new:\n");
		printmat(R2, V2, 3, stdout);
		printf("\n");

		printf("\nDom file:\n");
		printf(
		       "%s %s {%s}\n",
		       TR[trans1_3].f[id1_trans],
		       TR[trans1_3].id[id1_trans],
		       TR[trans1_3].des[id1_trans]
		       );

		printf(
		       "%s %s {%s}\n",
		       TR[trans2_4].f[id2_trans],
		       TR[trans2_4].id[id2_trans],
		       TR[trans2_4].des[id2_trans]
		       );
		printf(
		       "%s %s {%s\n",
		       TR[trans1_3].f[id3_trans],
		       TR[trans1_3].id[id3_trans],
		       TR[trans1_3].des[id3_trans]
		       );
		printmat(R2, V2, 3, stdout);
		printf("}\n");
		printf(
		       "%s %s {%s\n",
		       TR[trans2_4].f[id4_trans],
		       TR[trans2_4].id[id4_trans],
		       TR[trans2_4].des[id4_trans]
		       );
		printmat(R2, V2, 3, stdout);
		printf("}\n");
		*/

		// Transform the centres of mass of 1 and 2 by R1 and V1
		// Transform the centres of mass of 1 and 2 by R2 and V2
		for(i = 0; i < 7; ++i) {
		    for(j = 0; j < 3; ++j) {
			coords1[i][j]     = C[nid1].coord[i][j];
			coords2[i][j]     = C[nid1].coord[i][j];

			coords1[i + 7][j] = C[nid2].coord[i][j];
			coords2[i + 7][j] = C[nid2].coord[i][j];
		    }
		}

		/*
		for(i = 0; i < 14; ++i) {
		    printf("Original %8.3f %8.3f %8.3f\n", coords1[i][0], coords1[i][1], coords1[i][2]);
		}
		*/

		transform(coords1, 14, R1, V1);
		transform(coords2, 14, R2, V2);

		iRMSD = RMSD(coords1, coords2, 14);
		printf("iRMSD %8.5f", iRMSD);
	    } 
	}
	else {
	    printf("ERR: domain indexing");
	}
	printf("\n");
       
	for(i = 0; i < ntokens; ++i) {
	    free(tokens[i]);
	}
	free(tokens);
    }
}
