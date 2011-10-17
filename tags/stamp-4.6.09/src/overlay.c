#include <stdio.h>
#define max(A,B) ((A) > (B) ? (A) : (B))
#define leq(A,B) ((A) <= (B) ? (A) : (B))
#define max4(A,B,C,D) ((A)>(B))?(((A)>(C))?(((A)>(D))?(A):(D)):(((C)>(D))?(C):(D))):(((B)>(C))?(((B)>(D))?(B):(D)):(((C)>(D))?(C):(D)))

/* This was just a little piddly routine I wrote to figure something out
 *  for STAMP */
int main(int argc, char *argv[]) {

	int M,N,Step;
	int P,count;
	int Q1over,Q2over,D1over,D2over;
	int Q1fit,Q2fit,D1fit,D2fit;
	int i,Min_Size;
	int add;

	float Trunc;

	if(argc<5) { 
	   fprintf(stderr,"error: format overlay <M> <N> <Step> <Trunc>\n");
	   exit(-1);
	}
	
	sscanf(argv[1],"%d",&M);
	sscanf(argv[2],"%d",&N);
	sscanf(argv[3],"%d",&Step);
	sscanf(argv[4],"%f",&Trunc);

	P=0-N+Step;

	count=1;
	if(N>=M) Min_Size=(int)(Trunc*(float)M);
	else Min_Size=(int)(Trunc*(float)N);
	printf("M = %4d, N = %4d, Step = %4d, Trunc = %5.3f, Min_Size = %4d\n", 
	   M,N,Step,Trunc,Min_Size);
	while(P<M) {
	  D1over=max(0,P);
	  D2over=leq(M,P+N);
	  Q1over=max(0,-1*P);
	  Q2over=leq(N,M-P);

	  /* the range of positions in the overlay */
	  printf("Overhang: %3d P %4d Q %4d - %4d; D %4d - %4d\n",
	    count,P,Q1over,Q2over,D1over,D2over);
	  /* now we need to get the fraction to use 
	   *  in the fit */
	  Q1fit=Q1over; Q2fit=Q2over;
	  D1fit=D1over; D2fit=D2over;
	  if((Q2fit-Q1fit)<Min_Size) { 
	     add=(int)((float)(Min_Size-(Q2fit-Q1fit))/2);
	     printf("Q2fit-Q1fit: %4d; Min_Size: %4d; add: %4d\n",
		Q2fit-Q1fit,Min_Size,add);
	     Q1fit-=(add+1); Q2fit+=add;
	     if(Q1fit<0) { 
	       Q2fit=Q2fit-Q1fit; 
	       Q1fit=0;
	     }
	     if(Q2fit>(N)) {
	       Q1fit=Q1fit-(Q2fit-(N));
	       Q2fit=N;
	     }
	     if(Q1fit<0) Q1fit=0;
	   }

	  
          if((D2over-D1over)<Min_Size) { 
             add=(int)((float)(Min_Size-(D2over-D1over))/2);
             D1fit-=(add+1); D2fit+=add;
             if(D1fit<0) { 
               D2fit-=D1fit; 
               D1fit=0;
             }
             if(D2fit>(M)) {
               D1fit-=(D2fit-(M));
               D2fit=M;
             }
	     if(D1fit<0) D1fit=0;
           }
	  printf("To fit:   %3d P %4d Q %4d - %4d; D %4d - %4d\n",
		      count,P,Q1fit,Q2fit,D1fit,D2fit);
	  for(i=0; i<P; ++i) if(i%10==0) printf(" ");
	  for(i=0; i<N; ++i) if(i%10==0) printf("q");
	  printf("\n");
	  for(i=P; i<0; ++i) if(i%10==0) printf(" ");
	  for(i=0; i<M; ++i) if(i%10==0) printf("d");
	  printf("\n");
	  P+=Step; count++;
	}

	exit(0);
}

