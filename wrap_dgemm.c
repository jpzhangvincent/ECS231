#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*void square_dgemm(const int M, const double *A, const double *B, double *C)
{
  extern int dgemm_(const char* transA, const char* transB,
                    const int* M, const int *N, const int* K,
                    const double* alpha, 
                    const double* A, const int* ldA,
                    const double* B, const int* ldB,
                    const double* beta,
                    double* C, const int* ldC,
                    int* ltransA, int* ltransB);

  char trans = 'N';
  double one = 1.0;
  double zero = 0.0;
  int ione = 1;

  dgemm_(&trans, &trans,   
         &M, &M, &M,
         &one,
         A, &M,   B, &M,
         &one,
         C, &M, 
         &ione, &ione);
}*/

#define DGEMM dgemm_
extern void DGEMM (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 

const char* dgemm_desc = "Reference dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are N-by-N matrices stored in column-major format.
 * On exit, A and B maintain their input values.    
 * This function wraps a call to the BLAS-3 routine DGEMM, via the standard FORTRAN interface - hence the reference semantics. */
void square_dgemm (int N, double* A, double* B, double* C)
{
  char TRANSA = 'N';
  char TRANSB = 'N';
  int M = N;
  int K = N;
  double ALPHA = 1.;
  double BETA = 1.;
  int LDA = N;
  int LDB = N;
  int LDC = N;
  DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
} 

int main(int argc, char *argv[])
{
  int N = atoi(argv[1]);
  double *A = (double *)malloc(N * N * sizeof(double));
  double *B = (double *)malloc(N * N * sizeof(double));
  double *C = (double *)malloc(N * N * sizeof(double));
  unsigned i, j;
  
  int seed = time(NULL);
  srand(seed);

  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
        A[i*N+j] = rand();
        B[i*N+j] = rand();
    }
  }

  clock_t t;
  t = clock();
  square_dgemm(N,A,B,C);
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;

  double perf = 2*pow((double)N,3)/time_taken;
  printf("%f\n",perf);

  return 0;
}

