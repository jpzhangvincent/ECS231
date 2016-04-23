#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void square_dgemm (const unsigned M, const double *A, const double *B, double *C)
{
  unsigned i, j, k;

  for (i = 0; i < M; ++i) {
       for (j = 0; j < M; ++j) {
            const double *Ai_ = A + i;
            const double *B_j = B + j*M;

            double cij = *(C + j*M + i);

            for (k = 0; k < M; ++k) {
                 cij += *(Ai_ + k*M) * *(B_j + k);
            }

            *(C + j*M + i) = cij;
       }
  }
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