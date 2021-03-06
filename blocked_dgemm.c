#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
/*
  In case you're wondering, dgemm stands for Double-precision, GEneral
  Matrix-Matrix multiplication.
*/

/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/

void
basic_dgemm (const unsigned lda,
             const unsigned M, const unsigned N, const unsigned K,
             const double *A, const double *B, double *C)
{
  unsigned i, j, k;

  /*
    To optimize this, think about loop unrolling and software
    pipelining.  Hint:  For the majority of the matmuls, you
    know exactly how many iterations there are (the block size)...
  */

  for (i = 0; i < M; ++i) {
       const double *Ai_ = A + i;
       for (j = 0; j < N; ++j) {
            const double *B_j = B + j*lda;

            double cij = *(C + j*lda + i);

            for (k = 0; k < K; ++k) {
                 cij += *(Ai_ + k*lda) * *(B_j + k);
            }

            *(C + j*lda + i) = cij;
       }
  }
}

/* You'll definitely change this... */
/*#define BLOCK_SIZE ((unsigned) 5) */

void
do_block (const unsigned lda,
          const double *A, const double *B, double *C,
          const unsigned i, const unsigned j, const unsigned k, int BLOCK_SIZE)
{
     /*
       Remember that you need to deal with the fringes in each
       dimension.

       If the matrix is 7x7 and the blocks are 3x3, you'll have 1x3,
       3x1, and 1x1 fringe blocks.

             xxxoooX
             xxxoooX
             xxxoooX
             oooxxxO
             oooxxxO
             oooxxxO
             XXXOOOX

       You won't get this to go fast until you figure out a `better'
       way to handle the fringe blocks.  The better way will be more
       machine-efficient, but very programmer-inefficient.
     */
     const unsigned M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
     const unsigned N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
     const unsigned K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);

     basic_dgemm (lda, M, N, K,
                  A + i + k*lda, B + k + j*lda, C + i + j*lda);
}

void
square_dgemm (const unsigned M, 
              const double *A, const double *B, double *C, int bs)
{
     const unsigned n_blocks = M / bs + (M%bs? 1 : 0);
     unsigned bi, bj, bk;

     for (bi = 0; bi < n_blocks; ++bi) {
          const unsigned i = bi * bs;
          
          for (bj = 0; bj < n_blocks; ++bj) {
               const unsigned j = bj * bs;

               for (bk = 0; bk < n_blocks; ++bk) {
                    const unsigned k = bk * bs;
                    
                    do_block (M, A, B, C, i, j, k,bs);
               }
          }
     }
}



int main(int argc, char *argv[])
{
  int N = atoi(argv[1]);
  int block_size = atoi(argv[2]);
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
  square_dgemm(N,A,B,C,block_size);
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC;

  double perf = 2*pow((double)N,3)/time_taken;
  printf("%d,%d,%f\n",block_size, N, perf);

  return 0;
}

