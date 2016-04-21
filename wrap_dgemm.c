
void
square_dgemm(const unsigned M, 
             const double *A, const double *B, double *C)
{
  extern int dgemm_(const char* transA, const char* transB,
                    const unsigned* M, const unsigned *N, const unsigned* K,
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
}

