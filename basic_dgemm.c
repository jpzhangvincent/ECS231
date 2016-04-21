
void
square_dgemm (const unsigned M, 
              const double *A, const double *B, double *C)
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
