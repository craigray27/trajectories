#ifndef FINITE_DIFERENCE_H_INCLUDED
#define FINITE_DIFERENCE_H_INCLUDED

double* allovec(int row);
void devec(double *a);
double* MultiplyMatrixVec2(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double multiplyvecvec(int n, double *b,double *c);
double *multiplyconstvec(int n, double alpha,double*b);
double *plusvecvec(int n, double *b,double *c);
double *minusvecvec(int n, double *b,double *c);
double *algebra(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
                int n, double * matrix_entries, double * col_n,double * row_start,double *b,double *x0,double tol);

void generatesAF(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y));
void printvec(int n,double *a);


#endif // FINITE_DIFERENCE_H_INCLUDED
