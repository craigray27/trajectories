#ifndef LINEAR_ALGEBRA_Q4_H_INCLUDED
#define LINEAR_ALGEBRA_Q4_H_INCLUDED

struct csr_matrix;
struct csr_matrix
{
 int row;
 double *matrix_entries;
 double *col_n;
 double *row_start;
 bool p;
};
double* allovec(int row);
void devec(double *a);
double** allomatrix(int m);
void dematrix(int m,double **A);
double* MultiplyMatrixVec1(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double* MultiplyMatrixVec2(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double multiplyvecvec(int n, double *b,double *c);
double *multiplyconstvec(int n, double alpha,double*b);
double *plusvecvec(int n, double *b,double *c);
double *minusvecvec(int n, double *b,double *c);
void multiplyPA(int n,double * matrix_entries, double * col_n,double * row_start,double *b);
double *algebraq4(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
                int n, double * matrix_entries, double * col_n,double * row_start,double *b,double *x0,double tol);

void generatesAF2(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y));
void printvec(int n,double *a);

csr_matrix CSR(double ** A, int nrow, bool p);
double **inversecsr(int n,double * matrix_entries, double * col_n,double * row_start);


#endif // LINEAR_ALGEBRA_Q4_H_INCLUDED
