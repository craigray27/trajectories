#ifndef CSR_MATRIX_H_INCLUDED
#define CSR_MATRIX_H_INCLUDED

struct csr_matrix;
struct csr_matrix
{
 int row;
 double *matrix_entries;
 double *col_n;
 double *row_start;
 bool p;//if p==0, then matrix is stored as a not symmetric format. if p==1, it is stored as a symmetric format;
};
csr_matrix CSR(double ** A, int nrow, bool p);
bool judgesymmetric(int m,double **A);
double** allomatrix(int m);
void dematrix(int m,double **A);
double* allovec(int row);
void devec(double *a);
double* MultiplyMatrixVec1(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double* MultiplyMatrixVec2(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
void printvec(int n,double *a);


#endif // CSR_MATRIX_H_INCLUDED
