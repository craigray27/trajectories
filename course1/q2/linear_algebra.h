#ifndef LINEAR_ALGEBRA_H_INCLUDED
#define LINEAR_ALGEBRA_H_INCLUDED

#include <fstream>

double* allovec(int row);
void devec(double *a);
double* MultiplyMatrixVec1(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double* MultiplyMatrixVec2(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow);
double multiplyvecvec(int n, double *b,double *c);
double *multiplyconstvec(int n, double alpha,double*b);
double *plusvecvec(int n, double *b,double *c);
double *minusvecvec(int n, double *b,double *c);
double* ReadVector(std::string filename, int& noRows);
void ReadMatrix(std::string filename, double*& matrix_entries,double*& col_n,double*& row_start,bool& p,int &nr,int &row,int& nm);
double *algebra(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
                int n, double * matrix_entries, double * col_n,double * row_start,double *b,double *x0,double tol);
void printvec(int n,double *a);


#endif // LINEAR_ALGEBRA_H_INCLUDED
