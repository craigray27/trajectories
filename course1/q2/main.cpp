#include <iostream>
#include <vector>
#include "linear_algebra.h"
#include <fstream>

int main()
{
 //used for read from files
 double* matrix_entries;
 double* col_n;
 double* row_start;
 std::vector<double> xi;
 double *x;
 int noRows;
 int nrows,nm,nr;
 bool p;
 //read from files
 x=ReadVector("vector2.dat",noRows);
 ReadMatrix("matrix2.dat",matrix_entries,col_n,row_start,p,nr,nrows,nm);

 double tol=1e-10;
 double *xhat;
 double *x0;
 double *b;
 //allocate
 x0=allovec(nrows);
 b=allovec(nrows);
 xhat=allovec(nrows);
 xi=std::vector<double> (nrows);

 if(p){//p decides what function should be used to compute matrix stored as CSR format multiplies a vector;
 b=MultiplyMatrixVec2(matrix_entries,col_n,row_start,x,nrows);
 x0=&xi[0];
 xhat=algebra(MultiplyMatrixVec2,nrows,matrix_entries,col_n,row_start,b,x0,tol);//if symmetric, the function pointer should be MultiplyMatrixVec2
 std::cout<<"The approximation to x (xhat) is"<<std::endl;
 printvec(nrows,xhat);
 std::cout<<std::endl;
 //compute error
 double *delta;
 double e;
 delta=allovec(nrows);
 delta=minusvecvec(nrows,x,xhat);
 e=multiplyvecvec(nrows,delta,delta);
 std::cout<<"the error is "<<e<<std::endl;

 devec(matrix_entries);
 devec(col_n);
 devec(row_start);
 devec(x);
 devec(delta);
 }

 else{
 b=MultiplyMatrixVec1(matrix_entries,col_n,row_start,x,nrows);
 x0=&xi[0];
 xhat=algebra(MultiplyMatrixVec1,nrows,matrix_entries,col_n,row_start,b,x0,tol);
 std::cout<<"The approximation to x (xhat) is"<<std::endl;
 printvec(nrows,xhat);
 std::cout<<std::endl;
 double *delta;
 double e;
 delta=allovec(nrows);
 delta=minusvecvec(nrows,x,xhat);
 e=multiplyvecvec(nrows,delta,delta);
 std::cout<<"the error is "<<e<<std::endl;
 devec(matrix_entries);
 devec(col_n);
 devec(row_start);
 devec(x);
 devec(delta);}

  return 0;
}
