#include <iostream>
#include "linear_algebra.h"
#include <cassert>
#include <fstream>

double* allovec(int row){
    double* a;
    a=new double[row];
    return a;
}

void devec(double* a){
 delete[] a;
}

//the next two function is from q1, so there is no comment here
double* MultiplyMatrixVec1(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow){
 double *b;
 b=allovec(nrow);
 for(int i=0;i<nrow;i++){
    b[i]=0.0;
 }
 int k=0;
 int q=1;
 while(k<nrow){
          double n=row_start[k];
          double m=row_start[q];
          for(int i=(int) n;i<(int) m;i++){
              b[k]+=matrix_entries[i]*x[(int) col_n[i]]; //A[k][(int) col_n[i]]*x[(int) col_n[i]]
              }
          k++;
          q++;
 }
 return b;
}


double* MultiplyMatrixVec2(double *matrix_entries, double *col_n, double *row_start,double * x,int nrow){
 double *b;
 b=allovec(nrow);
 for(int i=0;i<nrow;i++){
    b[i]=0.0;
 }
 int k=0;
 int q=1;
 while(k<nrow){
          double n=row_start[k];
          double m=row_start[q];
          for(int i=(int) n;i<(int) m;i++){
              b[k]+=matrix_entries[i]*x[(int) col_n[i]];
              if((int) col_n[i]!=k){
              b[(int) col_n[i]]+=matrix_entries[i]*x[k];}//A[(int) col_n[i]][k]*x[k], which is on the lower part of A and will not be illustrated in CSR format
              }
          k++;
          q++;
 }
 return b;
}

//compute result of b*c(b,c are two vectors with same length n), where two-norm is the special case when b==c;
double multiplyvecvec(int n, double *b,double *c){
 double alpha=0;
 for(int i=0;i<n;i++){
    alpha+=b[i]*c[i];
 }
 return alpha;
}

//compute result of alpha*b(b is a vectors and alpha is a constant number)
double *multiplyconstvec(int n, double alpha,double*b){
 double*c;
 c=allovec(n);
 for(int i=0;i<n;i++){
    c[i]=b[i]*alpha;
 }
 return c;
}

//compute b+c(b,c are two vectors with same length n)
double *plusvecvec(int n, double *b,double *c){
 double *x;
 x=allovec(n);
 for(int i=0;i<n;i++){
    x[i]=b[i]+c[i];
 }
 return x;
}

//compute b-c(b,c are two vectors with same length n)
double *minusvecvec(int n, double *b,double *c){
 double *x;
 x=allovec(n);
 for(int i=0;i<n;i++){
    x[i]=b[i]-c[i];
 }
 return x;
}

//this function is used to read vector
double* ReadVector(std::string filename, int& noRows)
{
  // Open the file
  std::ifstream readFile(filename);
  assert(readFile.is_open());

// Read dimensions
  std::string lineAsString;
  std::string val;
  std::getline(readFile,lineAsString);
  readFile >> noRows;

  double* v = new double[noRows];

  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  for (int k=0;k<noRows;k++)
  {
    readFile >> v[k];
  }

  readFile.close();

  return v;

}

//this function is used to read matrix stored as csr format
void ReadMatrix(std::string filename, double*& matrix_entries,double*& col_n,double*& row_start,bool& p,int &nr,int &row,int& nm)
{
  std::ifstream readFile(filename);
  assert(readFile.is_open());

  std::string lineAsString;
  std::string val;
  std::getline(readFile,lineAsString);
  readFile >> p;//read symmetric or not

//dimension of matrix
  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  readFile >> row;

  nr=row+1;//dimension of row_start

  //dimension of matrix_entries
  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  readFile >> nm;

  matrix_entries = allovec(nm);
  col_n = allovec(nm);
  row_start = allovec(nr);

  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  for (int k=0;k<nr;k++)
  {
    readFile >> row_start[k];
  }

  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  for (int k=0;k<nm;k++)
  {
    readFile >> col_n[k];
  }

  std::getline(readFile,lineAsString);
  std::getline(readFile,lineAsString);
  for (int k=0;k<nm;k++)
  {
    readFile >> matrix_entries[k];
  }


  readFile.close();

}

//this function is used to compute xhat which is approximation to x(Ax=b; x=A^(-1)b)
double *algebra(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
                int n, double * matrix_entries, double * col_n,double * row_start,double *b,double *x0,double tol){

 int k=0;
 double *x;//it is xhat actually;
 double *r;
 double *r1;
 double *p;
 double *a1;
 double *b1;
 double alpha,beta,m;
 x=allovec(n);
 r=allovec(n);
 r1=allovec(n);
 p=allovec(n);
 a1=allovec(n);
 b1=allovec(n);

 x=x0;//x0 should set to be equal to zero vector
 a1=multiplymatrixvec( matrix_entries,  col_n, row_start,x0,n);
 r=minusvecvec(n,b,a1);
 p=r;
 double t=multiplyvecvec(n,r,r);//2-norm of gamma0
 std::cout<<"2-norm of gamma for step "<<k+1<<" "<<"is "<<t<<std::endl;//step=k+1
 while(t>tol){
    b1=multiplymatrixvec(matrix_entries,  col_n, row_start,p,n);//the function pointer depends on p(whether A is symmetric or not)
    m=multiplyvecvec(n,p,b1);
    alpha=t/m;
    x=plusvecvec(n,x,multiplyconstvec(n,alpha,p));
    r1=minusvecvec(n,r,multiplyconstvec(n,alpha,b1));
    beta=multiplyvecvec(n,r1,r1)/t;
    p=plusvecvec(n,r1,multiplyconstvec(n,beta,p));
    r=r1;
    t=multiplyvecvec(n,r,r);//2-norm of gamma(k)
    k++;
    std::cout<<"2-norm of gamma for step "<<k+1<<" "<<"is "<<t<<std::endl;

 }
 std::cout<<std::endl;
 std::cout<<"The iteration counter number is "<<k<<std::endl;
 std::cout<<std::endl;
 return x;
}

void printvec(int n,double *a)//print a vector(pointer) to the screen
{
    for (int k=0;k<n;k++)
    {
        std::cout << a[k] << " ";
    }
    std::cout<<std::endl;
}

