#include <iostream>
#include "csr_matrix.h"

double** allomatrix(int m){ //Dynamically allocate a matrix;
 double** a;
 a=new double*[m];
 for(int i=0;i<m;i++){
    a[i]=new double[m];
 }
 return a;
}


void dematrix(int m,double **A){//free the storage for a matrix;
 for(int i=0;i<m;i++){
    delete[] A[i];
 }
 delete[] A;
}

double* allovec(int row){ //Dynamically allocate a vector;
    double* a;
    a=new double[row];
    return a;
}

void devec(double* a){ //free the storage for a vector;
 delete[] a;
}

//this function is used to print CSR format of A and output the dimension of A.
csr_matrix CSR(double ** A, int nrow, bool p)
{
 csr_matrix m;
 m.row=nrow;
 m.p=p;

 if(p==0){
 int w=0;
 for(int i=0;i<nrow;i++){
    for(int j=0;j<nrow;j++){
        if(A[i][j]!=0){
                w++;}}}//count the size of matrix_entries and col_n;
 m.matrix_entries=allovec(w);
 m.col_n=allovec(w);
 m.row_start=allovec(nrow+1);
 int q=0;
 int g;
 m.row_start[0]=0;
 for(int i=0;i<nrow;i++){
    for(int j=0;j<nrow;j++){
        if(A[i][j]!=0){
                q++;//once A[i][j] is not zero, q plus 1
                g=q-1;
                m.matrix_entries[g]=A[i][j];
                m.col_n[g]=((double) j);
                }
    }
    m.row_start[i+1]=(double) q;
 }
 return m;
}

 else {
 int w=0;
 for(int i=0;i<nrow;i++){
    for(int j=i;j<nrow;j++){//j=i,regardless of lower triangle part
        if(A[i][j]!=0){
                w++;}}}
 m.matrix_entries=allovec(w);
 m.col_n=allovec(w);
 m.row_start=allovec(nrow+1);
 int q=0;
 int g;
 m.row_start[0]=0;
 for(int i=0;i<nrow;i++){
    for(int j=i;j<nrow;j++){
        if(A[i][j]!=0){
                q++;
                g=q-1;
                m.matrix_entries[g]=A[i][j];
                m.col_n[g]=((double) j);
                }
    }
    m.row_start[i+1]=(double) q;
 }
 return m;
}
}



bool judgesymmetric(int m,double **A){//the function is used to judge whether a matrix is symmetric or not, if so return 1;
 bool p=1;
 for(int i=0;i<m;i++){
    for(int j=i;j<m;j++){
        if(A[i][j]!=A[j][i]){
               p=0;
               return p;}
        else{continue;}
                }
    }
   return p;
 }

//the function is used to compute Ax=b which A was stored as not symmetric CSR format
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
              b[k]+=matrix_entries[i]*x[(int) col_n[i]]; //A[k][(int) col_n[i]] * x[(int) col_n[i]]
              }
          k++;
          q++;
 }
 return b;
}

//the function is used to compute Ax=b which A was stored as symmetric CSR format
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

void printvec(int n,double *a)//print a vector(pointer) to the screen
{
    for (int k=0;k<n;k++)
    {
        std::cout << a[k] << " ";
    }
    std::cout<<std::endl;
}
