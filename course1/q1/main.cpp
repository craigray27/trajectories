#include <iostream>
#include "csr_matrix.h"

void setmatrix(double **A,double *a);//set matrix A(used to be stored as CSR format)

int main(){
 double** A=allomatrix(4);
 double* a=allovec(4);
 double *b=allovec(4);
 setmatrix(A,a);
 bool p=judgesymmetric(4,A);//if p==1, you can also store A as a symmetric format;
 csr_matrix m=CSR(A,4,0);
 int w=0;
 for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
        if(A[i][j]!=0){
                w++;}}}
 std::cout<<"matrix_entries is ";
 printvec(w,m.matrix_entries);
 std::cout<<"column_no is ";
 printvec(w,m.col_n);
 std::cout<<"row_start is ";
 printvec(5,m.row_start);
 std::cout<<std::endl;
 b=MultiplyMatrixVec1(m.matrix_entries,m.col_n,m.row_start,a,4);
 std::cout<<"The result of Ax(using nonsymmetric format) is"<<std::endl;
 for(int i=0;i<4;i++){
        std::cout<<b[i]<<std::endl;
 }
 std::cout<<std::endl;
 if(p==0){
      std::cout<<"matrix A is not symmetric"<<std::endl;
 }
 else{
    bool print;
    std::cout<<"A is symmetric, please input 1 if you want to use symmetric format"<<std::endl;
    std::cin>>print;
    if(print){
        csr_matrix m=CSR(A,4,print);
        int w=0;
        for(int i=0;i<4;i++){
             for(int j=i;j<4;j++){
               if(A[i][j]!=0){
                 w++;}}}
        std::cout<<"matrix_entries is ";
        printvec(w,m.matrix_entries);
        std::cout<<"column_no is ";
        printvec(w,m.col_n);
        std::cout<<"row_start is ";
        printvec(5,m.row_start);
        b=MultiplyMatrixVec2(m.matrix_entries,m.col_n,m.row_start,a,4);
        std::cout<<std::endl;
        std::cout<<"The result of Ax(using symmetric format) is"<<std::endl;
        for(int i=0;i<4;i++){
          std::cout<<b[i]<<std::endl;
          }
 }
 }

 dematrix(4,A);
 devec(a);
 devec(b);

 return 0;
}

void setmatrix(double **A,double *a){
 for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
        A[i][j]=0.0;
    }
  }
 A[0][0]=8.0;
 A[0][3]=2.0;
 A[1][1]=3.0;
 A[1][3]=1.0;
 A[2][2]=4.0;
 A[3][0]=2.0;
 A[3][1]=1.0;
 A[3][3]=7.0;
 a[0]=6.0;
 a[1]=8.0;
 a[2]=2.0;
 a[3]=5.0;

}
