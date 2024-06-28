#include <iostream>
#include <cmath>
#include "linear_algebra q4.h"

//following functions are the same as q3 except algebra and generateAF;


double* allovec(int row){
    double* a;
    a=new double[row];
    return a;
}

void devec(double* a){
 delete[] a;
}

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

//the function is used to compute Ax=b which A was stored as not symmetric format
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

//the function is used to compute Ax=b which A was stored as symmetric format
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


double *plusvecvec(int n, double *b,double *c){
 double *x;
 x=allovec(n);
 for(int i=0;i<n;i++){
    x[i]=b[i]+c[i];
 }
 return x;
}

double *minusvecvec(int n, double *b,double *c){
 double *x;
 x=allovec(n);
 for(int i=0;i<n;i++){
    x[i]=b[i]-c[i];
 }
 return x;
}

//The difference between generatesAF in q3 and generatesAF2 is the latter one is used to generate A with asymmetric CSR format, which is the only format can be used to calculate p^-1 * A directly
void generatesAF2(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y)){

    double h=1.0/(double)(m+1);
    int nm=((m-1)*(3*m-1)+2*m-1)*2-m*m;//Note value of nm is different from generatesAF;
    matrix_entries=allovec(nm);
    col_n=allovec(nm);
    row_start=allovec(m*m+1);

    double x=h;
    double y=h;
    double q=0;
    f[0]=func(x,y);

    int k=0;
    int i=0;
    double element_row_start=3;//u(1,1) is same as the one in generatesAF
    row_start[0]=0;
    row_start[1]=3;
    matrix_entries[0]=4*pow((double)(m+1),2.0)+1.0;
    matrix_entries[1]=-pow((double)(m+1),2.0);
    matrix_entries[2]=-pow((double)(m+1),2.0);
    col_n[0]=q;
    col_n[1]=q+1.0;
    col_n[2]=q+(double) m;
    k++;
    q++;
    i=3;

    while(k<m*m){
    if(k%m!=0){
        x+=h;
        f[k]=func(x,y);
    }
    else{
    x=h;
    y+=h;
    f[k]=func(x,y);}

    if((k<m)&&((k+1)%m!=0)){//In first row of A except entry next to right boundary, u(1,j) is determined by 5 points including itself but u(0,j)=0
        element_row_start+=4;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+2]=-pow((double)(m+1),2.0);
        matrix_entries[i+3]=-pow((double)(m+1),2.0);
        col_n[i]=q-1;
        col_n[i+1]=q;
        col_n[i+2]=q+1;
        col_n[i+3]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=4;
        k++;
        q++;
    }
    else if((k<m)&&((k+1)%m==0)){//u(1,n-1) is determined by 3 points because u(1,n)=0 as well;
        element_row_start+=3;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+2]=-pow((double)(m+1),2.0);
        col_n[i]=q-1;
        col_n[i+1]=q;
        col_n[i+2]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=3;
        k++;
        q++;
    }

    else if((k>=m)&&(k%m==0) && (k<m*m-m)){//u(i,1) i!=n-1 is determined by 4 points as u(i,n)=0;
        element_row_start+=4;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+2]=-pow((double)(m+1),2.0);
        matrix_entries[i+3]=-pow((double)(m+1),2.0);
        col_n[i]=q-(double) m;
        col_n[i+1]=q;
        col_n[i+2]=q+1;
        col_n[i+3]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=4;
        k++;
        q++;
    }

    else if((k>=m)&&(k%m!=0) &&((k+1)%m!=0)&&(k<m*m-m)){//u(i,j) i!=n-1, j!=1,j!=n-1 is determined by u(i-1,j),u(i,j),u(i,j-1),u(i,j+1),u(i+1,j)
        element_row_start+=5;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        matrix_entries[i+2]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+3]=-pow((double)(m+1),2.0);
        matrix_entries[i+4]=-pow((double)(m+1),2.0);
        col_n[i]=q-(double) m;
        col_n[i+1]=q-1;
        col_n[i+2]=q;
        col_n[i+3]=q+1;
        col_n[i+4]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=5;
        k++;
        q++;
    }

     else if((k>=m)&&((k+1)%m==0)&&(k<m*m-m)){//u(i,n-1), i!=n-1; u(i,n)=0
       element_row_start+=4;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        matrix_entries[i+2]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+3]=-pow((double)(m+1),2.0);
        col_n[i]=q-(double) m;
        col_n[i+1]=q-1;
        col_n[i+2]=q;
        col_n[i+3]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=4;
        k++;
        q++;
    }

    else if((k%m==0)&&(k>=m*m-m)){//u(n-1,1);u(n,1)=0;u(n-1,0)=0
        element_row_start+=3;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+2]=-pow((double)(m+1),2.0);
        col_n[i]=q-(double) m;
        col_n[i+1]=q;
        col_n[i+2]=q+1;
        row_start[k+1]=element_row_start;
        i+=3;
        k++;
        q++;
    }

        else if((k>m*m-m)&&((k+1)%m!=0)){//u(n-1,j) j!=1
        element_row_start+=4;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        matrix_entries[i+2]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+3]=-pow((double)(m+1),2.0);
        col_n[i]=q-(double) m;
        col_n[i+1]=q-1;
        col_n[i+2]=q;
        col_n[i+3]=q+1;
        row_start[k+1]=element_row_start;
        i+=4;
        k++;
        q++;
    }

    else if(((k+1)%m==0) && (k==m*m-1)){//u(n-1,n-1);
        element_row_start+=3;
        matrix_entries[i]=-pow((double)(m+1),2.0);
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        matrix_entries[i+2]=4*pow((double)(m+1),2.0)+1.0;
        col_n[i]=q-(double) m;
        col_n[i+1]=q-1;
        col_n[i+2]=q;
        row_start[k+1]=element_row_start;
        i+=3;
        k++;
        q++;
    }
    }
}

//This function is used to calculate p^-1 * A and A should be stored in asymmetric format because matrix p^-1 * A would be asymmetric once any two diagonal entries of A is different.
void multiplyPA(int n,double * matrix_entries, double * col_n,double * row_start,double *b){
 double *p;
 p=allovec(n);
 int d=0;
 while(d<n){
     for(int i=(int) row_start[d];i<(int) row_start[d+1];i++){
        if((int) col_n[i]==d){//find diagonal entry of row d, A[d][d]
            p[d]=matrix_entries[i];
            b[d]=b[d]/matrix_entries[i];//P^-1 * b

        }
     }
       d++;
 }
 d=0;
 while(d<n){
     for(int i=(int) row_start[d];i<(int) row_start[d+1];i++){
            matrix_entries[i]=matrix_entries[i]/p[d];//In row d, A[d][j]=A[d][j]/A[d][d]
     }
                 d++;
 }
}

//Algebraq4 function is same as q2, the only difference is function pointer should be MultiplyMatrixVec1 because A is stored as asymmetric format
double *algebraq4(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
                int n, double * matrix_entries, double * col_n,double * row_start,double *b,double *x0,double tol){

 int k=0;
 double *x;
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

 x=x0;
 a1=multiplymatrixvec( matrix_entries,  col_n, row_start,x0,n);
 r=minusvecvec(n,b,a1);
 p=r;
 double t=multiplyvecvec(n,r,r);
 std::cout<<"2-norm of gamma for step "<<k+1<<" "<<"is "<<t<<std::endl;

 while(t>tol){
    b1=multiplymatrixvec(matrix_entries,  col_n, row_start,p,n);
    m=multiplyvecvec(n,p,b1);
    alpha=t/m;
    x=plusvecvec(n,x,multiplyconstvec(n,alpha,p));
    r1=minusvecvec(n,r,multiplyconstvec(n,alpha,b1));
    beta=multiplyvecvec(n,r1,r1)/t;
    p=plusvecvec(n,r1,multiplyconstvec(n,beta,p));
    r=r1;
    t=multiplyvecvec(n,r,r);
    k++;
    std::cout<<"2-norm of gamma for step "<<k+1<<" "<<"is "<<t<<std::endl;

 }
 std::cout<<std::endl;
 std::cout<<"The iteration counter number is "<<k<<std::endl;
 std::cout<<std::endl;
 return x;
}


void printv(int n,double *a)
{
    for (int k=0;k<n;k++)
    {
        std::cout << a[k] << " ";
    }
    std::cout<<std::endl;
}



//following functions will not used in main.cpp of q4, they are used for more general situation that the diagonal entries of A are totally not equal and A is stored as symmetric CSR format.
//In this situation, matrix P^-1 * A will no longer be a symmetric matrix, so the only method is find the A of asymmetric CSR format(A is symmetric but it should be regarded as asymmetric)
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
                w++;}}}
 m.matrix_entries=allovec(w);
 m.col_n=allovec(w);
 m.row_start=allovec(nrow+1);
 int q=0;
 int g;
 m.row_start[0]=0;
 for(int i=0;i<nrow;i++){
    for(int j=0;j<nrow;j++){
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

 else {
 int w=0;
 for(int i=0;i<nrow;i++){
    for(int j=i;j<nrow;j++){
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

//This function is used to obtain matrix A from its symmetric CSR format
double **inversecsr(int n,double * matrix_entries, double * col_n,double * row_start){
 double **A;
 A=allomatrix(n);
 for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        A[i][j]=0;
    }
 }
 int k=0;
 int col;
 while(k<n){
    for(int i=(int) row_start[k];i<(int) row_start[k+1];i++){//In row k, A's upper triangle part has several entries form matrix_entries[row_start[k]] to matrix_entries[row_start[k+1]]
        col=(int) col_n[i];
        A[k][col]=matrix_entries[i];
        if(col!=k){
            A[col][k]=matrix_entries[i];//entries in lower triangle part A[i][j]=A[j][i];
        }
    }
    k++;
 }
 return A;
}
//After inversecsr(), using CSR() function to get A stored as asymmetric CSR format, then continue with algebraq4() and MultiplyMatrixVec1();
