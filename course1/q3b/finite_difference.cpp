#include <iostream>
#include <cmath>
#include "finite_difference.h"

//following functions are from q3 except generatesAF.
double* allovec(int row){
    double* a;
    a=new double[row];
    return a;
}

void devec(double* a){
 delete[] a;
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
              b[(int) col_n[i]]+=matrix_entries[i]*x[k];}
              }
          k++;
          q++;
 }
 return b;
}


double multiplyvecvec(int n, double *b,double *c){
 double alpha=0;
 for(int i=0;i<n;i++){
    alpha+=b[i]*c[i];
 }
 return alpha;
}


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

//The function pointer would be MultiplyMatrixVec2 because  it is easy find out that A is a symmetric matrix, but situation is totally different in q4;
double *algebra(double *(*multiplymatrixvec)(double *matrix_entries, double *col_n,double *row_start,double *tu,int nrow),
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

//This function is generating A stored as symmetric CSR format and RHS f;
void generatesAF(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y)){//m=n-1 where n is steps and m is dimension of matrix A

    double h=1.0/(double)(m+1);
    int nm=(m-1)*(3*m-1)+2*m-1;//It is the relationship between the size of matrix entries which is stored as symmetric format and the dimension of A,
                                //it can be observed easily by reading following code
    matrix_entries=allovec(nm);
    col_n=allovec(nm);
    row_start=allovec(m*m+1);
    //initialize RHS f;
    double x=h;
    double y=h;
    double q=0;
    f[0]=func(x,y);

    int k=0;
    int i=0;
    double element_row_start=3;//In the first entries of first row of A: u(1,1) is determined by 5 entries but 2 are zero(boundary) so only 3 will be displayed on CSR
    row_start[0]=0;
    row_start[1]=3;
    matrix_entries[0]=4*pow((double)(m+1),2.0)+1.0;
    matrix_entries[1]=-pow((double)(m+1),2.0);
    matrix_entries[2]=-pow((double)(m+1),2.0);
    col_n[0]=q;//u(1,1) itself
    col_n[1]=q+1.0;//u(1,2)
    col_n[2]=q+(double) m;//u(1,8)
    k++;
    q++;
    i=3;
    //compute RHS f
    while(k<m*m){
    if(k%m!=0){
        x+=h;
        f[k]=func(x,y);//f(k*h,h)
    }
    else{
    x=h;
    y+=h;//if k%m=0, row+1,so x should be set to h again and y should plus h
    f[k]=func(x,y);}
    //generates CSR format of A, you can see the other format in q4;
    if(((k+1)%m!=0) && (k<m*m-m)){
        element_row_start+=3;//Every u(i,j) except those next to the right boundary (u(i,m))) and the last row(u(n-1,j)), will display 3 number in CSR format because the u(i-1,j) and u(i,j-1) will be ignored in symmetric format
        matrix_entries[i]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        matrix_entries[i+2]=-pow((double)(m+1),2.0);
        col_n[i]=q;//u(i,j) itself
        col_n[i+1]=q+1;//u(i,j+1)
        col_n[i+2]=q+(double) m;//u(i+1,j)
        row_start[k+1]=element_row_start;
        i+=3;
        k++;
        q++;
    }
    else if(((k+1)%m==0) && (k<m*m-m)){
        element_row_start+=2;//every u(i,m) next to right boundary except last row, will display 2 number in CSR because u(i,m+1)=u(i,n)=0;
        matrix_entries[i]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        col_n[i]=q;
        col_n[i+1]=q+(double) m;
        row_start[k+1]=element_row_start;
        i+=2;
        k++;
        q++;

    }
    else if((((k+1)%m!=0) && (k>m*m-m)) or (k==m*m-m)){
        element_row_start+=2;//every u(m,j)=u(n-1,j) will display 2 number because u(m+1,j)=u(n,j)=0;
        matrix_entries[i]=4*pow((double)(m+1),2.0)+1.0;
        matrix_entries[i+1]=-pow((double)(m+1),2.0);
        col_n[i]=q;
        col_n[i+1]=q+1.0;
        row_start[k+1]=element_row_start;
        i+=2;
        k++;
        q++;
    }
    else if(((k+1)%m==0) && (k==m*m-1)){
        element_row_start+=1;//u(m-1,m-1) is only determined by itself as the u(m-2,m-1) and u(m-1,m-2) will not shown, the other two are zero
        matrix_entries[i]=4*pow((double)(m+1),2.0)+1.0;
        col_n[i]=q;
        row_start[k+1]=element_row_start;
        k++;
    }
    }



}

void printv(int n,double *a)
{
    for (int k=0;k<n;k++)
    {
        std::cout << a[k] << " ";
    }
    std::cout<<std::endl;
}
