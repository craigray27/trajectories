
#include <iostream>
#include <vector>
#include <cmath>
#include "mex.h"

const double pi=3.1415927;
double* allovec(int row);
void generatesAF(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y));
double Func(double x,double y);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{

    int m=8;
    int nm=(m-1)*(3*m-1)+2*m-1;
    double* matrix_entries;
    double* col_n;
    double* row_start;
    double *f;
    plhs[0] = mxCreateDoubleMatrix( 1, nm, mxREAL); //第一个输出是一个5*6的矩阵
     //获得矩阵的第一个元素的指
    matrix_entries = mxGetPr(plhs[0]);
        f=allovec(m*m);
    generatesAF(m, matrix_entries,col_n,row_start,f,Func);

}

void generatesAF(int m,double*& matrix_entries, double*& col_n,double*& row_start,
                 double *f,double(* func)(double x,double y)){//m=n-1 where n is steps and m is dimension of matrix A

    double h=1.0/(double)(m+1);
    int nm=(m-1)*(3*m-1)+2*m-1;//It is the relationship between the size of matrix entries which is stored as symmetric format and the dimension of A,
                                //it can be observed easily by reading following code
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


double* allovec(int row){
    double* a;
    a=new double[row];
    return a;
}

double Func(double x,double y){
 double r;
 r=(2.0*pow(pi,2.0)+1.0)*sin(pi*x)*sin(pi*y);
 return r;
}
