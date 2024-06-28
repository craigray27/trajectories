#include <iostream>
#include <vector>
#include <cmath>
#include "finite_difference.h"

const double pi=3.14159265359;//pi

double Func(double x,double y);//RHS function

int main()
{
    int n;
    std::cout<<"please input n= ";
    std::cin>>n;
    int m=n-1;
    double* matrix_entries;
    double* col_n;
    double* row_start;
    double *f;
    f=allovec(m*m);
    generatesAF(m, matrix_entries,col_n,row_start,f,Func);

    double *u;
    u=allovec(m*m);
    int k=0;
    double h=1.0/(double)(m+1);
    double x=h;
    double y=h;
    //Following code is generating the accurate solution to pde(sin(pi*x)*sin(pi*y))
    u[0]=sin(pi*x)*sin(pi*y);
    k++;
    while(k<m*m){
    if(k%m!=0){
        x+=h;
        u[k]=sin(pi*x)*sin(pi*y);
        k++;
    }
    else{x=h;
    y+=h;
    u[k]=sin(pi*x)*sin(pi*y);
    k++;}}

    double *uhat;
    uhat=allovec(m*m);
    double *x0;
    x0=allovec(m*m);
    std::vector<double> xi(m*m);
    x0=&xi[0];
    double tol=1e-10;
    uhat=algebra(MultiplyMatrixVec2,m*m,matrix_entries,col_n,row_start,f,x0,tol);
    //Observe max value of uhat;
    double umax=uhat[0];
    for(int i=1;i<m*m;i++){
        if (uhat[i]>umax){
            umax=uhat[i];
        }
    }
    std::cout<<"max value of uhat is "<<umax<<std::endl;
    std::cout<<std::endl;
    //calculate the error between accurate solution and approximation
    double *e;
    e=allovec(m*m);
    e=minusvecvec(m*m,u,uhat);//vector e=u-uhat
    double ee;
    ee=multiplyvecvec(m*m,e,e);//e^T * e
    std::cout<<"the error between accurate solution and approximation is "<<ee<<std::endl;

    devec(matrix_entries);
    devec(col_n);
    devec(row_start);
    devec(f);
    devec(u);
    devec(uhat);
    devec(x0);
    devec(e);

    return 0;

}

double Func(double x,double y){
 double r;
 r=(2.0*pow(pi,2.0)+1.0)*sin(pi*x)*sin(pi*y);
 return r;
}
