#include <iostream>
#include <cmath>
#include "NBodyODE.hpp"

NBodyODE::NBodyODE(){
    N=0;
}

void NBodyODE::ComputeF( const double t, const Vector& x,  Vector& f, const int j ) const
{
    Vector z(3);
    Vector m(N);

//This Vector is used to renew f(t,x(t))=0 in every step in StoermerVerletSolver because we will sum them in code below so that f of last steps will affect this step
    Vector initialf(3);

    m=*mass;
    double znorm;
    for(int i=0; i<N;i++){
    if(i!=j){
    z=(*location[i])-x;

    znorm=z.CalculateNorm(2);//distance between body i and j

    initialf=initialf+z*(G*m[i]/pow(znorm,3.0));
    }
    }
    f=initialf;
}

//set Vector of mass and matrix of initial location, number=N
void NBodyODE:: setm(Vector& NBodymass,Vector** &Nlocation,double number){
 mass=new Vector(NBodymass);
 N=number;
 location=setmatrix(N);
 location=Nlocation;
}

void NBodyODE::setlocation(Vector& p,const int j){
(*location[j])=p;//renew the location of body j
}

Vector** NBodyODE::setmatrix(int N){
 Vector **A;
 A=new Vector*[N];
 for(int i=0;i<N;i++){
        A[i]=new Vector(3);
 }
 return A;
}

