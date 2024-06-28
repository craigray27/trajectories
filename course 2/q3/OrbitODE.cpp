#include <iostream>
#include <cmath>
#include "OrbitODE.hpp"

// Override default constructor
OrbitODE::OrbitODE(){
 mass=0.0;
 xp=new Vector(3);
 for(int i=1; i<=3;i++){
    (*xp)(i)=0.0;
 }
}

void OrbitODE::ComputeF( const double t, const Vector& x,  Vector& f ) const
{
    Vector z(3);
    double znorm;
    z=(*xp)-x;
    znorm=z.CalculateNorm(2);//distance between earth and moon
    f=z*(G*mass/pow(znorm,3.0));
}

void OrbitODE::ComputeAnalyticSolution( const double t, Vector& x ) const
{
}

//set members
void OrbitODE:: setm(double a,Vector& p){
 mass=a;
 xp=new Vector(p);
}


