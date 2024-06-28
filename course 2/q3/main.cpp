#include <iostream>
#include <cmath>
#include <iomanip>
#include "Vector.hpp"
#include "OrbitODE.hpp"
#include "ForwardEulerSolver.hpp"
#include "SymplecticEulerSolver.hpp"
#include "StoermerVerletSolver.hpp"

const double G=6.674e-11;
const double masse=5.972e24;
const double massm=7.342e22;
const double r0=3.844e8;
const double rearth=6.378e6;
const double rmoon=1.737e6;
const double pi=3.141592653589793;


int main()
{
 Vector xp(3);
 Vector x0(3);
 Vector v0(3);
 double v;
 double T;

 v=G*masse/r0;
 v=pow(v,0.5);
 T=2*pi*pow(r0,1.5)/(pow(G*masse,0.5));

 OrbitODE w;
 w.setm(masse,xp);
 x0(1)=r0;
 v0(2)=v;

 ForwardEulerSolver *e=new ForwardEulerSolver(w,x0,v0,0,T,10000,"c.dat",1,1);
 e->Solve();
 SymplecticEulerSolver *e3=new SymplecticEulerSolver(w,x0,v0,0,T,10000,"c3.dat",1,1);
 e3->Solve();
 StoermerVerletSolver *e4=new StoermerVerletSolver(w,x0,v0,0,T,10000,"c4.dat",1,1);
 e4->Solve();


 v0(2)=0;
 ForwardEulerSolver *e2=new ForwardEulerSolver(w,x0,v0,0,T,10000,"c2.dat",1,1);
 e2->setcollisionr(rearth+rmoon,xp);
 double ct=0;
 Vector cx(3);
 e2->Directcollision(ct,cx);
 std::cout<<std::setprecision(5)<<ct<<std::endl;

}


