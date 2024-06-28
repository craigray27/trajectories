#include <iostream>
#include <cmath>
#include "OscillatorODE.hpp"
#include "SymplecticEulerSolver.hpp"
#include "StoermerVerletSolver.hpp"



int main(int argc, char* argv[]){
 double a=1.5;
 OscillatorODE t;
 t.seta(a);

 //SymplecticEulerSolver
 SymplecticEulerSolver *z=new SymplecticEulerSolver(t,0,1.5,0,30,3000,"SymplecticEulerSolver.dat");
 std::cout<<"The results of SymplecticEulerSolver"<<std::endl;
 z->Solve();
 //compute E(h)
 SymplecticEulerSolver *z2=new SymplecticEulerSolver(t,0,1.5,0,1,10000,"SymplecticEulerSolverEh.dat");
 z2->eh();

 //StoermerVerletSolver
 std::cout<<std::endl;
 std::cout<<"The results of StoermerVerletSolver"<<std::endl;
 StoermerVerletSolver *q=new StoermerVerletSolver(t,0,1.5,0,30,3000,"StoermerVerletSolver.dat");
 q->Solve();

 StoermerVerletSolver *q2=new StoermerVerletSolver(t,0,1.5,0,1,10000,"StoermerVerletSolverEh.dat");
 q2->eh();



}
