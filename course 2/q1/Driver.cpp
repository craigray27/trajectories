#include <iostream>
#include <cmath>
#include "OscillatorODE.hpp"
#include "ForwardEulerSolver.hpp"



int main(int argc, char* argv[]){
 double a=1.5;
 OscillatorODE t;
 t.seta(a);
 ForwardEulerSolver *z=new ForwardEulerSolver(t,0,1.5,0,30,3000,"ForwardEulerSolver.dat");
 z->Solve();
 ForwardEulerSolver *z2=new ForwardEulerSolver(t,0,1.5,0,1,10000,"Eh.dat");
 z2->computeEh();
}
