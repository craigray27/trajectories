#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cassert>
#include "ForwardEulerSolver.hpp"

//set private members
 ForwardEulerSolver:: ForwardEulerSolver ( ODEInterface & anODESystem ,const Vector & initialState ,const Vector & initialVelocity ,
const double initialTime ,const double finalTime ,const double stepSize ,const std :: string fileName ,const int saveGap  ,const int printGap ) {

   mFinalTime=finalTime ;
   mInitialTime=initialTime ;
   mpODESystem= &anODESystem;
   mStepSize=stepSize;
   minitialState=new Vector(initialState);
   minitialVelocity=new Vector(initialVelocity);
   mfileName=fileName;
   msaveGap=saveGap;
   mprintGap=printGap;
}

//set collisionr an xp
void ForwardEulerSolver::setcollisionr(double r, Vector& p){
 collisionr=r;
 xp=new Vector(p);
}

void  ForwardEulerSolver::Solve(){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector v(3);
   Vector x(3);
   Vector f(3);
   v=*minitialVelocity;
   x=*minitialState;
   double t=mInitialTime;

  std::ofstream output_file;
  output_file.setf(std::ios::scientific,std::ios::floatfield);
  // Print format header
  output_file.open(mfileName.c_str());
  assert(output_file.is_open());
    output_file << std::setw(24) << "x"<<std::setw(1)<<" ";
    output_file << std::setw(24) << "y"<<std::setw(1)<<" ";
    output_file << std::setw(24) << "z"<<std::endl;
    output_file << std::setw(14) << x[0]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[1]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[2]<<std::endl;

   //ForwardEulerSolver method
   for(int i=0;i<mStepSize;i++){
    mpODESystem->ComputeF(t,x,f);//compute f(t,x)
    x=x+v*h;
    v=v+f*h;
    t=t+h;
    output_file << std::setw(14) << x[0]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[1]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[2]<<std::endl;
}
}

//find out whether collisin happens or not and get the time and the state of the system at the end of the time-step during which the collision took place.
void  ForwardEulerSolver::Detectcollision(double& collisiontime, Vector& collisionx){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector v(3);
   Vector x(3);
   Vector f(3);
   Vector u(3);
   v=*minitialVelocity;
   x=*minitialState;
   double t=mInitialTime;

   u=*xp-x;
   double normx=u.CalculateNorm(2);

   int k=1;
   while(normx>collisionr){
    mpODESystem->ComputeF(t,x,f);
    x=x+v*h;
    v=v+f*h;
    t=t+h;
    u=*xp-x;
    normx=u.CalculateNorm(2);//the distance between earth and moon in time tk
    k++;

    if(k>mStepSize){//if collision has not happened during the period, break loop
        break;
    }
   }

   if(k>mStepSize){

   if(h>60){//set h=60s to find nearest minute
   t=t-h;
   v=v-f*h;
   x=x-v*h;
   h=60;
   u=*xp-x;
   normx=u.CalculateNorm(2);
   //continue loop
   while(normx>collisionr){
    mpODESystem->ComputeF(t,x,f);
    x=x+v*h;
    v=v+f*h;
    t=t+h;
    u=*xp-x;
    normx=u.CalculateNorm(2);
   }
   collisiontime=t/60;//collision time
   collisionx=x;//state of system when collision happens
   }
   //if h<60s, then output collision time directly
   else{
   collisiontime=t/60;
   collisionx=x;
   }}

   else{std::cout<<"The collision will not happen within period"<<std::endl;}
}


