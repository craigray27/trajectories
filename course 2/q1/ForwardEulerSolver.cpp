#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "ForwardEulerSolver.hpp"

//set private members
 ForwardEulerSolver::ForwardEulerSolver ( ODEInterface & anODESystem ,const double initialState ,const double initialVelocity ,
const double initialTime ,const double finalTime ,const double stepSize,const std::string fileName,const int saveGap ,
const int printGap) {

   mFinalTime=finalTime ;
   mInitialTime=initialTime ;
   mpODESystem= &anODESystem;
   mStepSize=stepSize;
   minitialState=initialState;
   minitialVelocity=initialVelocity;
   mfileName=fileName;
   msavegap=saveGap;
   mprintgap=printGap;
}

void  ForwardEulerSolver::Solve(){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   int printgap=100;
   double v=minitialVelocity;
   double x=minitialState;
   double t=mInitialTime;
   double f=0;
   // Open file
  std::ofstream output_file;
  output_file.setf(std::ios::scientific,std::ios::floatfield);
  // Print format header
  output_file.open(mfileName.c_str());
  assert(output_file.is_open());

  //print t0, v(0), x(0)
  output_file <<std::setw(12) << "t"<<std::setw(1)<<" ";
  output_file << std::setw(24) << "v"<<std::setw(1)<<" ";
  output_file << std::setw(24) << "x"<<std::endl;
  std::cout<<std::setw(5)<<"t"<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<"v"<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<"x"<<std::endl;
   std::cout<<std::setw(5)<<t<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<v<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<x<<std::endl;
    output_file << std::setw(14) << t<<std::setw(1)<<" ";
    output_file << std::setw(14) << v<<std::setw(1)<<" ";
    output_file << std::setw(14) << x<<std::endl;

   //ForwardEulerSolver method
   int k=0;
   for(int i=0;i<mStepSize;i++){
    x=x+h*v;
    mpODESystem->ComputeF(t,x,f);//compute f(t,x)
    v=v+h*f;
    t=t+h;
    k=k+1;
    if( k%printgap==0){//Print t,v,x to screen every 100 steps. 100=time gap/h.
    std::cout<<std::setw(5)<<t<<std::setw(1)<<" ";
    std::cout<<std::setw(10)<<v<<std::setw(1)<<" ";
    std::cout<<std::setw(10)<<x<<std::endl;}
    output_file << std::setw(14) << t<<std::setw(1)<<" ";
    output_file << std::setw(14) << v<<std::setw(1)<<" ";
    output_file << std::setw(14) << x<<std::endl;
    }

}

//This function is used to compute E(h) with different h.
void  ForwardEulerSolver::computeEh(){
   double h;
   double v,x,t,f,realx,e,maxe;
  std::ofstream output_file;
  output_file.setf(std::ios::scientific,std::ios::floatfield);

  output_file.open(mfileName.c_str());
  assert(output_file.is_open());
  output_file<<std::setw(12)<<"h"<<std::setw(1)<<" ";
  output_file<<std::setw(24)<<"log(h)"<<std::setw(1)<<" ";
  output_file<<std::setw(24)<<"E(h)"<<std::endl;

   while(mStepSize>100){//while h is smaller than 1/100=0.01, and h0=1/10000=0.0001
   h=(mFinalTime-mInitialTime)/mStepSize;
    v=minitialVelocity;
    x=minitialState;
    t=mInitialTime;
    f=0;
    realx=0;//realx is analytic solution to the ODE;
    maxe=0;

   for(int i=0;i<mStepSize;i++){
    mpODESystem->ComputeF(t,x,f);
    x=x+h*v;
    v=v+h*f;
    t=t+h;
    mpODESystem->ComputeAnalyticSolution(t,realx);
    e=x-realx;
    if(fabs(e)>maxe){
        maxe=fabs(e);//maxe is E(h)
    }
    }
    output_file<<std::setw(5)<<h<<std::setw(1)<<" ";
    output_file << std::setw(10) << log(h)<<std::setw(1)<<" ";
    output_file << std::setw(10) << maxe<<std::setw(1)<<std::endl;
    mStepSize=mStepSize-10;//reduce stepsize and increase value of h.
   }
}





