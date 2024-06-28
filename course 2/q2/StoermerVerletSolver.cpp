#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "StoermerVerletSolver.hpp"

//set private members
StoermerVerletSolver::StoermerVerletSolver ( ODEInterface & anODESystem ,
const double initialState ,
const double initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize,
const std::string fileName,const int saveGap ,
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

void  StoermerVerletSolver::Solve(){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   double v=minitialVelocity;
   double x=minitialState;
   double t=mInitialTime;
   double f=0;
   double newv, newx,newf;
   // Open file
  std::ofstream output_file;
  output_file.setf(std::ios::scientific,std::ios::floatfield);
  // Print format header
  output_file.open(mfileName.c_str());
  assert(output_file.is_open());

  output_file <<std::setw(12) << "t"<<std::setw(1)<<" ";
  output_file << std::setw(24) << "v"<<std::setw(1)<<" ";
  output_file << std::setw(24) << "x"<<std::endl;
  std::cout<<std::setw(5)<<"t"<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<"v"<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<"x"<<std::endl;

 //print t0,v(0),x(0)
  std::cout<<std::setw(5)<<t<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<v<<std::setw(1)<<" ";
  std::cout<<std::setw(10)<<x<<std::endl;
    output_file << std::setw(14) << t<<std::setw(1)<<" ";
    output_file << std::setw(14) << v<<std::setw(1)<<" ";
    output_file << std::setw(14) << x<<std::endl;

    //StoermerVerletSolver method
    int k=0;
   for(int i=0;i<mStepSize;i++){//Print t,v,x to screen every 100 steps. 100=time gap/h.
    mpODESystem->ComputeF(t,x,f);

    newx=x+v*h+0.5*pow(h,2.0)*f;//compute x(t+h)

    mpODESystem->ComputeF(t+h,newx,newf);//compute f(t+h,x(t+h)) by using newx

    newv=v+(h/2)*(f+newf);//compute v(t+h) by using v(t), f(t,x(t)) and f(t+h,x(t+h))

    t=t+h;
    x=newx;
    v=newv;
    k=k+1;
    if(k%100==0){//Print t,v,x to screen every 100 steps. 100=time gap/h.
    std::cout<<std::setw(5)<<t<<std::setw(1)<<" ";
    std::cout<<std::setw(10)<<v<<std::setw(1)<<" ";
    std::cout<<std::setw(10)<<x<<std::endl;}

    output_file << std::setw(14) << t<<std::setw(1)<<" ";
    output_file << std::setw(14) << v<<std::setw(1)<<" ";
    output_file << std::setw(14) << x<<std::endl;
    }

}

//This function is used to compute E(h) with different h.
void  StoermerVerletSolver::eh(){
   double h;
   double v,x,t,f,realx,e,maxe, newv, newx,newf;
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
    realx=0;//Analytic solution
    e=0;
    maxe=0;

   for(int i=0;i<mStepSize;i++){
    mpODESystem->ComputeF(t,x,f);
    newx=x+v*h+0.5*pow(h,2.0)*f;
    mpODESystem->ComputeF(t+h,newx,newf);
    newv=v+(h/2)*(f+newf);
    t=t+h;
    x=newx;
    v=newv;
    mpODESystem->ComputeAnalyticSolution(t,realx);
    e=x-realx;

    if(fabs(e)>maxe){
        maxe=fabs(e);//find E(h)
    }
    }
    output_file<<std::setw(5)<<h<<std::setw(1)<<" ";
    output_file << std::setw(10) << log(h)<<std::setw(1)<<" ";
    output_file << std::setw(10) << maxe<<std::setw(1)<<std::endl;

    mStepSize=mStepSize-10;//reduce stepsize and increase h
   }
}
