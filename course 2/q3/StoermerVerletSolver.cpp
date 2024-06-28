#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "StoermerVerletSolver.hpp"

//set private members
StoermerVerletSolver::StoermerVerletSolver( ODEInterface & anODESystem ,
const Vector & initialState ,
const Vector & initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize ,
const std :: string fileName ,
const int saveGap  ,
const int printGap) {

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

void StoermerVerletSolver::setcollisionr(double r, Vector& p){
 collisionr=r;
 xp=new Vector(p);
}

//it decides whether print to screen. It will not print data when p=0;
void StoermerVerletSolver::setprint(bool print){
 mprint=print;
}

void  StoermerVerletSolver::Solve(){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector v(3);
   Vector x(3);
   Vector f(3);

   Vector newv(3);//v(t+h)
   Vector newx(3);//x(t+h)
   Vector newf(3);//f(t+h,(x(t+h))

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
    //whether print to screen
    if(mprint==1){
    std::cout << std::setw(14) << "t"<<std::setw(1)<<" ";
    std::cout << std::setw(14) << "x"<<std::setw(1)<<" ";
    std::cout<< std::setw(14) << "y"<<std::setw(1)<<" ";
    std::cout << std::setw(14) << "z"<<std::endl;
    std::cout << std::setw(14) << t<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[0]<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[1]<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[2]<<std::endl;}

    output_file << std::setw(14) << x[0]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[1]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[2]<<std::endl;

    //StoermerVerletSolver
   int k=0;
   for(int i=0;i<mStepSize;i++){
    mpODESystem->ComputeF(t,x,f);//f(t,x(t))

    newx=x+v*h+f*0.5*pow(h,2.0);//x(t+h)

    mpODESystem->ComputeF(t+h,newx,newf);//f(t+h,x(t+h))

    newv=v+(f+newf)*(h*0.5);//v(t+h)

    t=t+h;

    x=newx;
    v=newv;

    k=k+1;
    if(mprint==1){
    if(k%200==0){//print data to screen, time gap =200h
    std::cout << std::setw(14) << t<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[0]<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[1]<<std::setw(1)<<" ";
    std::cout << std::setw(14) << x[2]<<std::endl;
    }}
    output_file << std::setw(14) << x[0]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[1]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[2]<<std::endl;
}
}

void  StoermerVerletSolver::Detectcollision(double& collisiontime, Vector& collisionx){
   double h;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector v(3);
   Vector x(3);
   Vector f(3);
   Vector u(3);
   Vector newv(3);
   Vector newx(3);
   Vector newf(3);
   v=*minitialVelocity;
   x=*minitialState;
   double t=mInitialTime;

   u=*xp-x;
   double normx=u.CalculateNorm(2);

   int k=1;
   while(normx>collisionr){
    mpODESystem->ComputeF(t,x,f);
    newx=x+v*h+f*0.5*pow(h,2.0);
    mpODESystem->ComputeF(t+h,newx,newf);
    newv=v+(f+newf)*(h*0.5);
    t=t+h;
    x=newx;
    v=newv;
    u=*xp-x;
    normx=u.CalculateNorm(2);//the distance between earth and moon in time tk
    k++;

    if(k>mStepSize){//if collision has not happened during the period, break loop
        break;
    }
   }

   if(k<mStepSize){
   if(h>60){//set h=60s to find nearest minute
   t=t-h;
   v=v-(f+newf)*(h*0.5);
   x=x-v*h-f*0.5*pow(h,2.0);;
   h=60;
   u=*xp-x;
   normx=u.CalculateNorm(2);
    //continue loop
   while(normx>collisionr){
    mpODESystem->ComputeF(t,x,f);
    newx=x+v*h+f*0.5*pow(h,2.0);
    mpODESystem->ComputeF(t+h,newx,newf);
    newv=v+(f+newf)*(h*0.5);
    t=t+h;
    x=newx;
    v=newv;
    u=*xp-x;
    normx=u.CalculateNorm(2);
   }

   collisiontime=t/60;
   collisionx=x;
   }

   else{
   collisiontime=t/60;
   collisionx=x;
   }}
   else{std::cout<<"The collision will not happen within period"<<std::endl;}
}
