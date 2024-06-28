#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "StoermerVerletSolver.hpp"

//set private members
StoermerVerletSolver::StoermerVerletSolver( ODEInterface & anODESystem ,
 Vector**  initialState ,
 Vector**  initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize ,
const int N,
const std :: string fileName ,
const int saveGap  ,
const int printGap) {

   mFinalTime=finalTime ;
   mInitialTime=initialTime ;
   mpODESystem= &anODESystem;
   mStepSize=stepSize;
   minitialState=initialState;
   minitialVelocity=initialVelocity;
   mfileName=fileName;
   msaveGap=saveGap;
   mprintGap=printGap;
   mN=N;
}

void  StoermerVerletSolver::Solve(){
   double h,N;
   N=mN;
   h=(mFinalTime-mInitialTime)/mStepSize;

   Vector **V;
   Vector **X;
   V=minitialVelocity;
   X=minitialState;

   Vector f(3);
   Vector v(3);
   Vector x(3);
   Vector newv(3);
   Vector newx(3);
   Vector newf(3);
   double t=mInitialTime;

  std::ofstream output_file;
  output_file.setf(std::ios::scientific,std::ios::floatfield);
  // Print format header
  output_file.open(mfileName.c_str());
  assert(output_file.is_open());


int m=0;
for(int i=0;i<mStepSize;i++){
    m=m+1;
    for(int j=0;j<N;j++){
    //get current state and velocity of body j
    x=(*X[j]);
    v=(*V[j]);

    mpODESystem->ComputeF(t,x,f,j);//f(t,x(t))

    newx=x+v*h+f*0.5*pow(h,2.0);//x(t+h)

    mpODESystem->ComputeF(t+h,newx,newf,j);//f(t+h,x(t+h))

    newv=v+(f+newf)*(h*0.5);//v(t+h)

    (*X[j])=newx;//renew state of body j
    mpODESystem->setlocation(newx,j);//renew matrix location in NBodyODE.cpp
    (*V[j])=newv;////renew velocity of body j

    output_file << std::setw(14) << x[0]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[1]<<std::setw(1)<<" ";
    output_file << std::setw(14) << x[2]<<std::setw(1)<<" ";
}
 t=t+h;
 output_file << std::endl;
}
}

//This function will not output data to files.
void  StoermerVerletSolver::Solveforcomputetime(){
   double h,N;
   N=mN;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector **V;
   Vector **X;
   V=minitialVelocity;
   X=minitialState;

   Vector f(3);
   Vector v(3);
   Vector x(3);
   Vector newv(3);
   Vector newx(3);
   Vector newf(3);
   double t=mInitialTime;

for(int i=0;i<mStepSize;i++){
    for(int j=0;j<N;j++){
    x=(*X[j]);
    v=(*V[j]);
    mpODESystem->ComputeF(t,x,f,j);

    newx=x+v*h+f*0.5*pow(h,2.0);
    mpODESystem->ComputeF(t+h,newx,newf,j);
    newv=v+(f+newf)*(h*0.5);
    (*X[j])=newx;
    mpODESystem->setlocation(newx,j);
    (*V[j])=newv;
}
 t=t+h;
}
}


void StoermerVerletSolver::Detectcollision(double coliisiontime, Vector& ritself, Vector **state){
   double h,N;
   N=mN;
   h=(mFinalTime-mInitialTime)/mStepSize;
   Vector **V;
   Vector **X;
   V=minitialVelocity;
   X=minitialState;

   Vector f(3);
   Vector v(3);
   Vector x(3);
   Vector newv(3);
   Vector newx(3);
   Vector newf(3);
   double t=mInitialTime;

   //Compute distance of N bodies
    Vector zN(3);
    double znormN;
    bool p=1;
for(int i=0;i<mStepSize;i++){
    for(int j=0;j<N;j++){

    x=(*X[j]);
    v=(*V[j]);
    mpODESystem->ComputeF(t,x,f,j);

    for(int k=0; k<N;k++){
    if(k!=j){

    zN=(*X[k])-x;//distance between body k and j
    znormN=zN.CalculateNorm(2);

    if(znormN<ritself[k]+ritself[j]){//if distance between body k and j is smaller than the sum of their radius, first collision has happened
        std::cout<<"Collision happens"<<std::endl;
        coliisiontime=t;
        state=X;
        std::cout<<"Collision time is"<<" "<<coliisiontime<<"s"<<std::endl;
        //set p=0 to jump out outer loop.
        p=0;

        break;//jump out inner loop
    }
    }
    }

    if(p==0){break;}//jump out second loop.

    newx=x+v*h+f*0.5*pow(h,2.0);
    mpODESystem->ComputeF(t+h,newx,newf,j);
    newv=v+(f+newf)*(h*0.5);
    (*X[j])=newx;
    mpODESystem->setlocation(newx,j);
    (*V[j])=newv;
}
 if(p==0){break;}//jump out first loop.
 t=t+h;
}
  if(p==1){std::cout<<"Collision has not happened"<<std::endl;}
}

Vector** StoermerVerletSolver::setmatrix(int N){
 Vector **A;
 A=new Vector*[N];
 for(int i=0;i<N;i++){
        A[i]=new Vector(3);
 }
 return A;
}

void StoermerVerletSolver::deletematrix(Vector** A,int N){
 for(int i=0;i<N;i++){
        delete[] A[i];
 }
 delete[] A;
}
