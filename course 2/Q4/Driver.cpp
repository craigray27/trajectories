#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include "Vector.hpp"
#include "NBodyODE.hpp"
#include "StoermerVerletSolver.hpp"

const double G=6.674e-11;

const double masssun=1.9891e30;//mass of Sun

const double masse=5.972e24;//mass of Earth

const double massm=7.342e22;//mass of Moon

const double r0=3.844e8;//distance between earth and moon

const double rearth=6.378e6;//radius of Earth

const double rmoon=1.737e6;//radius of Moon

const double pi=3.141592653589793;

//Allocate matrix
Vector** setmatrix(int N);
//Deallocate matrix
void deletematrix(Vector** A,int N);
//Set parameters of 5 bodies system
void setfivebodies( Vector **X,Vector **V);
//Set parameters of N bodies system
void setNbodies( Vector **X,Vector **V,int N);

int main(int argc, char* argv[])
{
 double collisiontime=0;

 //Earth and Moon system
 int N=2;

 Vector **X;//matrix of initial state includes N Vectors of N bodies
 Vector **V;//matrix of initial velocity includes N Vectors of N bodies
 X=setmatrix(N);
 V=setmatrix(N);

 Vector mass(N);//Vector of mass includes N bodies' mass
 Vector rself(N);//Vector of N bodies' own radius
 rself[0]=rearth;
 rself[1]=rmoon;
 mass[0]=masse;
 mass[1]=massm;

 Vector xearth(3);//initial state of Earth
 Vector xmoon(3);
 Vector veearth(3);
 Vector vmoon(3);//initial velocity of Moon

 double v;
 double T;
 v=G*masse/r0;
 v=pow(v,0.5);
 T=2*pi*pow(r0,1.5)/(pow(G*masse,0.5));
 xmoon(1)=r0;
 vmoon(2)=v;//initial velocity of Moon

 (*X[0])=xearth;
 (*X[1])=xmoon;
 (*V[0])=veearth;
 (*V[1])=vmoon;

 NBodyODE w;
 w.setm(mass,X,N);
 StoermerVerletSolver *e4=new StoermerVerletSolver(w,X,V,0,T,10000,N,"earthandmoon.dat",1,1);
 std::cout<<"In Earth and Moon system"<<std::endl;

 e4->Detectcollision(collisiontime,rself,X);//check whether collision has happened
 e4->Solve();


 //Five bodies system
 N=5;
 T=31536000;//3153600s=365days, period of Earth

 Vector mass5(N);
 Vector ritself5(N);
 Vector **X5;
 Vector **V5;
 X5=setmatrix(N);
 V5=setmatrix(N);
 setfivebodies(X5,V5);

 double rsun=6.966e8;//radius of Sun

 double massmercury=3.3011e23;//mass of Mercury

 double rmercury=2.44e6;//radius of Mercury

 double massvenus=4.8675e24;//mass of Venus

 double rvenus=1.2103e7;

 double massmars=6.4171e23;//mass of Mars

 double rmars=3.3895e6;

 mass5[0]=masssun;
 mass5[1]=massmercury;
 mass5[2]=massvenus;
 mass5[3]=masse;
 mass5[4]=massmars;

 ritself5[0]=rsun;
 ritself5[1]=rmercury;
 ritself5[2]=rvenus;
 ritself5[3]=rearth;
 ritself5[4]=rmars;

 NBodyODE w5;
 w5.setm(mass5,X5,N);
 StoermerVerletSolver *e5=new StoermerVerletSolver(w5,X5,V5,0,T,10000,N,"5bodies.dat",1,1);
 std::cout<<std::endl;
 std::cout<<"In five bodies system"<<std::endl;
 e5->Detectcollision(collisiontime,ritself5,X);
 e5->Solve();

 //N bodies to compute CPU time
 std::cout<<std::endl;
 std::cout<<"Please input N to compute CPU time"<<std::endl;
 std::cin>>N;
 Vector **XN;
 Vector **VN;
 Vector massN(N);
 XN=setmatrix(N);
 VN=setmatrix(N);
 setNbodies(XN,VN,N);
 T=2*pi*pow(r0,1.5)/(pow(G*masse,0.5));

 for(int i=0; i<N-1; i+=2){//set mass of N bodies
    massN[i]=masse;//massN[i] is Earth

    massN[i+1]=massm;//massN[i+1] is Moon
 }

 std::chrono::high_resolution_clock::time_point t0=std::chrono::high_resolution_clock::now() ;
/* Run the section of code which you want timings for */

 NBodyODE wN;
 wN.setm(massN,XN,N);
 StoermerVerletSolver *e6=new StoermerVerletSolver(wN,XN,VN,0,T,100000,N,"Nbodies.dat",1,1);
 e6->Solveforcomputetime();
  // Store the final time point in t1
std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now();
// Report t1 -t0 (in milliseconds )
std::cout<<std::endl << " Time elapsed : "<<std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()<<" milliseconds "<<std::endl;

//Free memory
 deletematrix(X,2);
 deletematrix(V,2);
 deletematrix(X5,5);
 deletematrix(XN,N);

 mass.~Vector();
 rself.~Vector();
 xearth.~Vector();
 xmoon.~Vector();
 veearth.~Vector();
 vmoon.~Vector();
 mass5.~Vector();
 ritself5.~Vector();
 massN.~Vector();

}

Vector** setmatrix(int N){
 Vector **A;
 A=new Vector*[N];
 for(int i=0;i<N;i++){
        A[i]=new Vector(3);
 }
 return A;
}

void deletematrix(Vector** A,int N){
 for(int i=0;i<N;i++){
        delete[] A[i];
 }
 delete[] A;
}

void setfivebodies( Vector **X,Vector **V){
 double Rmercury=5.791e10; //Revolution radius of Mercury in Solar system
 double vmercury=G*masssun/Rmercury;
 vmercury=pow(vmercury,0.5);//initial velocity of Mercury

 double Rvenus=1.08e11;//Revolution radius
 double vvenus=G*masssun/Rvenus;
 vvenus=pow(vvenus,0.5);

 double Rearth=1.5e11;
 double vearth=G*masssun/Rearth;
 vearth=pow(vearth,0.5);

 double Rmars=2.25e11;
 double vmars=G*masssun/Rmars;
 vmars=pow(vmars,0.5);

 //Set initial state and initial velocity with order Sun, Mercury, Venus, Earth and Mars
 (*X[1])[0]=Rmercury;
 (*X[2])[1]=Rvenus;
 (*X[3])[0]=-Rearth;
 (*X[4])[1]=-Rmars;

 (*V[1])[1]=vmercury;
 (*V[2])[0]=vvenus;
 (*V[3])[1]=-vearth;
 (*V[4])[0]=vmars;
}

//I set N bodies system is insisted of K(N=2K) Earth and Moon system, and distance between Earth j and Earth j+1 is 1.5e11.
void setNbodies( Vector **X,Vector **V,int N){
 for(int i=0; i<N-1; i+=2){

    (*X[i])[1]=1.5e11*i;//initial State of Earth i+1
    //initial State of moon i+1
    (*X[i+1])[0]=r0;
    (*X[i+1])[1]=1.5e11*i;

    double velocity;
    velocity=G*masse/r0;
    velocity=pow(velocity,0.5);
    //initial velocity of moon i+1
    (*V[i+1])[1]=velocity;
 }
}

