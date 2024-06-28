#ifndef FORWARDEULERSOLVER_HPP_INCLUDED
#define FORWARDEULERSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"



class ForwardEulerSolver: public AbstractODESolver
{
public:
  ForwardEulerSolver ( ODEInterface & anODESystem ,
const Vector & initialState ,
const Vector & initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize ,
const std :: string fileName =" output . dat ",
const int saveGap = 1 ,
const int printGap = 1) ;

void setcollisionr(double r, Vector& p);//r is the distance between earth and moon when collision happens and p is location of earth

void Solve () ;
void Detectcollision(double& collisiontime, Vector& collisionx);

private:
    std::string mfileName;
    int msaveGap;
    int mprintGap;
    double collisionr;//distance between earth and moon when collision happens
    Vector* xp;//location of earth


};

#endif // FORWARDEULERSOLVER_HPP_INCLUDED
