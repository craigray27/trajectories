#ifndef STOERMERVERLETSOLVER_HPP_INCLUDED
#define STOERMERVERLETSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"


class StoermerVerletSolver: public AbstractODESolver
{
public:
    StoermerVerletSolver( ODEInterface & anODESystem ,
const Vector & initialState ,
const Vector & initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize ,
const std :: string fileName =" output . dat ",
const int saveGap = 1 ,
const int printGap = 1)  ;

void setcollisionr(double r, Vector& p);
void setprint(bool print);//r is the distance between earth and moon when collision happens and p is location of earth
void Solve () ;
void Detectcollision(double& collisiontime, Vector& collisionx);

private:
    std::string mfileName;
    int msaveGap;
    int mprintGap;
    Vector* xp;//location of earth
    double collisionr;//distance between earth and moon when collision happens
    bool mprint; //it decides whether print to screen. It will not print data when p=0;

};

#endif // STOERMERVERLETSOLVER_HPP_INCLUDED
