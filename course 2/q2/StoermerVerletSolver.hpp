#ifndef STOERMERVERLETSOLVER_HPP_INCLUDED
#define STOERMERVERLETSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"



class StoermerVerletSolver: public AbstractODESolver
{
public:
    StoermerVerletSolver( ODEInterface & anODESystem ,
const double initialState ,
const double initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize,
const std::string fileName="output.dat",const int saveGap=1,
const int printGap=1) ;

void Solve () ;
void eh();

private:
    std::string mfileName;
    int msavegap;
    int mprintgap;



};


#endif // STOERMERVERLETSOLVER_HPP_INCLUDED
