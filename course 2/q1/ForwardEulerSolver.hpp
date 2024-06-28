#ifndef FORWARDEULERSOLVER_HPP_INCLUDED
#define FORWARDEULERSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"



class ForwardEulerSolver: public AbstractODESolver
{
public:
    ForwardEulerSolver( ODEInterface & anODESystem ,
const double initialState ,
const double initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize,
const std::string fileName="output.dat",const int saveGap = 1 ,
const int printGap = 1) ;

void Solve () ;
//computeEh() is used to save h, E(h) to a file
void computeEh();

private:
    std::string mfileName;
    int msavegap;
    int mprintgap;


};

#endif // FORWARDEULERSOLVER_HPP_INCLUDED
