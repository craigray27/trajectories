#ifndef SYMPLECTICEULERSOLVER_HPP_INCLUDED
#define SYMPLECTICEULERSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"



class SymplecticEulerSolver: public AbstractODESolver
{
public:
    SymplecticEulerSolver( ODEInterface & anODESystem ,
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

#endif // SYMPLECTICEULERSOLVER_HPP_INCLUDED
