#ifndef STOERMERVERLETSOLVER_HPP_INCLUDED
#define STOERMERVERLETSOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "AbstractODESolver.hpp"


class StoermerVerletSolver: public AbstractODESolver
{
public:
    StoermerVerletSolver( ODEInterface & anODESystem ,
 Vector**  initialState ,
 Vector**  initialVelocity ,
const double initialTime ,
const double finalTime ,
const double stepSize ,
const int N,
const std :: string fileName =" output . dat ",
const int saveGap = 1 ,
const int printGap = 1);

void Solve();

void Solveforcomputetime();//used to minimize CPU time(it will not output data to a file)

void Detectcollision(double coliisiontime, Vector& ritself, Vector **state);

Vector** setmatrix(int N);
void deletematrix(Vector** A,int N);

private:
    int mN;
    std::string mfileName;
    int msaveGap;
    int mprintGap;

};

#endif // STOERMERVERLETSOLVER_HPP_INCLUDED
