#ifndef ABSTRACTODESOLVER_HPP_INCLUDED
#define ABSTRACTODESOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "Vector.hpp"

class AbstractODESolver
{
  public:

  virtual void Solve () = 0;


 protected:

   double mFinalTime ;
   double mInitialTime ;
   ODEInterface * mpODESystem ;
   double mStepSize ;
   Vector **minitialState;
   Vector **minitialVelocity;



};

#endif
