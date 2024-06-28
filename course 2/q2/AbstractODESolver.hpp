#ifndef ABSTRACTODESOLVER_HPP_INCLUDED
#define ABSTRACTODESOLVER_HPP_INCLUDED

#include "ODEInterface.hpp"

class AbstractODESolver
{
  public:

  virtual void Solve () = 0;
  virtual void eh()=0;

 protected:

   double mFinalTime ;
   double mInitialTime ;
   ODEInterface * mpODESystem ;
   double mStepSize ;
   double minitialState;
   double minitialVelocity;



};



#endif // ABSTRACTODESOLVER_HPP_INCLUDED
