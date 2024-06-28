#ifndef ODEINTERFACEHEADERDEF
#define ODEINTERFACEHEADERDEF

#include "Vector.hpp"

class ODEInterface
{

  public:

    // Compute right-hand side (pure virtual)
    virtual void ComputeF( const double t, const Vector& x, Vector& f ) const = 0;

    // Compute analytical solution
    virtual void ComputeAnalyticSolution( const double t,  Vector& x ) const;

};

#endif
