#ifndef ODEINTERFACE_HPP_INCLUDED
#define ODEINTERFACE_HPP_INCLUDED

class ODEInterface
{

  public:

    // Compute right-hand side (pure virtual)
    virtual void ComputeF( const double t, const double x, double& f ) const = 0;

    // Compute analytical solution
    virtual void ComputeAnalyticSolution( const double t, double &x ) const;

};


#endif // ODEINTERFACE_HPP_INCLUDED
