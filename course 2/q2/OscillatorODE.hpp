#ifndef OSCILLATORODE_HPP_INCLUDED
#define OSCILLATORODE_HPP_INCLUDED

#include "ODEInterface.hpp"

class OscillatorODE: public ODEInterface
{
  public:
      OscillatorODE();
      void ComputeF( const double t, const double x, double& f ) const;
      void ComputeAnalyticSolution( const double t, double &x ) const;
      void seta(double a);
  private:
    double ma;


};


#endif // OSCILLATORODE_HPP_INCLUDED
