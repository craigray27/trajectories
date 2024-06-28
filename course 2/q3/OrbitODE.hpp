#ifndef OSCILLATORODE_HPP_INCLUDED
#define OSCILLATORODE_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "Vector.hpp"


class OrbitODE: public ODEInterface
{
  public:
      OrbitODE();
      void ComputeF( const double t, const Vector& x,  Vector& f ) const;
      void ComputeAnalyticSolution( const double t,  Vector& x ) const;

      void setm(double a,Vector& p);//set mass and xp

  private:
    double mass;//mass of earth
    Vector* xp;//location of earth
    const double G=6.674e-11;

};


#endif // OSCILLATORODE_HPP_INCLUDED


