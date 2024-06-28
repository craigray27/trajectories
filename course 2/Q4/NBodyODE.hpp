#ifndef NBODYODE_HPP_INCLUDED
#define NBODYODE_HPP_INCLUDED

#include "ODEInterface.hpp"
#include "Vector.hpp"


class NBodyODE: public ODEInterface
{
  public:
      NBodyODE();
      void ComputeF( const double t, const Vector& x,  Vector& f, const int j ) const;

      void setm(Vector& NBodymass,Vector** &Nlocation,double number);//set Vector of mass and matrix of initial location, number=N

      void setlocation(Vector& p,const int j);//change the location of body j

      Vector** setmatrix(int N);

  private:
    double N;
    Vector** location;//initial location of N bodies
    Vector* mass;//Vector of mass
    const double G=6.674e-11;

};

#endif // NBODYODE_HPP_INCLUDED
