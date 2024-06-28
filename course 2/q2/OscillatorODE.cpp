#include <iostream>
#include <cmath>
#include "OscillatorODE.hpp"

OscillatorODE::OscillatorODE(){
 ma=0.0;
}

void OscillatorODE:: ComputeF( const double t, const double x, double& f ) const
{
    f=-pow(ma,2.0)*x;

}

void OscillatorODE::ComputeAnalyticSolution( const double t, double &x ) const
{
    x=sin(ma*t);
}

void OscillatorODE:: seta(double a){
 ma=a;
}
