#include <iostream>
#include <cmath>
#include <iomanip>
#include "Vector.hpp"
#include "OrbitODE.hpp"
#include "ForwardEulerSolver.hpp"
#include "SymplecticEulerSolver.hpp"
#include "StoermerVerletSolver.hpp"

const double G=6.674e-11;
const double masse=5.972e24;//mass of earth

const double massm=7.342e22;//mass of moon

const double r0=3.844e8;//Distance between earth and moon

const double rearth=6.378e6;//radius of earth
const double rmoon=1.737e6;//radius of moon

const double pi=3.141592653589793;




