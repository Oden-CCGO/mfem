#ifndef CASCADIA_WAVE_SOL
#define CASCADIA_WAVE_SOL

#include "mfem.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/// Define global parameters
namespace Cascadia
{
   // Selected dimensional scales
   const double t0 = 1.0;   // time [s]
   const double l0 = 1.0e3; // length [m]
   const double p0 = 1.0e9; // pressure [Pa]
   
   // Derived dimensional scales
   const double u0 = l0 / t0; // velocity [m/s]
   
   // Physical constants
   const double rho = 1025;    // seawater density [kg/m^3]
   const double grv = 9.80665; // gravity [m/s^2]
   const double K   = 2.34e9;  // Bulk modulus [Pa]
   
   // Non-dimensional scalar coefficients
   const double c1_const = (p0*t0*t0) / (rho*l0*l0);
   const double c2_const = K / p0;
   const double c3_const = (rho*grv*l0) / p0;
   
   // Acoustic impedance variables for absorbing BC
   // note: non-dim. wave speed = c_acoustic / u0 = sqrt(c1*c2)
   const double c_acoustic = sqrt(K / rho);    // dim. wave speed [m/s]
   const double Z_acoustic = rho * c_acoustic; // dim. imp. const. [Pa s/m]
//   const double Z0 = (u0/p0) * Z_acoustic;     // non-dim. imp. const.
   
   extern void PrintScales();
   extern void PrintParams();
}

/// Define various manufactured solutions
namespace WaveSolution
{
   extern int problem;       // type of solution or forcing
   extern bool manufactured; // known/unknown solution
   extern bool stationary;   // stationary/time-dependent solution
   
   extern double t_final; // total simulation time
   extern int n_steps;    // number of time steps
   extern double dt;      // size of each time step
   
   extern double c1,c2,c3; // non-dimensional scalar coefficients
   
   extern double xmin,xmax,ymin,ymax,zmin,zmax; // mesh dimensions
   extern int nx,ny,nz; // number of elements in uniform grid
   
   // Functions returning initial values at a point x
   extern void uInitial(const mfem::Vector &x, mfem::Vector &u);
   extern double pInitial(const mfem::Vector &x);

   // Define the analytical solution and forcing terms / boundary conditions
   extern void uExact(const mfem::Vector &x, double t, mfem::Vector &u);
   extern double pExact(const mfem::Vector &x, double t);
   extern double uDivExact(const mfem::Vector &x, double t);
   extern void pGradExact(const mfem::Vector &x, double t, mfem::Vector &p_grad);
   extern void uDtExact(const mfem::Vector &x, double t, mfem::Vector &u_dt);
   extern double pDtExact(const mfem::Vector &x, double t);
   extern void fLoad(const mfem::Vector &x, double t, mfem::Vector &f);
   extern double gLoad(const mfem::Vector &x, double t);
   extern double fNatural(const mfem::Vector &x, double t);
   extern void gNatural(const mfem::Vector &x, double t, mfem::Vector &g);
   extern double mParameter(const mfem::Vector &x, double t);
   
   extern void Init(int problem_); // Must be called to initialize
   
   extern bool IsKnown();
   extern bool IsUnknown();
   
   extern bool IsStationary();
   extern bool IsTimeDependent();
   
   extern void PrintInfo();
   
   /// Define a linear function for testing purposes
   extern double LinearFunction(const mfem::Vector &x);
};

}

#endif
