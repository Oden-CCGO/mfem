
#include "cascadia.hpp"

namespace mfem
{

using namespace std;

namespace Cascadia
{

void PrintScales()
{
   std::cout << std::endl << "Cascadia::PrintScales" << std::endl;
   
   std::cout << "- Dimensional scales:" << std::endl;
   std::cout << "   t0 = " << t0 << " s" << std::endl;
   std::cout << "   l0 = " << l0 << " m" << std::endl;
   std::cout << "   p0 = " << p0 << " Pa" << std::endl;
   std::cout << "   u0 = " << u0 << " m/s" << std::endl;
}

void PrintParams()
{
   std::cout << std::endl << "Cascadia::PrintParams" << std::endl;
   
   std::cout << "- Physical constants:" << std::endl;
   std::cout << "   Seawater density (rho) = " << rho << " kg/m^3" << std::endl;
   std::cout << "   Gravity          (grv) = " << grv << " m/s^2" << std::endl;
   std::cout << "   Bulk modulus     (K)   = " << K   << " Pa" << std::endl << std::endl;
   
   std::cout << "- Impedance variables:" << std::endl;
   std::cout << "   Acoustic wave speed     (c_acoustic) = " << c_acoustic << " m/s" << std::endl;
   std::cout << "   Acoustic wave impedance (Z_acoustic) = " << Z_acoustic << " Pa s/m" << std::endl;
}

} // close namespace Cascadia

namespace WaveSolution
{

/// declarations
int problem, n_steps;
bool manufactured, stationary;
double t_final, dt, c1, c2, c3;
double xmin,xmax,ymin,ymax,zmin,zmax;

/// Implementation of WaveSolution:
/// - initial conditions
/// - analytical solutions
/// - boundary conditions
/// - forcing terms

/// Initialization of WaveSolution parameters
void Init(int problem_)
{
   problem = problem_;
   switch(problem)
   {
      case 1: case 2: case 3: // stationary, known solution
         manufactured = true;
         stationary = true;
         c1 = 1.0; c2 = 1.0; c3 = 1.0;
         break;
      case 10: case 20: case 30: case 40: case 50: // time-dependent, known solution
         manufactured = true;
         stationary = false;
         c1 = 1.0; c2 = 1.0; c3 = 1.0;
         break;
      case 100: // unknown solution (forcing only)
         manufactured = false;
         stationary = false;
         c1 = Cascadia::c1_const;
         c2 = Cascadia::c2_const;
         c3 = Cascadia::c3_const;
         Cascadia::PrintScales();
         Cascadia::PrintParams();
         break;
      default:
         MFEM_ABORT("init : undefined problem");
   }
   PrintInfo();
}

/// Initial condition (velocity u)
void uInitial(const Vector &x, Vector &u)
{
   if (manufactured)
   {
      uExact(x, 0.0, u);
   }
   else
   {
      u(0) = 0.0;
      u(1) = 0.0;
      u(2) = 0.0;
   }
}

/// Initial condition (pressure p)
double pInitial(const Vector &x)
{
   if (manufactured)
   {
      return pExact(x, 0.0);
   }
   else
   {
      return 0.0;
   }
}

/// Manufactured exact solution (velocity u)
void uExact(const Vector &x, double t, Vector &u)
{
   MFEM_VERIFY(manufactured,
               "uExact only defined for manufactured solution.");

   double xi(x(0));
   double yi(x(1));
   double zi(x(2));

   switch(problem)
   {
      case 1:
         u(0) = - sin(M_PI*xi);
         u(1) = - sin(M_PI*yi);
         u(2) = - sin(M_PI*zi);
         break;
      case 2:
         u(0) = xi;
         u(1) = yi;
         u(2) = zi;
         break;
      case 3:
         u(0) = xi*xi;
         u(1) = yi*yi;
         u(2) = zi*zi;
         break;
      case 10:
         u(0) = - sin(M_PI*xi)*exp(-t);
         u(1) = - sin(M_PI*yi)*exp(-t);
         u(2) = - sin(M_PI*zi)*exp(-t);
         break;
      case 20:
         u(0) = xi*exp(-t);
         u(1) = yi*exp(-t);
         u(2) = zi*exp(-t);
         break;
      case 30:
         u(0) = xi*xi*exp(-t);
         u(1) = yi*yi*exp(-t);
         u(2) = zi*zi*exp(-t);
         break;
      case 40:
         u(0) = xi*xi*(t+1);
         u(1) = yi*yi*(t+1);
         u(2) = zi*zi*(t+1);
         break;
      case 50:
         u(0) = xi*xi*(t*t+1);
         u(1) = yi*yi*(t*t+1);
         u(2) = zi*zi*(t*t+1);
         break;
      default:
         MFEM_ABORT("uExact : undefined problem");
   }
}

/// Manufactured exact solution (pressure p)
double pExact(const Vector &x, double t)
{
   MFEM_VERIFY(manufactured,
               "pExact only defined for manufactured solution.");

   double xi(x(0));
   double yi(x(1));
   double zi(x(2));
   
   double p;
   switch(problem)
   {
      case 1:
         p = sin(M_PI*xi)*sin(M_PI*yi)*sin(M_PI*zi);
         break;
      case 2:
         p = xi + yi + zi;
         break;
      case 3:
         p = xi*xi + yi*yi + zi*zi;
         break;
      case 10:
         p = sin(M_PI*xi)*sin(M_PI*yi)*sin(M_PI*zi)*exp(-t);
         break;
      case 20:
         p = (xi + yi + zi)*exp(-t);
         break;
      case 30:
         p = (xi*xi + yi*yi + zi*zi)*exp(-t);
         break;
      case 40:
         p = (xi*xi + yi*yi + zi*zi)*(t+1);
         break;
      case 50:
         p = (xi*xi + yi*yi + zi*zi)*(t*t+1);
         break;
      default:
         p = 0;
         MFEM_ABORT("pExact : undefined problem");
   }
   return p;
}

/// Manufactured exact solution (divergence of velocity u)
double uDivExact(const Vector &x, double t)
{
   double xi(x(0));
   double yi(x(1));
   double zi(x(2));
   
   double u_div;
   switch(problem)
   {
      case 1:
         u_div = -M_PI*(cos(M_PI*xi)+cos(M_PI*yi)+cos(M_PI*zi));
         break;
      case 2:
         u_div = 3;
         break;
      case 3:
         u_div = 2*(xi + yi + zi);
         break;
      case 10:
         u_div = -M_PI*(cos(M_PI*xi)+cos(M_PI*yi)+cos(M_PI*zi))*exp(-t);
         break;
      case 20:
         u_div = 3*exp(-t);
         break;
      case 30:
         u_div = 2*(xi + yi + zi)*exp(-t);
         break;
      case 40:
         u_div = 2*(xi + yi + zi)*(t+1);
         break;
      case 50:
         u_div = 2*(xi + yi + zi)*(t*t+1);
         break;
      default:
         u_div = 0;
         MFEM_ABORT("uDivExact : undefined problem");
   }
   return u_div;
}

/// Manufactured exact solution (gradient of pressure p)
void pGradExact(const Vector &x, double t, Vector &p_grad)
{
   double xi(x(0));
   double yi(x(1));
   double zi(x(2));
   
   switch(problem)
   {
      case 1:
         p_grad(0) = M_PI*(cos(M_PI*xi)*sin(M_PI*yi)*sin(M_PI*zi));
         p_grad(1) = M_PI*(sin(M_PI*xi)*cos(M_PI*yi)*sin(M_PI*zi));
         p_grad(2) = M_PI*(sin(M_PI*xi)*sin(M_PI*yi)*cos(M_PI*zi));
         break;
      case 2:
         p_grad = 1;
         break;
      case 3:
         p_grad(0) = 2*xi;
         p_grad(1) = 2*yi;
         p_grad(2) = 2*zi;
         break;
      case 10:
         p_grad(0) = M_PI*(cos(M_PI*xi)*sin(M_PI*yi)*sin(M_PI*zi))*exp(-t);
         p_grad(1) = M_PI*(sin(M_PI*xi)*cos(M_PI*yi)*sin(M_PI*zi))*exp(-t);
         p_grad(2) = M_PI*(sin(M_PI*xi)*sin(M_PI*yi)*cos(M_PI*zi))*exp(-t);
         break;
      case 20:
         p_grad = exp(-t);
         break;
      case 30:
         p_grad(0) = 2*xi*exp(-t);
         p_grad(1) = 2*yi*exp(-t);
         p_grad(2) = 2*zi*exp(-t);
         break;
      case 40:
         p_grad(0) = 2*xi*(t+1);
         p_grad(1) = 2*yi*(t+1);
         p_grad(2) = 2*zi*(t+1);
         break;
      case 50:
         p_grad(0) = 2*xi*(t*t+1);
         p_grad(1) = 2*yi*(t*t+1);
         p_grad(2) = 2*zi*(t*t+1);
         break;
      default:
         MFEM_ABORT("pGradExact : undefined problem");
   }
}

/// Manufactured exact solution (time derivative of velocity u)
void uDtExact(const Vector &x, double t, Vector &u_dt)
{
   switch(problem)
   {
      case 1: case 2: case 3:
         u_dt = 0.0;
         break;
      case 10: case 20: case 30:
         uExact(x,t,u_dt);
         u_dt.Neg();
         break;
      case 40:
         uExact(x,t,u_dt);
         u_dt /= (t+1);
         break;
      case 50:
         uExact(x,t,u_dt);
         u_dt *= (2*t)/(t*t+1);
         break;
      default:
         MFEM_ABORT("uDtExact : undefined problem");
   }
}

/// Manufactured exact solution (time derivative of pressure p)
double pDtExact(const Vector &x, double t)
{
   double p_dt;
   switch(problem)
   {
      case 1: case 2: case 3:
         p_dt = 0.0;
         break;
      case 10: case 20: case 30:
         p_dt = -pExact(x,t);
         break;
      case 40:
         p_dt = pExact(x,t)/(t+1);
         break;
      case 50:
         p_dt = pExact(x,t)*(2*t)/(t*t+1);
         break;
      default:
         p_dt = 0;
         MFEM_ABORT("pDtExact : undefined problem");
   }
   return p_dt;
}

/// Manufactured vector-load f (first equation)
void fLoad(const Vector &x, double t, Vector &f)
{
   Vector u_dt(3);
   switch(problem)
   {
      case 1: case 2: case 3:
         pGradExact(x,t,f);
         f *= c1;
         break;
      case 10: case 20: case 30: case 40: case 50:
         pGradExact(x,t,f);
         f *= c1;
         uDtExact(x,t,u_dt);
         f += u_dt;
         break;
      default:
         MFEM_ABORT("fLoad : undefined problem");
   }
}

/// Manufactured scalar-load g (second equation)
double gLoad(const Vector &x, double t)
{
   double g;
   switch(problem)
   {
      case 1: case 2: case 3:
         g = uDivExact(x,t);
         g *= c2;
         break;
      case 10: case 20: case 30: case 40: case 50:
         g = uDivExact(x,t);
         g *= c2;
         g += pDtExact(x,t);
         break;
      default:
         g = 0.0;
         MFEM_ABORT("gLoad : undefined problem");
   }
   return g;
}

/// Load for natural BC ( (-c1*p) , \tau \dot n )
double fNatural(const Vector &x, double t)
{
   return -c1*pExact(x,t);
}

/// Load for natural BC ( (-c2*u) \dot n , v )
void gNatural(const Vector &x, double t, Vector &g)
{
   // Manufactured solution boundary load
   if (IsKnown())
   {
      uExact(x,t,g);
      g *= -c2;
   }
   else
   {
      MFEM_ABORT("gNatural : solving for unknown solution.")
   }
}

/// Load for natural BC: c2*(m,v)
double mParameter(const Vector &x, double t)
{
   double m_load;
   // Unknown solution: boundary load given by parameter field m
   if (IsUnknown())
   {
      double xi(x(0));
      double yi(x(1));
   
      // source term not non-dimensionalized in time yet (assumes t0 = 1 s)
      double ar = 0.1; // Rise amplitude [scale of m]
      double xr = 1.0; // Rise width in x [scale of km]
      double yr = 2.0; // Rise width in y [scale of km]
      double tr = 2.0; // Rise time [scale of s]
      double xi_c = (xmax-xmin)/2; // Rise center in x
      double yi_c = (ymax-ymin)/2; // Rise center in y
      
      if (t==0)
      {
         cout << endl << "mParameter:" << endl;
         cout << " ar = " << ar << " m" << endl;
         cout << " xr = " << xr << " km" << endl;
         cout << " yr = " << yr << " km" << endl;
         cout << " xi_c = " << xi_c << " km" << endl;
         cout << " yi_c = " << yi_c << " km" << endl;
         cout << " tr = " << tr << " s" << endl << endl;
      }
      
      // Non-dimensionalize length scales
      ar /= Cascadia::l0;
      xr *= (1000/Cascadia::l0);
      yr *= (1000/Cascadia::l0);
      
      // Note: using dimensionalized time; function assumes tr in seconds
      if (t <= tr)
      {
         m_load = ar * exp(-pow((xi-xi_c)/xr,2)-pow((yi-yi_c)/yr,2)) * M_PI/(2*tr) * sin(M_PI*Cascadia::t0*t/tr);
      }
      else
      {
         m_load = 0;
      }
      //       c2*m <-- c2*(db/dt) = c2*(-u_n)
      //       coefficient c2 is accounted for in assembling the load (UpdateLoad)
   }
   else
   {
      m_load = 0;
      MFEM_ABORT("mParameter : solving for known solution.")
   }
   return m_load;
}

bool IsKnown()   { return  manufactured; };
bool IsUnknown() { return !manufactured; };

bool IsStationary()    { return  stationary; };
bool IsTimeDependent() { return !stationary; };
   
void PrintInfo()
{
   std::cout << std::endl << "WaveSolution::PrintInfo" << std::endl;
   std::cout << "- Problem configuration:" << std::endl;
   std::cout << "   problem number : " << problem << std::endl;
   std::cout << "   is manufactured: " << manufactured << std::endl;
   std::cout << "   is stationary  : " << stationary << std::endl << std::endl;
   
   std::cout << "- Non-dimensional coefficients:" << std::endl;
   std::cout << "        c1 = " << c1 << std::endl;
   std::cout << "        c2 = " << c2 << std::endl;
   std::cout << "        c3 = " << c3 << std::endl << std::endl;
   
   std::cout << "- Time-stepping parameters:" << std::endl;
   std::cout << "   t_final = " << t_final*(Cascadia::t0) << " s" << std::endl;
   std::cout << "   n_steps = " << n_steps << std::endl;
   std::cout << "        dt = " << dt*(Cascadia::t0) << " s" << std::endl << std::endl;
   
   std::cout << "- Mesh dimensions (km):" << std::endl;
   std::cout << "   (xmin,xmax) = (" << xmin*(Cascadia::l0/1000)
                         << "," << xmax*(Cascadia::l0/1000) << ")" << std::endl;
   std::cout << "   (ymin,ymax) = (" << ymin*(Cascadia::l0/1000)
                         << "," << ymax*(Cascadia::l0/1000) << ")" << std::endl;
   std::cout << "   (zmin,zmax) = (" << zmin*(Cascadia::l0/1000)
                         << "," << zmax*(Cascadia::l0/1000) << ")" << std::endl;
}
   
double LinearFunction(const mfem::Vector &x)
{
   return x(0);
}

} // close namespace WaveSolution

} // close namespace mfem
