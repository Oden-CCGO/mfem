#ifndef CASCADIA_WAVE_OP
#define CASCADIA_WAVE_OP

#include "mfem.hpp"

#include "../common/sum_operator.hpp"

#include "wave_sol.hpp"
#include "wave_map.hpp"
#include "wave_obs.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/** After spatial discretization, the first-order mixed acoustic wave model can be written as:
 *
 *     explicit:    d(up)/dt = M^{-1} (-K(up) + f(t))
 *
 *     implicit:    d(up)/dt = (M+dt K)^{-1} (-K(up) + f(t+dt))
 *
 *  where "up" is the block vector representing (velocity, pressure), M is the block mass matrix,
 *  and K is the first-order mixed wave operator with (non-dimensional) material parameters:
 *  c1, c2, c3; t is time, dt is the time increment, and f is the load vector.
 *
 *  WaveOperator represents the right-hand side of the above ODE.
 */
class WaveOperator : public TimeDependentOperator
{
protected:
   FiniteElementSpace &fes_u;   // velocity FE space
   FiniteElementSpace &fes_p;   // pressure FE space
   FiniteElementSpace &fes_m;   // parameter FE space
   WaveTransferMap &wave_map;   // defines mapping between FE spaces (param<->state)
   WaveObservationOp &wave_obs; // defines mapping between FE spaces (data <->state)
   
   Array<int> ess_bdr;       // this list marks boundary attributes (essential BCs)
   Array<int> nat_bdr;       // this list marks boundary attributes (natural BCs)
   Array<int> ess_tdof_list; // this list remains empty for natural BC
   
   Array<int> surfac_bdr; // sea surface boundary markers
   Array<int> bottom_bdr; // sea bottom boundary markers
   Array<int> absorb_bdr; // boundary markers for absorbing BCs

   Array<int> block_offsets; // offsets for block operator

   SparseMatrix *M_u, *M_p; // velocity, pressure mass matrix
   BlockOperator *BlockM;   // block operator for mass matrices
   BilinearForm *mVarf_u;   // var form for first equation mass term
   BilinearForm *mVarf_p;   // var form for second equation mass term

   SparseMatrix *K_u, *K_p;    // velocity, pressure stiffness matrix
   TransposeOperator *K_pt;    // transpose of pressure stiffness matrix
   BlockOperator *BlockK;      // block operator for stiffness matrices
   MixedBilinearForm *kVarf_u; // var form for first equation stiffness term
   MixedBilinearForm *kVarf_p; // var form for second equation stiffness term
   
   SparseMatrix *K_imp;     // impedance term stiffness matrix
   BilinearForm *kVarf_imp; // var form for impedance boundary term
   
   VectorFunctionCoefficient *f_coeff;
   FunctionCoefficient *g_coeff;
   VectorFunctionCoefficient *g_nat_coeff;
   FunctionCoefficient *tdf_coeff;
   ProductCoefficient *m_coeff;
   GridFunction *m_gf;
   GridFunction *p_gf;
   GridFunctionCoefficient *m_gf_coeff;
   
   LinearForm *fform; // linear form for first equation RHS
   LinearForm *gform; // linear form for second equation RHS
   BlockVector *rhs;  // block vector for RHS

   CGSolver M_solver; // Krylov solver for block mass matrix M
   BlockDiagonalPreconditioner *BlockM_prec; // Preconditioner for block mass matrix M
   Solver *invM_u, *invM_p; // Solvers for velocity and pressure mass matrices

   GMRESSolver T_solver; // Implicit solver for T
   SumOperator *T_oper;  // T := M + dt K

   mutable Vector z; // auxiliary vector
   
   double c1,c2,c3; // non-dimensional scalar coefficients
   bool init_load;  // indicates whether load has been set up
   bool store_load; // indicates whether load is stored by GridFunctionCoefficient
   mutable int cycle_load; // auxiliary value used for updating load with time step
   
   const int param_rate, param_steps; // parameter frequency
   const int obs_rate, obs_steps;     // observation frequency
   
   const bool adj; // indicates whether this is the adjoint operator
   
   GridFunction **parameters;
   Vector **data;

public:
   WaveOperator(FiniteElementSpace &fes_u_, FiniteElementSpace &fes_p_,
                FiniteElementSpace &fes_m_, WaveTransferMap &wave_map_,
                WaveObservationOp &wave_obs_,
                Array<int> &ess_bdr_, Array<int> &nat_bdr_,
                Array<int> &surfac_bdr_, Array<int> &bottom_bdr_, Array<int> &absorb_bdr_,
                bool lump_, int height_, double c1_, double c2_, double c3_,
                int param_rate_, int param_steps_, int obs_rate_, int obs_steps_,
                bool adj_=false);

   virtual void Mult(const Vector &up, Vector &dup_dt) const;

   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);
   
   /// Define space-time RHS directly via vector of GridFunctions
   void StoreLoad(const GridFunction *const parameters_[]);
   
   /// Define space-time RHS via a time-dependent function
   void TDFunctionLoad(std::function<double(const Vector &, double)> TDF);
   
   /// Advance cycle for stored loads
   void Cycle() const;
   
   /// Update time-dependent load (rhs)
   void UpdateLoad() const;

   /// Define space-time RHS for the adjoint directly via vector of GridFunctions
   void StoreAdjLoad(const Vector* const data_[]);
   
   /// Update time-dependent load for the adjoint (rhs)
   void UpdateAdjLoad() const;
   
   /// Reset internal variables for stored load
   void ResetLoad();

   int ParamRate() { return param_rate; }
   int ParamSteps() { return param_steps; }
   int ObsRate() { return obs_rate; }
   int ObsSteps() { return obs_steps; }
   bool IsAdjoint() { return adj; }

   FiniteElementSpace *VelocityFESpace() { return &fes_u; }
   const FiniteElementSpace *VelocityFESpace() const { return &fes_u; }
   
   FiniteElementSpace *PressureFESpace() { return &fes_p; }
   const FiniteElementSpace *PressureFESpace() const { return &fes_p; }
   
   FiniteElementSpace *ParamFESpace() { return &fes_m; }
   const FiniteElementSpace *ParamFESpace() const { return &fes_m; }

   virtual ~WaveOperator();
};

}

#endif
