#ifndef CASCADIA_WAVE_PRIOR
#define CASCADIA_WAVE_PRIOR

#include "mfem.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/** Action of the Hessian of the regularization operator:
 *
 *     Strong form:  alpha1  |m|^2 - alpha2 |dm/dx|^2 - alpha3 |dm/dt|^2
 *               dm/dx \cdot n = 0 on Gamma
 *               TODO: specify BC at t_k, k=1,N (need symmetry!) ; where N = param_steps
 *
 *     Weak form:  alpha1 (m_k, v) + alpha2 (grad m_k, grad v) + (alpha3/dt^2) [(-m_{k-1}, v) + 2(m_k, v) - (m_{k+1}, v)] , k = 2, ..., N-1
 *              TODO: ..., k = 1
 *              TODO: ..., k = N
 */
class WavePrior : public Operator
{
protected:
   FiniteElementSpace &fes_m;   // parameter FE space
   
   SparseMatrix *M;     // mass matrix (fe space)
   BilinearForm *m_var; // variational form (fe space)
   
   BlockMatrix *G;      // block mass matrix (space/time)
   Solver *M_prec;      // mass matrix preconditioner (fe space)
   BlockDiagonalPreconditioner *BlockM_prec; // block mass matrix preconditioner (space/time)
   CGSolver M_solver;   // mass matrix solver (space/time)
   
   SparseMatrix *K;     // stiffness matrix (fe space)
   BilinearForm *k_var; // variational form (fe space)
   
   SparseMatrix *C; // sparse matrix (auxiliary)
   SparseMatrix *L; // sparse matrix (auxiliary)
   
   SparseMatrix *D; // sparse matrix (auxiliary)
   BlockMatrix *B;  // sparse matrix (space/time)
   
   Array<int> block_offsets;
   
   // Regularization parameters
   ConstantCoefficient *alpha1; // |m|^2
   ConstantCoefficient *alpha2; // |dm/dx|^2
   ConstantCoefficient *alpha3; // |dm/dt|^2
   
   const int type; // type of prior (1: Laplacian; 2: Bi-Laplacian)
   
   const double time_step, dt;
   const int param_rate, param_steps; // parameter frequency; dt := time_step * param_rate
   
   mutable Vector z; // auxiliary vector
   
   /// Reindex csr matrix dof ordering from time->space to space->time (outer->inner)
   SparseMatrix* ReindexCSR(const SparseMatrix *R) const;
   /// Auxiliary routine for I/O
   void MatrixToFile(bool binary, bool mass);
public:
   WavePrior(FiniteElementSpace &fes_m_,
             int height_, int type_,
             double time_step_, int param_rate_, int param_steps_,
             double alpha1_=1.0, double alpha2_=1.0, double alpha3_=0.0);

   virtual void Mult(const Vector &x, Vector &y) const ;
   
   /// Apply mass matrix
   virtual void MultMass(const Vector &x, Vector &y) const ;
   /// Apply first part of prior/regularization only - alpha1
   virtual void MultReg1(const Vector &x, Vector &y) const ;
   /// Apply second part of prior/regularization only - alpha2
   virtual void MultReg2(const Vector &x, Vector &y) const ;
   /// Apply third part of prior/regularization only - alpha3
   virtual void MultReg3(const Vector &x, Vector &y) const ;
   
   /// Specifies whether PriorToFile re-indexes CSR matrix before
   /// writing to file, as is needed by the FFT matvec code
   static bool reindex;
   
   /// Write prior as a csr matrix to file
   void PriorToFile(bool binary=false);

   int ParamRate() { return param_rate; }
   int ParamSteps() { return param_steps; }
   
   FiniteElementSpace *ParamFESpace() { return &fes_m; }
   const FiniteElementSpace *ParamFESpace() const { return &fes_m; }

   virtual ~WavePrior();
};

}

#endif
