
#include "cascadia.hpp"

namespace mfem
{

using namespace std;

WaveOperator::WaveOperator(FiniteElementSpace &fes_u_, FiniteElementSpace &fes_p_,
                           FiniteElementSpace &fes_m_, WaveTransferMap &wave_map_,
                           WaveObservationOp &wave_obs_,
                           Array<int> &ess_bdr_, Array<int> &nat_bdr_,
                           Array<int> &surfac_bdr_, Array<int> &bottom_bdr_, Array<int> &absorb_bdr_,
                           bool lump_, int height_, double c1_, double c2_, double c3_,
                           int param_rate_, int param_steps_, int obs_rate_, int obs_steps_,
                           bool adj_)
   : TimeDependentOperator(height_, 0.0),
     fes_u(fes_u_), fes_p(fes_p_), fes_m(fes_m_), wave_map(wave_map_), wave_obs(wave_obs_),
     ess_bdr(ess_bdr_), nat_bdr(nat_bdr_),
     surfac_bdr(surfac_bdr_), bottom_bdr(bottom_bdr_), absorb_bdr(absorb_bdr_),
     block_offsets(3),
     M_u(nullptr), M_p(nullptr), BlockM(nullptr), mVarf_u(nullptr), mVarf_p(nullptr),
     K_u(nullptr), K_p(nullptr), K_pt(nullptr), BlockK(nullptr), kVarf_u(nullptr), kVarf_p(nullptr),
     K_imp(nullptr), kVarf_imp(nullptr),
     f_coeff(nullptr), g_coeff(nullptr), g_nat_coeff(nullptr),
     tdf_coeff(nullptr), m_coeff(nullptr), m_gf(nullptr), p_gf(nullptr), m_gf_coeff(nullptr),
     fform(nullptr), gform(nullptr), rhs(nullptr),
     BlockM_prec(nullptr), invM_u(nullptr), invM_p(nullptr), T_oper(nullptr),
     z(height_), c1(c1_), c2(c2_), c3(c3_),
     param_rate(param_rate_), param_steps(param_steps_),
     obs_rate(obs_rate_), obs_steps(obs_steps_),
     adj(adj_),
     parameters(nullptr), data(nullptr)
{
   StopWatch chrono;
   chrono.Start();
   
   store_load = false;
   
   if (WaveSolution::IsKnown() && adj)
   {
      MFEM_ABORT("WaveOperator adjoint only implemented for unknown solution.");
   }

   // 1a. Set block offsets
   block_offsets[0] = 0;
   block_offsets[1] = fes_u.GetVSize();
   block_offsets[2] = fes_p.GetVSize();
   block_offsets.PartialSum();
   
   // 1b. Convert marked boundary attributes (essential BCs)
   //     to a list of true boundary DOFs
   fes_u.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   MFEM_VERIFY(ess_tdof_list.Size() == 0,
               "WaveOperator: ess_tdof are not implemented here.");

   // 2. Assemble mass matrix block operator
   BlockM = new BlockOperator(block_offsets);
   
   mVarf_u = new BilinearForm(&fes_u);
   mVarf_p = new BilinearForm(&fes_p);

   if (lump_) { cout << "Mass lumping is enabled." << endl; }
   if (lump_) { mVarf_u->AddDomainIntegrator(new LumpedIntegrator(new VectorMassIntegrator())); }
   else       { mVarf_u->AddDomainIntegrator(new VectorMassIntegrator()); }
   mVarf_u->Assemble();
   mVarf_u->Finalize();
   M_u = &(mVarf_u->SpMat());
   
   if (lump_) { mVarf_p->AddDomainIntegrator(new LumpedIntegrator(new MassIntegrator())); }
   else       { mVarf_p->AddDomainIntegrator(new MassIntegrator()); }
   ConstantCoefficient c2c3(c2/c3);
   if (WaveSolution::IsUnknown())
   {
      // Add surface ODE-type BC for dp/dt
      if (lump_) { mVarf_p->AddBoundaryIntegrator(new LumpedIntegrator(new MassIntegrator(c2c3)), surfac_bdr); }
      else       { mVarf_p->AddBoundaryIntegrator(new MassIntegrator(c2c3), surfac_bdr); }
      // mVarf_p->AddBoundaryIntegrator(new BoundaryMassIntegrator(c2c3), surfac_bdr);
   }
   mVarf_p->Assemble();
   mVarf_p->Finalize();
   M_p = &(mVarf_p->SpMat());

   BlockM->SetDiagonalBlock(0, M_u);
   BlockM->SetDiagonalBlock(1, M_p);
   
   cout << "Timer: Assembly of Block M Matrix     : " << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();

   // 3. Create mass matrix block solver
   //    with diagonal preconditioner
   //      P = [ diag(M_u)  0         ]
   //          [ 0          diag(M_p) ]
   BlockM_prec = new BlockDiagonalPreconditioner(block_offsets);

   invM_u = new DSmoother(*M_u);
   invM_p = new DSmoother(*M_p);
   invM_u->iterative_mode = false;
   invM_p->iterative_mode = false;

   BlockM_prec->SetDiagonalBlock(0, invM_u);
   BlockM_prec->SetDiagonalBlock(1, invM_p);

   const double rel_tol = 1e-8;
   const double abs_tol = 1e-14;
   const int max_iter = 30;

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(abs_tol);
   M_solver.SetMaxIter(max_iter);
   M_solver.SetPrintLevel(0);
   M_solver.SetPreconditioner(*BlockM_prec);
   M_solver.SetOperator(*BlockM);

   cout << "Timer: Solver setup for Block M Matrix: " << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();

   // 4. Assemble stiffness matrix block operator */
   double cc1,cc2;
   if (adj)
   {
      cc1 =  c2;
      cc2 = -c1;
   }
   else
   {
      cc1 =  c1;
      cc2 = -c2;
   }
   ConstantCoefficient c1_coeff(cc1);
   ConstantCoefficient c2_coeff(cc2);
   
   BlockK = new BlockOperator(block_offsets);
   
   kVarf_u = new MixedBilinearForm(&fes_p, &fes_u);
   kVarf_p = new MixedBilinearForm(&fes_p, &fes_u);

   kVarf_u->AddDomainIntegrator(new GradientIntegrator(c1_coeff));
   kVarf_u->Assemble();
   kVarf_u->Finalize();
   K_u = &(kVarf_u->SpMat());

   kVarf_p->AddDomainIntegrator(new GradientIntegrator(c2_coeff));
   kVarf_p->Assemble();
   kVarf_p->Finalize();
   K_p = &(kVarf_p->SpMat());

   K_pt = new TransposeOperator(K_p);

   BlockK->SetBlock(0,1, K_u);
   BlockK->SetBlock(1,0, K_pt);
   
   // Boundary integral for impedance BC
   if (WaveSolution::IsUnknown())
   {
      double a = (adj) ? -1.0 : 1.0; // sign change in adjoint absorbing BC
      ConstantCoefficient c1c2(a*sqrt(c1*c2)); // non-dim. wavespeed
      kVarf_imp = new BilinearForm(&fes_p);
      kVarf_imp->AddBoundaryIntegrator(new MassIntegrator(c1c2), absorb_bdr);
      kVarf_imp->Assemble();
      kVarf_imp->Finalize();
      K_imp = &(kVarf_imp->SpMat());
      BlockK->SetBlock(1,1, K_imp);
   }
   
   cout << "Timer: Assembly of Block K Matrix     : " << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();
   
   // 5. If using an implicit solver, then we also need to set up solver for (M+dt*K)
   T_oper = new SumOperator(BlockM, BlockK, 1.0);
   
   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(abs_tol);
   T_solver.SetMaxIter(max_iter);
   T_solver.SetPrintLevel(1);
   T_solver.SetPreconditioner(*BlockM_prec);
   T_solver.SetOperator(*T_oper);
   
   cout << "Timer: Solver setup for Block K Matrix: " << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();

   // 6. Define time-dependent function coefficients for RHS
   //    and initialize the RHS block load vector
   rhs = new BlockVector(block_offsets); // TODO: add mem type
   *rhs = 0.0;
   
   if (WaveSolution::IsKnown())
   {
      f_coeff = new VectorFunctionCoefficient(3, WaveSolution::fLoad);
      g_coeff = new FunctionCoefficient(WaveSolution::gLoad);
      g_nat_coeff = new VectorFunctionCoefficient(3, WaveSolution::gNatural);
      
      fform = new LinearForm();
      fform->Update(&fes_u, rhs->GetBlock(0), 0);
      fform->AddDomainIntegrator(new VectorDomainLFIntegrator(*f_coeff));
   
      gform = new LinearForm();
      gform->Update(&fes_p, rhs->GetBlock(1), 0);
      gform->AddDomainIntegrator(new DomainLFIntegrator(*g_coeff));
      // If using natural BC, add the boundary loads
      // note, this integrator uses the normal component of the given vector coefficient
      gform->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(*g_nat_coeff), nat_bdr);
      
      init_load = true;
   }
   else
   {
      // Load must be initialized via "TDFunctionLoad" or "StoreLoad"
      init_load = false;
   }
   
   cout << "Timer: Setup for Block RHS Vector     : " << chrono.RealTime() << " seconds." << endl << endl;
}

void WaveOperator::Mult(const Vector &up, Vector &dup_dt) const
{
   // Compute:
   //    Eqn: M*(d(up)_dt) + K*(up) = f
   //
   //    d(up)_dt = M^{-1}*(-K*(up) + f)
   MFEM_VERIFY(init_load,
               "Load must be initialized before using WaveOperator::Mult.");
   
   // Update time-dependent load vector (RHS)
//   cout << "updating load; time = " << this->GetTime() << endl;
   if (!store_load)
   {
      this->UpdateLoad();
   }
   
   // Multiply block stiffness matrix
   BlockK->Mult(up, z);
   z.Neg(); // z = -z

   // Mult operator adds load before inverting mass matrix
   z += *rhs;

   // Solve block mass matrix system
   M_solver.Mult(z, dup_dt);
   
   if (WaveSolution::IsStationary())
   {
      // dup_dt should be near-zero
      cout << "Stationary solution: ||dup_dt||_inf = " << dup_dt.Normlinf() << endl;
   }
}

// NOTE: Implicit time integration is not supported with partial assembly
void WaveOperator::ImplicitSolve(const double dt,
                                 const Vector &up, Vector &dup_dt)
{
   // Solve the equation:
   //    Eqn: M*(d(up)_dt) + K*(up+dt*d(up)_dt) = f(t+dt)
   //         (M+dt*K)*(d(up)_dt) + K*(up) = f(t+dt)
   //
   //    d(up)_dt = (M+dt*K)^{-1}*[-K*(up) + f(t+dt)]
   
   // Update time-dependent load vector (RHS)
   // note: implicit ODE solver->step has set time to t+dt before calling ImplicitSolve
//   cout << "updating load; time = " << this->GetTime() << endl;
   if (!store_load)
   {
      this->UpdateLoad();
   }
      
   // Multiply block stiffness matrix
   BlockK->Mult(up, z);
   z.Neg(); // z = -z

   // Mult operator adds load before inverting mass matrix
   z += *rhs;

   // Solve block stiffness matrix system M+dt*K
   T_oper->SetFactor(dt);
   T_solver.Mult(z, dup_dt);
   
   if (WaveSolution::IsStationary())
   {
      // dup_dt should be near-zero
      cout << "Stationary solution: ||d(up)/dt||_inf = " << dup_dt.Normlinf() << endl;
   }
}

void WaveOperator::TDFunctionLoad(function<double(const Vector &, double)> TDF)
{
   MFEM_VERIFY(WaveSolution::IsUnknown(),
               "WaveOperator::TDFunctionLoad : use for unknown solution only.");
   
   delete gform;
   gform = new LinearForm();
   gform->Update(&fes_p, rhs->GetBlock(1), 0);
   
   // TD FunctionCoefficient m - works only if m given as a function of (x,t)
   delete m_coeff; delete tdf_coeff;
   tdf_coeff = new FunctionCoefficient(TDF);
   m_coeff = new ProductCoefficient(c2, *tdf_coeff);
   gform->AddBoundaryIntegrator(new BoundaryLFIntegrator(*m_coeff), bottom_bdr);
   
   init_load = true;
   store_load = false;
}

void WaveOperator::StoreLoad(const GridFunction *const parameters_[])
{
   MFEM_VERIFY(WaveSolution::IsUnknown(), "Load term for unknown solution only.");
   MFEM_VERIFY(!adj, "Use StoreAdjLoad for adjoint operator.");
   MFEM_VERIFY(!init_load, "ResetLoad first.");
   
   delete gform;
   gform = new LinearForm();
   gform->Update(&fes_p, rhs->GetBlock(1), 0);
   
   // GridFunctionCoefficient m - works for any input GridFunction m
   // Note: We define the GridFunctionCoefficient based on state GF (fes_p).
   //       GFCoefficient MUST be defined using GF on fes_p for gform->AddBoundaryIntegrator.
   //       In update_load, the state GF must be updated for each new time step by mapping
   //       the parameter GF (fes_m) to state GF (fes_p) using SubMesh-to-Mesh transfer

   p_gf = new GridFunction(&fes_p);
   m_gf_coeff = new GridFunctionCoefficient(p_gf);
   gform->AddBoundaryIntegrator(new BoundaryLFIntegrator(*m_gf_coeff), bottom_bdr);
   
   parameters = new GridFunction *[param_steps];
   for (int k = 0; k < param_steps; k++)
   {
      parameters[k] = (GridFunction *)parameters_[k];
   }
   
   init_load = true;
   store_load = true;
   cycle_load = 0;
}

void WaveOperator::StoreAdjLoad(const Vector* const data_[])
{
   MFEM_VERIFY(WaveSolution::IsUnknown(), "Load term for unknown solution only.");
   MFEM_VERIFY(adj, "Not an adjoint operator.");
   MFEM_VERIFY(!init_load, "ResetLoad first.");
   
   delete gform;
   gform = new LinearForm();
   gform->Update(&fes_p, rhs->GetBlock(1), 0);
   
   p_gf = new GridFunction(&fes_p);
   m_gf_coeff = new GridFunctionCoefficient(p_gf);
   gform->AddBoundaryIntegrator(new BoundaryLFIntegrator(*m_gf_coeff), bottom_bdr);
   
   // Note: data is assumed to be given forward-in-time, i.e.
   //       d[0] <-> t_0, d[1] <-> t_1, ...
   data = new Vector *[obs_steps];
   for (int k = 0; k < obs_steps; k++)
   {
      data[k] = (Vector *)data_[k];
   }
   
   init_load = true;
   store_load = true;
   cycle_load = 0;
}

void WaveOperator::ResetLoad()
{
   if (parameters) { delete parameters; parameters = nullptr; }
   if (data) { delete data; data = nullptr; }
   init_load = false;
   store_load = false;
}

void WaveOperator::Cycle() const
{
   if (store_load)
   {
      if (adj) { UpdateAdjLoad(); }
      else     { UpdateLoad(); }
      cycle_load++;
   }
}

void WaveOperator::UpdateLoad() const
{
   MFEM_VERIFY(!adj, "Use UpdateAdjLoad for adjoint operator.");
   
   double t = this->GetTime();
   
   // Load for manufactured solution
   if (WaveSolution::IsKnown())
   {
      // Update time-dependent function coefficients
      f_coeff->SetTime(t);
      g_coeff->SetTime(t);
      g_nat_coeff->SetTime(t);

      // Assemble load for the first equation
      fform->Assemble();
      fform->SyncAliasMemory(*rhs);

      // Assemble load for the second equation
      gform->Assemble();
      gform->SyncAliasMemory(*rhs);
   }
   
   // Boundary load for unknown solution
   else
   {
      if (store_load)
      {
         // When using GridFunctionCoefficient m_gf_coeff,
         // we must update the underlying GridFunction p_gf.
         // Problem: UpdateLoad may be called multiple times per time step and
         //          expects intermediate values in multi-stage ODE solvers.
         if (cycle_load % param_rate == 0)
         {
            int k = cycle_load / param_rate;
            wave_map.ParamToState(*(parameters[k]), *p_gf);
            (*p_gf) *= c2;
         }
         else
         {
            // TODO: could instead linearly interpolate here
            *p_gf = 0;
         }
      }
      else
      {
         // Update time-dependent function coefficient
         m_coeff->SetTime(t);
      }
      // Load for the second equation (natural BC)
      // depends on parameter field m
      gform->Assemble();
      gform->SyncAliasMemory(*rhs);
   }
}

void WaveOperator::UpdateAdjLoad() const
{
   MFEM_VERIFY(WaveSolution::IsUnknown(), "Adj load for unknown solution only.");
   MFEM_VERIFY(adj, "Not an adjoint operator.");
   MFEM_VERIFY(store_load, "Adjoint must be defined by data (StoreAdjLoad).");
   
//   double t = this->GetTime();
   
   // When using GridFunctionCoefficient m_gf_coeff,
   // we must update the underlying GridFunction p_gf.
   // Problem: UpdateAdjLoad may be called multiple times per time step and
   //          expects intermediate values in multi-stage ODE solvers.
   if (cycle_load % obs_rate == 0)
   {
      // adjoint is computed backwards-in-time,
      // but we assume data is given with ordering forward-in-time
      int k = cycle_load / obs_rate;
      int l = (obs_steps-1) - k;
      wave_obs.ObservationToState(*(data[l]), *p_gf);
   }
   else
   {
      // TODO: could instead linearly interpolate here
      *p_gf = 0;
   }
   // Load for the second equation (natural BC)
   // depends on data field d
   gform->Assemble();
   gform->SyncAliasMemory(*rhs);
}

WaveOperator::~WaveOperator()
{
   delete invM_u;
   delete invM_p;
   delete BlockM_prec;
   
   delete mVarf_u;
   delete mVarf_p;
   delete BlockM;
   
   delete kVarf_u;
   delete kVarf_p;
   delete kVarf_imp;
   delete K_pt;
   delete BlockK;
   
   delete T_oper;
   
   delete f_coeff;
   delete g_coeff;
   delete g_nat_coeff;
   delete tdf_coeff;
   delete m_coeff;
   delete m_gf;
   delete p_gf;
   delete m_gf_coeff;
   
   delete fform;
   delete gform;
   
   delete rhs;
}

}
