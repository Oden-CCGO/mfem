
#include "cascadia.hpp"
#include "hdf5.h"

#include <random>

namespace mfem
{

using namespace std;

WaveParamToObs::WaveParamToObs(WaveOperator *wave_fwd_, WaveOperator *wave_adj_,
                               WaveObservationOp &wave_obs_, WaveTransferMap &wave_map_,
                               WaveVis &wave_vis_,
                               ODESolver &ode_solver_, int n_steps_, double dt_)
   : wave_fwd(wave_fwd_), wave_adj(wave_adj_),
     wave_obs(wave_obs_), wave_map(wave_map_),
     wave_vis(wave_vis_),
     ode_solver(ode_solver_), n_steps(n_steps_), dt(dt_)
{
   WaveOperator *tmp = nullptr;
   if      (wave_fwd) { tmp = wave_fwd; }
   else if (wave_adj) { tmp = wave_adj; }
   
   MFEM_VERIFY(tmp, "Either fwd or adj operator must be defined.");
   
   n_param = tmp->ParamFESpace()->GetNDofs();
   dim = tmp->PressureFESpace()->GetMesh()->Dimension();
   
   param_rate  = tmp->ParamRate();
   param_steps = tmp->ParamSteps();
   obs_rate  = tmp->ObsRate();
   obs_steps = tmp->ObsSteps();
   
   n_obs = wave_obs.GetNrSensors();
   
   vis_fwd = 0;
   vis_adj = 0;
}

void WaveParamToObs::GetObs(Vector** &obs) const
{
   StopWatch chrono;
   chrono.Start();
   
   // 1. Define the BlockStructure of the problem, i.e. define the array of
   //    offsets for each variable. The last component of the Array is the sum
   //    of the dimensions of each block.
   Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = wave_fwd->VelocityFESpace()->GetVSize();
   block_offsets[2] = wave_fwd->PressureFESpace()->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "WaveParamToObs::Mult \n";
   std::cout << "dim(DG)    = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(CG)    = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(DG+CG) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n";
   
   // 2. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
   //MemoryType mt = device.GetMemoryType();
   BlockVector x(block_offsets); //, mt);
   
   // 3. Set the initial conditions for (u,p)
   
   // Initial condition (homogeneous)
   x = 0.0;
   
   // For verification purposes (when using manufactured solution),
   // the initial condition must satisfy exact manufactured solution
   GridFunction u_gf(wave_fwd->VelocityFESpace());
   GridFunction p_gf(wave_fwd->PressureFESpace());
   
   if (WaveSolution::IsKnown())
   {
      VectorFunctionCoefficient ucoeff(dim, WaveSolution::uInitial);
      FunctionCoefficient pcoeff(WaveSolution::pInitial);
   
      u_gf.ProjectCoefficient(ucoeff);
      p_gf.ProjectCoefficient(pcoeff);
      u_gf.GetTrueDofs(x.GetBlock(0)); // extracts true dofs from GF; inserts them into true dof vector u [x(0)]
      p_gf.GetTrueDofs(x.GetBlock(1)); // extracts true dofs from GF; inserts them into true dof vector p [x(1)]
   }
   
   // 4. Initiate visualization data in the ParaView format and
   //    write the initial solution to file
   ParaViewDataCollection *paraview_dc;
   if (wave_vis.IsVis())
   {
      ostringstream collection_oss;
      collection_oss << "WaveSolution_" << setfill('0') << setw(6) << vis_fwd++;
      string collection = collection_oss.str();
      wave_vis.NewCollection(collection, wave_fwd->PressureFESpace()->GetMesh());
      paraview_dc = wave_vis.ParaviewDC();
      u_gf.MakeRef(wave_fwd->VelocityFESpace(), x.GetBlock(0), 0);
      p_gf.MakeRef(wave_fwd->PressureFESpace(), x.GetBlock(1), 0);
      // paraview_dc->RegisterField("velocity",&u_gf);
      paraview_dc->RegisterField("pressure",&p_gf);
      paraview_dc->SetCycle(0);
      paraview_dc->SetTime(0.0);
      chrono.Clear();
      paraview_dc->Save();
      cout << endl << "p2o Mult: paraview I/O took " << chrono.RealTime() << " seconds." << endl;
   }
   
   // 5. Initialize ODE solver with forward wave operator
   ode_solver.Init(*wave_fwd);
   
   // 6. Allocate vectors for observations
   if (!obs)
   {
      obs = new Vector*[obs_steps];
      for (int k = 0; k < obs_steps; k++) { obs[k] = nullptr; }
   }
   for (int k = 0; k < obs_steps; k++)
   {
      if (obs[k])
      {
         obs[k]->SetSize(n_obs);
      }
      else
      {
         obs[k] = new Vector(n_obs);
      }
   }
   
   // 7. Perform time-stepping
   chrono.Clear();

   double t_curr = 0.0;
   for (int k = 1; k <= n_steps; k++)
   {
      // perform time step (increments t_curr)
      double dtc = dt;
      wave_fwd->Cycle(); // necessary for updating GridFunctionCoefficient (if store_load)
      ode_solver.Step(x, t_curr, dtc);

      if (abs(dtc-dt) > 1.0e-14)
      {
         MFEM_WARNING("abs(dtc-dt) = " << abs(dtc-dt));
      }
      
      double t_dim = t_curr * Cascadia::t0;
      if (k % 5 == 0)
      {
         cout << "=========================================" << endl
              << "time step = " << k << ", t = "
                       << setprecision(2) << t_dim << " s"   << endl
              << "=========================================" << endl;
         if (k % 10 == 0)
         {
            cout << endl << " Time-stepping until now took " << chrono.RealTime() << " seconds." << endl;
            cout <<         " Average time/step until now: " << chrono.RealTime()/k << " seconds." << endl;
         }
      }
      
      // Extract observations every "obs_rate"-th step
      if (k % obs_rate == 0)
      {
         int l = k / obs_rate;
         p_gf.MakeRef(wave_fwd->PressureFESpace(), x.GetBlock(1), 0);
         wave_obs.StateToObservation(p_gf, *(obs[l-1]));
      }
      
      // Visualization
      if (wave_vis.IsVis())
      {
         // Save data in the ParaView format every "VisSteps"-th step
         if (k % wave_vis.VisSteps() == 0 || k == n_steps)
         {
            paraview_dc->SetCycle(k);
            paraview_dc->SetTime(t_dim);
            // TODO: avoid writing mesh data at every cycle
            paraview_dc->Save();
         }
      }
   }
   cout << endl << "p2o Mult: Time-stepping took " << chrono.RealTime() << " seconds." << endl;
   cout <<         "p2o Mult: Average time/step: " << chrono.RealTime()/n_steps << " seconds." << endl;
   
   if (wave_vis.IsVis())
   {
      // paraview_dc->DeregisterField("velocity");
      paraview_dc->DeregisterField("pressure");
   }
   
   // 8. Compute the error if the solution is known
   if (WaveSolution::IsKnown())
   {
      double t_final = n_steps * dt;
      if (fabs(t_final-t_curr) > 1.0e-10)
      {
         MFEM_WARNING("t_final = " << t_final << " but t_curr = " << t_curr);
      }
      u_gf.MakeRef(wave_fwd->VelocityFESpace(), x.GetBlock(0), 0);
      p_gf.MakeRef(wave_fwd->PressureFESpace(), x.GetBlock(1), 0);
   
      int order_quad = max(2, 2*wave_fwd->PressureFESpace()->GetMaxElementOrder());
      const IntegrationRule *irs[Geometry::NumGeom];
      for (int i = 0; i < Geometry::NumGeom; i++)
      {
         irs[i] = &(IntRules.Get(i, order_quad));
      }
   
      VectorFunctionCoefficient u_coeff(dim, WaveSolution::uExact);
      FunctionCoefficient p_coeff(WaveSolution::pExact);
      VectorFunctionCoefficient p_grad_coeff(dim, WaveSolution::pGradExact);
      u_coeff.SetTime(t_final);
      p_coeff.SetTime(t_final);
      p_grad_coeff.SetTime(t_final);
      
      Mesh *mesh = wave_fwd->ParamFESpace()->GetMesh();
      
      double err_u_l2  = u_gf.ComputeL2Error(u_coeff, irs);
      double norm_u_l2 = ComputeLpNorm(2., u_coeff, *mesh, irs);
      
      double err_p_l2      = p_gf.ComputeL2Error(p_coeff, irs);
      double err_p_grad_l2 = p_gf.ComputeGradError(&p_grad_coeff, irs);
      double err_p_h1 = sqrt(err_p_l2*err_p_l2 + err_p_grad_l2*err_p_grad_l2);
      
      double norm_p_l2      = ComputeLpNorm(2., p_coeff, *mesh, irs);
      double norm_p_grad_l2 = ComputeLpNorm(2., p_grad_coeff, *mesh, irs);
      double norm_p_h1 = sqrt(norm_p_l2*norm_p_l2 + norm_p_grad_l2*norm_p_grad_l2);
      
      cout << endl << "Relative error:" << endl;
      cout << "|| u_h - u_ex || / || u_ex || = " << err_u_l2 / norm_u_l2 << endl;
      cout << "|| p_h - p_ex || / || p_ex || = " << err_p_h1 / norm_p_h1 << endl << endl;
   }
}

void WaveParamToObs::Mult(GridFunction **param, Vector** &obs) const
{
   MFEM_VERIFY(wave_fwd != nullptr, "No forward op available.");

   // Store parameters as GridFunctionCoefficients in wave operator
   wave_fwd->StoreLoad(param);
   
   // Perform Mult operator
   GetObs(obs);
}

void WaveParamToObs::Mult(function<double(const Vector &, double)> TDF, Vector** &obs) const
{
   MFEM_VERIFY(wave_fwd != nullptr, "No forward op available.");
   
   // Assign parameter time-dependent function to wave operator
   wave_fwd->TDFunctionLoad(TDF);
   
   // Perform Mult operator
   GetObs(obs);
}

void WaveParamToObs::GetAdj(GridFunction** &adj) const
{
   StopWatch chrono;
   chrono.Start();
   
   // 1. Define the BlockStructure of the problem, i.e. define the array of
   //    offsets for each variable. The last component of the Array is the sum
   //    of the dimensions of each block.
   Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = wave_adj->VelocityFESpace()->GetVSize();
   block_offsets[2] = wave_adj->PressureFESpace()->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "WaveParamToObs::MultTranspose \n";
   std::cout << "dim(DG)    = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(CG)    = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(DG+CG) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n";
   
   // 2. Allocate memory (x, rhs) for the analytical solution and the right hand
   //    side.  Define the GridFunction u,p for the finite element solution and
   //    linear forms fform and gform for the right hand side.  The data
   //    allocated by x and rhs are passed as a reference to the grid functions
   //    (u,p) and the linear forms (fform, gform).
   //MemoryType mt = device.GetMemoryType();
   BlockVector x(block_offsets); //, mt);
   
   // 3. Set the initial conditions for (u,p)
   
   // Initial condition (homogeneous)
   x = 0.0;

   // For verification purposes (when using manufactured solution),
   // the initial condition must satisfy exact manufactured solution
   GridFunction u_gf(wave_adj->VelocityFESpace());
   GridFunction p_gf(wave_adj->PressureFESpace());

   // 4. Initiate visualization data in the ParaView format and
   //    write the initial solution to file
   ParaViewDataCollection *paraview_dc;
   if (wave_vis.IsVis())
   {
      ostringstream collection_oss;
      collection_oss << "AdjointSolution_" << setfill('0') << setw(6) << vis_adj++;
      string collection = collection_oss.str();
      wave_vis.NewCollection(collection, wave_adj->PressureFESpace()->GetMesh());
      paraview_dc = wave_vis.ParaviewDC();
      u_gf.MakeRef(wave_adj->VelocityFESpace(), x.GetBlock(0), 0);
      p_gf.MakeRef(wave_adj->PressureFESpace(), x.GetBlock(1), 0);
      // paraview_dc->RegisterField("velocity_adj",&u_gf);
      paraview_dc->RegisterField("pressure_adj",&p_gf);
      paraview_dc->SetCycle(0);
      paraview_dc->SetTime(n_steps*dt * Cascadia::t0);
      chrono.Clear();
      paraview_dc->Save();
      cout << endl << "p2o MultTranspose: paraview I/O took " << chrono.RealTime() << " seconds." << endl;
   }

   // 5. Initialize ODE solver with adjoint wave operator
   ode_solver.Init(*wave_adj);
   
   // 6. Allocate GridFunctions for output
   if (!adj)
   {
      adj = new GridFunction*[param_steps];
      for (int k = 0; k < param_steps; k++) { adj[k] = nullptr; }
   }
   for (int k = 0; k < param_steps; k++)
   {
      if (adj[k]) { delete adj[k]; }
      adj[k] = new GridFunction(wave_adj->ParamFESpace());
   }
   
   // TODO: verify
//   GridFunction tmp(wave_adj->ParamFESpace());
//   BilinearForm m_var(wave_adj->ParamFESpace());
//   m_var.AddDomainIntegrator(new MassIntegrator(*(new ConstantCoefficient(-WaveSolution::c2))));
//   m_var.Assemble();
//   m_var.Finalize();
//   SparseMatrix *M = &(m_var.SpMat());
   
   // 7. Perform time-stepping
   chrono.Clear();

//   double t_curr = 0.0;
   double t_curr = n_steps*dt; // adjoint steps backward in time, starting at t=t_final
   for (int k = 1; k <= n_steps; k++)
   {
      // perform time step (increments t_curr)
      double dtc = -dt;
      wave_adj->Cycle(); // necessary for updating GridFunctionCoefficient (if store_load)
      ode_solver.Step(x, t_curr, dtc);

      if (abs(dtc+dt) > 1.0e-14)
      {
         MFEM_WARNING("abs(dtc-dt) = " << abs(dtc-dt));
      }
      
      double t_dim = t_curr * Cascadia::t0;
      if (k % 5 == 0)
      {
         cout << "=========================================" << endl
              << "time step = " << k << ", t = "
                       << setprecision(2) << t_dim << " s"   << endl
              << "=========================================" << endl;
      }
      
      // Extract output every "param_rate"-th step
      if (k % param_rate == 0)
      {
         int l = k / param_rate;
         int m = param_steps - l; // (param_steps-1) - (l-1);
         p_gf.MakeRef(wave_adj->PressureFESpace(), x.GetBlock(1), 0);
         wave_map.StateToParam(p_gf, *(adj[m]));
         // TODO: multiply adj[m] with "-c2"
         // TODO: multiply parameter space mass matrix with adj[m]
         // TODO: adj[m] = M * adj[m]
         *(adj[m]) *= -WaveSolution::c2;
//         M->Mult(tmp, *(adj[m]));
//         *(adj[m]) = tmp;
      }
      
      // Visualization
      if (wave_vis.IsVis())
      {
         // Save data in the ParaView format every "VisSteps"-th step
         if (k % wave_vis.VisSteps() == 0 || k == n_steps)
         {
            paraview_dc->SetCycle(k);
            paraview_dc->SetTime(t_dim);
            // TODO: avoid writing mesh data at every cycle
            paraview_dc->Save();
         }
      }
      
      if (k % 10 == 0)
      {
         cout << endl << " Time-stepping until now took " << chrono.RealTime() << " seconds." << endl;
         cout <<         " Average time/step until now: " << chrono.RealTime()/k << " seconds." << endl;
      }
   }
   cout << endl << "p2o MultTranspose: Time-stepping took " << chrono.RealTime() << " seconds." << endl;
   cout <<         "p2o MultTranspose: Average time/step: " << chrono.RealTime()/n_steps << " seconds." << endl;
   
   if (wave_vis.IsVis())
   {
      // paraview_dc->DeregisterField("velocity_adj");
      paraview_dc->DeregisterField("pressure_adj");
   }
}

/// Forward and adjoint system are equivalent except for the data (load),
/// thus forward solver can be used to solve the adjoint problem;
/// currently, the adjoint system is scaled s.t. u = c1*tau, p = c2*v
void WaveParamToObs::MultTranspose(Vector **data, GridFunction** &adj) const
{
   MFEM_VERIFY(wave_adj != nullptr, "No adjoint op available.");
   
   // Store data as GridFunctionCoefficients in wave operator
   wave_adj->StoreAdjLoad(data);
   
   // Perform Mult operator
   GetAdj(adj);
}

/// Evaluate misfit part of cost functional (between model observations and data)
double WaveParamToObs::EvalMisfit(Vector** &obs, Vector** &data) const
{
   double misfit = 0;
   
   for (int k = 0; k < obs_steps; k++)
   {
      Vector &obs_vec = *(obs[k]);
      Vector &data_vec = *(data[k]);
      Vector misfit_vec(obs_vec);
      misfit_vec -= data_vec;
      misfit += InnerProduct(misfit_vec, misfit_vec);
   }
   
   return misfit;
}
   
/// Add noise to observations
double WaveParamToObs::AddNoise(Vector** &obs, double rel_noise) const
{
   double max_val = 0.0;
   for (int k = 0; k < obs_steps; k++)
   {
      Vector &obs_vec = *(obs[k]);
      max_val = max(max_val, obs_vec.Normlinf());
   }
   double noise_std_dev = rel_noise * max_val;
   
   std::random_device rnd_dev;
   std::default_random_engine rnd_gen(rnd_dev());
   std::normal_distribution<double> distr(0.0, noise_std_dev);
   
   for (int k = 0; k < obs_steps; k++)
   {
      Vector &obs_vec = *(obs[k]);
      for (int i = 0; i < n_obs; i++)
      {
         obs_vec[i] += distr(rnd_gen);
      }
   }
   
   return noise_std_dev * noise_std_dev;
}

GridFunction** WaveParamToObs::ParamToGF(function<double(const Vector &, double)> TDF) const
{
   MFEM_VERIFY(wave_fwd != nullptr, "No forward operator has been defined.");
   
   GridFunction** params = new GridFunction*[param_steps];
   FunctionCoefficient coeff(TDF);
   double t = 0;
   for (int k = 0; k < param_steps; k++)
   {
      t = k*param_rate*dt;
      coeff.SetTime(t);
      params[k] = new GridFunction(wave_fwd->ParamFESpace());
      params[k]->ProjectCoefficient(coeff);
   }
   return params;
}

} // close namespace mfem
