#ifndef CASCADIA_WAVE_P2O
#define CASCADIA_WAVE_P2O

#include "mfem.hpp"

#include "wave_sol.hpp"
#include "wave_map.hpp"
#include "wave_obs.hpp"
#include "wave_op.hpp"
#include "wave_vis.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/// Define mapping from parameter field to observables (p2o map)
/// and its adjoint
class WaveParamToObs
{
protected:
   WaveOperator *wave_fwd;      // Forward wave operator
   WaveOperator *wave_adj;      // Adjoint wave operator
   WaveObservationOp &wave_obs; // Observation operator
   WaveTransferMap &wave_map;   // Transfer map (param<->state)
   WaveVis &wave_vis;           // Visualization object
   ODESolver &ode_solver;
   
   int n_steps; // number of time steps
   double dt;   // time step size
   int n_param; // parameters per time step
   int n_obs;   // number of sensors
   int dim;     // dimension
   
   int param_rate, param_steps; // parameter frequency
   int obs_rate, obs_steps;     // observation frequency
   
   mutable int vis_fwd;  // counter for vis directory (forward vis)
   mutable int vis_adj;  // counter for vis directory (adjoint vis)

public:
   WaveParamToObs(WaveOperator *wave_fwd_, WaveOperator *wave_adj_,
                  WaveObservationOp &wave_obs_, WaveTransferMap &wave_map_,
                  WaveVis &wave_vis_,
                  ODESolver &ode_solver_, int n_steps_, double dt_);
   
   /// Auxiliary method for operator mult, used after parameter has been set in wave operator
   void GetObs(Vector** &obs) const;
   
   /// Forward operator using GridFunction parameter input
   void Mult(GridFunction **param, Vector** &obs) const;
   
   /// Forward operator using Function parameter input
   void Mult(std::function<double(const Vector &, double)> TDF, Vector** &obs) const;
   
   /// Auxiliary method for operator mult transpose, used after data has been set in wave operator
   void GetAdj(GridFunction** &adj) const;
   
   /// Adjoint operator using data input
   void MultTranspose(Vector **data, GridFunction** &adj) const;
   
   /// Evaluate misfit part of cost functional (between model observations and data)
   double EvalMisfit(Vector** &obs, Vector** &data) const;
   
   /// Add noise to observations; returns noise covariance
   double AddNoise(Vector** &obs, double rel_noise) const;
   
   /// Project a TD function onto GridFunctions for each time step
   GridFunction** ParamToGF(std::function<double(const Vector &, double)> TDF) const;

   ~WaveParamToObs() { }
};

}

#endif
