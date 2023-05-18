#ifndef CASCADIA_WAVE_OBS
#define CASCADIA_WAVE_OBS

#include "mfem.hpp"

#include "wave_sol.hpp"
#include "wave_map.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/// Define mappings between state space and observables (data space).
/// Note that observations are made on the boundary of the domain,
/// thus the state space variable is first mapped to the parameter space,
/// which is defined on that boundary (sea bottom).
/// Further, note that we make the following assumptions:
///   the parameter space is a linear nodal space;
///   the observation points coincide with mesh vertices.
class WaveObservationOp
{
protected:
   FiniteElementSpace &fes_param; // parameter FE space
   WaveTransferMap &wave_map;     // defines maps between state and parameter space
   DenseMatrix &sensor_pts;       // sensor locations (physical coords)
   Array<int> sensor_dofs;        // dof index corresponding to each sensor
   int n_obs;                     // number of observers
public:
   WaveObservationOp(FiniteElementSpace &fes_param_,
                     WaveTransferMap &wave_map_,
                     DenseMatrix &sensor_pts_);
   
   /// Apply observation operator
   void StateToObservation(const GridFunction &gf_state, Vector &obs) const;
   
   /// Apply transpose of observation operator
   void ObservationToState(const Vector &obs, GridFunction &gf_state) const;
   
   /// Get number of sensors
   int GetNrSensors() const { return n_obs; };
   
   ~WaveObservationOp() { };
};

}

#endif
