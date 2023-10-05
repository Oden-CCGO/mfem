
#include "cascadia.hpp"

namespace mfem
{

using namespace std;

/// To set up the observation operator, we need to find out which nodal dof each sensor corresponds to.
/// Then, we can directly get the dof subvector of the GF to obtain the observational data,
///       and conversely set the subvector of the GF for the transpose of the observation operator.
/// Note: this implementation assumes that we get observations from a linear nodal FE space.
WaveObservationOp::WaveObservationOp(FiniteElementSpace &fes_param_,
                                     WaveTransferMap &wave_map_,
                                     DenseMatrix &sensor_pts_)
   : fes_param(fes_param_), wave_map(wave_map_),
     sensor_pts(sensor_pts_), n_obs(sensor_pts_.Width())
{
   MFEM_VERIFY(fes_param.FEColl()->GetOrder() == 1,
               "WaveObservationOp only supports linear nodal FE space.");

   const double eps = 1.0e-10;
   
   cout << endl << "Setting up Observation Operator ..." << endl;

   // Option 1:
   // Find element + integration point for each sensor location via mesh->FindPoints
   StopWatch chrono;
   chrono.Start();
   
   Array<int> sensor_elem_ids;
   Array<IntegrationPoint> sensor_ips;
   int n_found = fes_param.GetMesh()->FindPoints(sensor_pts, sensor_elem_ids, sensor_ips);
   MFEM_VERIFY(n_found == n_obs,
               "Found " << n_found << " / " << n_obs << " points.");
   
   sensor_dofs.SetSize(n_obs);
   sensor_dofs = -1;
   Array<int> elem_dofs;
   for (int i = 0; i < n_obs; i++)
   {
      // Get nodal FE in which sensor (i) is located
      auto fe_i = fes_param.GetFE(sensor_elem_ids[i]);
      
      // Get element dofs
      fes_param.GetElementDofs(sensor_elem_ids[i], elem_dofs);
      for (int j = 0; j < elem_dofs.Size(); j++)
      {
         // Get j-th integration point
         const IntegrationPoint &ip = fe_i->GetNodes().IntPoint(j);

         if (fabs(ip.x-sensor_ips[i].x) < eps &&
             fabs(ip.y-sensor_ips[i].y) < eps )
         {
            MFEM_VERIFY(fabs(ip.x) < eps || fabs(ip.x-1) < eps,
                        "Integration point is not at nodal point.");
            MFEM_VERIFY(fabs(ip.y) < eps || fabs(ip.y-1) < eps,
                        "Integration point is not at nodal point.");
            
            // This integration point corresponds to the sensor i-th location.
            // Since element dof ordering is the same as integration point ordering
            // for nodal FE, this nodal dof corresponds to the sensor location.
            sensor_dofs[i] = elem_dofs[j];
         }
      }
   }
   MFEM_VERIFY(sensor_dofs.Min() >= 0,
               "Did not find nodal point for a sensor.");
   
   cout << "   Option 1: " << chrono.RealTime() << " seconds." << endl;
   
   // Verification of Option 1
   GridFunction param_gf(&fes_param);
   FunctionCoefficient linear_coeff(WaveSolution::LinearFunction);
   param_gf.ProjectCoefficient(linear_coeff);
   
   double sensor_data, aux_data;
   for (int i = 0; i < n_obs; i++)
   {
      sensor_data = param_gf[sensor_dofs[i]];
      aux_data = param_gf.GetValue(sensor_elem_ids[i], sensor_ips[i]);
      MFEM_VERIFY(fabs(aux_data-sensor_data) < eps,
                  "Verification for sensor data failed.");
   }
   
   // Option 2:
   // Iterate through the mesh elements + integration points
   // and match IP to sensor locations directly
   
   chrono.Clear();
   
   Array<int> aux_dofs(n_obs);
   Vector coord;
   for (int i = 0; i < fes_param.GetMesh()->GetNE(); i++)
   {
      // Get i-th nodal FE
      auto fe_i = fes_param.GetFE(i);
   
      // Get element transformation
      auto elem_trans = fes_param.GetElementTransformation(i);
      
      // Get element dofs
      fes_param.GetElementDofs(i, elem_dofs);
      
      for (int j = 0; j < elem_dofs.Size(); j++)
      {
         // Get j-th integration point
         const IntegrationPoint &ip = fe_i->GetNodes().IntPoint(j);
         
         // Get physical coordinates of integration point
         elem_trans->Transform(ip, coord);
//         cout << "j = " << j << ":  "; coord.Print();
         
         // Search if the IP coordinates match any sensor location
         for (int k = 0; k < n_obs; k++)
         {
            if (fabs(coord[0]-sensor_pts.Elem(0,k)) < eps &&
                fabs(coord[1]-sensor_pts.Elem(1,k)) < eps)
            {
               MFEM_VERIFY(fabs(ip.x) < eps || fabs(ip.x-1) < eps,
                        "Integration point is not at nodal point.");
               MFEM_VERIFY(fabs(ip.y) < eps || fabs(ip.y-1) < eps,
                        "Integration point is not at nodal point.");
               
               // This IP coincides with the k-th sensor location.
               // Since element dof ordering is the same as integration point ordering
               // for nodal FE, this nodal dof corresponds to the sensor location.
               aux_dofs[k] = elem_dofs[j];
            }
         }
      }
   }
   
   cout << "   Option 2: " << chrono.RealTime() << " seconds." << endl;
   
   // Verification of Option 2
   for (int i = 0; i < n_obs; i++)
   {
      MFEM_VERIFY(sensor_dofs[i] == aux_dofs[i],
                  "Verification for sensor dofs option 2 failed.");
   }
   
   cout << "WaveObservationOp constructor done. All verifications passed." << endl << endl;
}


void WaveObservationOp::StateToObservation(const GridFunction &state_gf, Vector &obs) const
{
   obs.SetSize(n_obs);
   GridFunction param_gf(&fes_param);
   wave_map.StateToParam(state_gf, param_gf);
   param_gf.GetSubVector(sensor_dofs, obs);
}

void WaveObservationOp::ObservationToState(const Vector &obs, GridFunction &state_gf) const
{
   MFEM_VERIFY(obs.Size() == n_obs,
               "obs vector size must equal the number of observation points.");
   GridFunction param_gf(&fes_param);
   param_gf = 0;
   param_gf.SetSubVector(sensor_dofs, obs);
   wave_map.ParamToState(param_gf, state_gf);
}

}
