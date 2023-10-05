#ifndef CASCADIA_WAVE_IO
#define CASCADIA_WAVE_IO

#include "mfem.hpp"

#include "wave_obs.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/// Define functions and variables for file I/O
class WaveIO
{
protected:
   FiniteElementSpace* param_space;
   WaveObservationOp &wave_obs;

   int n_steps; // number of time steps
   double dt;   // time step size
   int n_param; // parameters per time step
   int n_obs;   // number of sensors
   int dim;     // dimension
   
   int param_rate, param_steps; // parameter frequency
   int obs_rate, obs_steps;     // observation frequency
   
   std::string rel_path;
   std::string prefix_fwd;
   std::string prefix_adj;
   bool init_meta_fwd;
   bool init_meta_adj;
   int count_fwd_text;   // counter for h5 file suffix (forward map)
   int count_fwd_binary; // counter for txt file suffix (forward map)
   bool memcpy; // specifies whether FwdToFile/AdjToFile copy sub-vectors
                // into a global vector before writing to dataset
   bool reverse_order; // specifies whether adj_vec are written in block-reverse
                       // order as is needed by the FFT matvec code
   
   bool binary;
   
   /// Write metadata to file
   void MetaToFile(bool adj);

public:
   WaveIO(FiniteElementSpace* param_space_, WaveObservationOp &wave_obs_,
          int n_steps_, double dt_, int n_param_, int n_obs_,
          int param_rate_, int obs_rate_, bool binary_=false);
   
   /// Write output of forward map to file (data space)
   void FwdToFile(Vector **obs);
   
   /// Specifies if writing adjoint vectors in block-reverse ordering
   static bool adj_reverse_order;
   
   /// Write output of adjoint map to file (parameter space)
   void AdjToFile(GridFunction **adj, int adjvec=0);
   
   /// Write observations to file
   void ObsToFile(const std::string &filename, Vector **obs,
                  double rel_noise=0.0, double noise_cov=1.0);
   
   /// Read observations from file
   Vector** ObsFromFile(const std::string &filename);
   
   /// Write parameter to file
   void ParamToFile(const std::string &filename, GridFunction **param);
   
   /// Read parameter from file
   GridFunction** ParamFromFile(const std::string &filename);
   
   /// Relative path to parent folder for all outputs
   static std::string output_dir;
   
   /// Helper function
   static int CreateDirectory(const std::string &dir_name,
                              const Mesh *mesh, int myid);

   ~WaveIO() { }
};

}

#endif
