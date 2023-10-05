#ifndef CASCADIA_WAVE_MAP
#define CASCADIA_WAVE_MAP

#include "mfem.hpp"

#include "wave_sol.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

/// Define mappings between parameter-space GF and state-space GF
class WaveTransferMap
{
protected:
   FiniteElementSpace &fes_state; // state FE space
   FiniteElementSpace &fes_param; // parameter FE space
   int order_state;
   int order_param;
   Mesh *state_mesh;
   SubMesh *param_mesh;
   FiniteElementCollection *fec_aux; // auxiliary FE collection
   FiniteElementSpace *fes_aux;      // auxiliary FE space
   TransferMap *state_to_aux;        // mesh -> submesh
   TransferMap *aux_to_state;        // submesh -> mesh
public:
   WaveTransferMap(FiniteElementSpace &fes_p_, FiniteElementSpace &fes_m_);

   /// Map parameter-space GF to state-space GF
   void ParamToState(const GridFunction &param_gf, GridFunction &state_gf);

   /// Map state-space GF to parameter-space GF
   void StateToParam(const GridFunction &state_gf, GridFunction &param_gf);
   
   ~WaveTransferMap();
};

}

#endif
