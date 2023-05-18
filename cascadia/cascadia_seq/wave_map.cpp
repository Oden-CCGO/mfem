
#include "cascadia.hpp"

namespace mfem
{

using namespace std;

WaveTransferMap::WaveTransferMap(FiniteElementSpace &fes_state_,
                                 FiniteElementSpace &fes_param_)
   :  fes_state(fes_state_), fes_param(fes_param_),
      state_mesh(NULL), param_mesh(NULL),
      fec_aux(NULL), fes_aux(NULL),
      state_to_aux(NULL), aux_to_state(NULL)
{
   MFEM_VERIFY(!fes_state.IsVariableOrder(),
               "WaveTransferMap requires uniform order.");
   MFEM_VERIFY(!fes_param.IsVariableOrder(),
               "WaveTransferMap requires uniform order.");
   
   // Get order of approximation of each FE space
   order_state = fes_state.FEColl()->GetOrder();
   order_param = fes_param.FEColl()->GetOrder();
   
   // Get meshes that each FE space is defined on
   state_mesh = fes_state.GetMesh();
   param_mesh = dynamic_cast<SubMesh *>(fes_param.GetMesh());
   
   MFEM_VERIFY(param_mesh->GetParent() == state_mesh,
               "WaveTransferMap: Parameter Mesh must be SubMesh of State Mesh.");
   
   // Create auxiliary FE space if needed
   if (order_param != order_state)
   {
      fec_aux = new H1_FECollection(order_state, param_mesh->Dimension());
      fes_aux = new FiniteElementSpace(param_mesh, fec_aux);
   }
   else
   {
      fes_aux = &fes_param;
   }
   
   // Create GridFunction TransferMaps
   GridFunction state_gf(&fes_state);
   GridFunction aux_gf(fes_aux);
   
   aux_to_state = new TransferMap(aux_gf, state_gf);
   state_to_aux = new TransferMap(state_gf, aux_gf);
}

void WaveTransferMap::ParamToState(const GridFunction &param_gf, GridFunction &state_gf)
{
   GridFunction aux_gf(fes_aux);
   aux_gf.ProjectGridFunction(param_gf);
   aux_to_state->Transfer(aux_gf, state_gf);
}

void WaveTransferMap::StateToParam(const GridFunction &state_gf, GridFunction &param_gf)
{
   GridFunction aux_gf(fes_aux);
   state_to_aux->Transfer(state_gf, aux_gf);
   param_gf.ProjectGridFunction(aux_gf);
}

WaveTransferMap::~WaveTransferMap()
{
   if (order_param != order_state)
   {
      delete fec_aux;
      delete fes_aux;
   }
   delete state_to_aux;
   delete aux_to_state;
}

}
