
#include "cascadia.hpp"

namespace mfem
{

using namespace std;

WaveVis::WaveVis(Mesh *mesh_, bool visualization_, int vis_steps_, int vis_order_,
                 string &collection_name_)
   : visualization(visualization_), vis_steps(vis_steps_), vis_order(vis_order_),
     mesh(mesh_), paraview_dc(nullptr)
{
   if (visualization)
   {
      NewCollection(collection_name_, mesh);
   }
}

void WaveVis::NewCollection(string &collection_name_, Mesh *mesh_)
{
   if (visualization)
   {
      delete paraview_dc;
      if (mesh_) { mesh = mesh_; }
      collection_name = collection_name_;
      paraview_dc = new ParaViewDataCollection(collection_name, mesh);
      string pv = WaveIO::output_dir + "/ParaView";
      paraview_dc->SetPrefixPath(pv);
      paraview_dc->SetLevelsOfDetail(vis_order);
      paraview_dc->SetDataFormat(VTKFormat::BINARY);
      paraview_dc->SetHighOrderOutput(true);
   }
}

}
