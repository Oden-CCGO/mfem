#ifndef CASCADIA_WAVE_VIS
#define CASCADIA_WAVE_VIS

#include "mfem.hpp"

#include "wave_sol.hpp"
#include "wave_map.hpp"
#include "wave_obs.hpp"
#include "wave_op.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace mfem
{

class WaveVis
{
protected:
   const bool visualization;
   const int vis_steps;
   const int vis_order;
   
   Mesh *mesh;
   ParaViewDataCollection *paraview_dc;
   
   std::string collection_name;
public:
   WaveVis(Mesh *mesh_, bool visualization_, int vis_steps_, int vis_order_,
           std::string &collection_name_);
   
   /// Create a new DataCollection (deleting previous one if any)
   void NewCollection(std::string &collection_name_, Mesh *mesh_=nullptr);
   
   bool IsVis() const { return visualization; }
   int VisSteps() const { return vis_steps; }
   ParaViewDataCollection* ParaviewDC() const { return paraview_dc; }
   
   ~WaveVis() { delete paraview_dc; }
};

}

#endif
