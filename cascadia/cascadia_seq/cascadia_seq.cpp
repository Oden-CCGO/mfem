//                          MFEM Cascadia Code
//                         (sequential version)
//
// Compile with: make cascadia_seq
//
// Sample runs:  cascadia_seq
//
// Description:  This example code solves a 3D mixed acoustic wave problem
//               corresponding to the following mixed formulation:
//
//                           du/dt + c1 * grad p = f
//                           dp/dt + c2 *  div u = g
//
//               with natural boundary condition -p = <given pressure>.
//               Here, we use a given exact solution (u,p) and compute the
//               corresponding r.h.s. (f,g).  We discretize continuous (H1)
//               finite elements (scalar-valued pressure p) and piecewise
//               discontinuous (L2) polynomials (vector-valued velocity u).
//
//
//  Boundary attributes for different parts of the boundary:
//  (using convention from Mesh::Make3D)
//
//                                       top (zmax)             back (ymax)
//                                       bdr attr 6             bdr attr 4
//                                       [surfac_bdr]           [absorb_bdr]
//                                              |              /
//                                              |             /
//                                              v
//                                    +---------------------+
//                                   /                     /|
//                                  /                     / |
//            left (xmin)   --->   +---------------------+  +      <--- right (xmax)
//            bdr attr 5           |                     | /            bdr attr 3
//            [absorb_bdr]         |                     |/             [absorb_bdr]
//                                 +---------------------+
//                                              ^
//                                 /            |
//                                /             |
//                      front (ymin)       bottom (zmin)
//                      bdr attr 2         bdr attr 1
//                      [absorb_bdr]       [bottom_bdr]
//

#include "cascadia.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace mfem;

/// Temporary
double EvalMisfit(Vector** &obs, Vector** &data, int obs_steps);
double AddNoise(Vector** &obs, double rel_noise, int obs_steps);

int main(int argc, char *argv[])
{
   StopWatch chrono;
   chrono.Start();

   // 1. Parse command-line options.
   const char *mesh_file = "";
   const char *output_dir = ".";
   int problem = 1;
   int fwd = 1;
   int adj = 0;
   int revadj = 0;
   int adjvec = -1;
   int prior = 0;
   int indprior = 0;
   double alpha1 = 1.0;
   double alpha2 = 1.0;
   double alpha3 = 0.0;
   int order = 2;
   int ode_solver_type = 11;
   bool lump = false;
   int ref_levels = 1;
   double t_final = 1.0;
   int n_steps = 1;
   int param_rate = 1;
   int obs_rate = 1;
   int nx_obs = 3;
   int ny_obs = 3;
   bool observations = true;
   bool visualization = true;
   int vis_steps = 5;
   bool hdf = true;
   bool pa = false;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&output_dir, "-d", "--directory",
                  "Relative path of output directory for I/O.");
   args.AddOption(&problem, "-p", "--problem",
                  "Problem: 1   - known stationary solution: sine-wave,\n\t"
                  "         2   - known stationary solution: lin polynomial,\n\t"
                  "         3   - known stationary solution: quadr polynomial,\n\t"
                  "         10  - known TD solution: exp decaying sine-wave,\n\t"
                  "         20  - known TD solution: exp decaying lin polynomial,\n\t"
                  "         30  - known TD solution: exp decaying quadr polynomial,\n\t"
                  "         40  - known TD solution: lin growing quadr polynomial,\n\t"
                  "         50  - known TD solution: quadr growing quadr polynomial),\n\t"
                  "         100 - unknown solution: forcing with single (x,y) Gaussian deformation,\n\t"
                  "  [N/A]  200 - unknown solution: forcing with superimposed Gaussian deformations)");
   args.AddOption(&fwd, "-fwd", "--forward",
                  "Forward solve: 0 - Disable forward operator,\n\t"
                  "               1 - Enable forward operator (one solve),\n\t"
                  "               2 - Use forward operator to export p2o map.");
   args.AddOption(&adj, "-adj", "--adjoint",
                  "Adjoint solve: 0 - Disable adjoint operator,\n\t"
                  "               1 - Enable adjoint operator (one solve),\n\t"
                  "               2 - Use adjoint operator to export adjoint p2o map.");
   args.AddOption(&revadj, "-revadj", "--adjoint-reverse-order",
                  "Adjoint reverse order: 0 - Write adjoint vectors in standard ordering,\n\t"
                  "                       1 - Write adjoint vectors in block-reverse ordering.");
   args.AddOption(&adjvec, "-adjvec", "--adjoint-vec",
                  "Adjoint p2o vectors:  -1 - Compute all adjoint vectors (if adj=2),\n\t"
                  "                     >=0 - Compute this particular adjoint vector (if adj=2).");
   args.AddOption(&prior, "-prior", "--prior",
                  "Prior: 0 - Do not assemble prior,\n\t"
                  "       1 - Laplacian prior (assemble),\n\t"
                  "      11 - Laplacian prior (assemble + write to file),\n\t"
                  "       2 - Bi-Laplacian prior (assemble),\n\t"
                  "      22 - Bi-Laplacian prior (assemble + write to file).");
   args.AddOption(&indprior, "-indprior", "--prior-reindex",
                  "Prior reindex: 0 - Write prior with indexing as time(outer)->space(inner),\n\t"
                  "               1 - Write prior re-indexed as space(outer)->time(inner).");
   args.AddOption(&alpha1, "-alpha1", "--alpha1",
                  "Regularization parameter for scaling |m|.");
   args.AddOption(&alpha2, "-alpha2", "--alpha2",
                  "Regularization parameter for scaling |grad m|.");
   args.AddOption(&alpha3, "-alpha3", "--alpha3",
                  "Regularization parameter for scaling |dm/dt|.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&ode_solver_type, "-ode", "--ode-solver",
                  "ODE solver: 1 - Forward Euler,\n\t"
                  "            2 - RK2, 3 - RK3 SSP,\n\t"
                  "            4 - RK4, 6 - RK6,\n\t"
                  "            11 - Backward Euler");
   args.AddOption(&lump, "-lump", "--lump-mass", "-no-lump", "--no-lump-mass",
                  "Enable or disable lumping mass matrices.");
   args.AddOption(&ref_levels, "-ref", "--ref-levels",
                  "Number of uniform h-refinements.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&n_steps, "-nt", "--n-steps",
                  "Number of time steps.");
   args.AddOption(&param_rate, "-pr", "--param-rate",
                  "Parameter is defined for every n-th time step.");
   args.AddOption(&obs_rate, "-or", "--obs-rate",
                  "Observations are defined for every n-th time step.");
   args.AddOption(&nx_obs, "-nxo", "--nx-obs",
                  "Number of observers in x-direction.");
   args.AddOption(&ny_obs, "-nyo", "--ny-obs",
                  "Number of observers in y-direction.");
   args.AddOption(&observations, "-obs", "--observations", "-no-obs", "--no-observations",
                  "Enable or disable writing observations files.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization",
                  "Enable or disable writing paraview visualization files.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Write visualization files every n-th timestep.");
   args.AddOption(&hdf, "-hdf", "--hdf-format", "-no-hdf", "--no-hdf-format",
                  "Enable or disable HDF5 binary output format.");
                  
   // PA and Device options not yet supported in this code
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa", "--no-partial-assembly",
                  "Enable or disable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   cout << endl;
   args.PrintOptions(cout);
   
   // Set global variables
   WaveSolution::t_final = t_final;
   WaveSolution::n_steps = n_steps;
   WaveSolution::dt = t_final/n_steps;
   
   WaveIO::output_dir = output_dir;
   if (revadj) { WaveIO::adj_reverse_order = true; }
   if (indprior) { WavePrior::reindex = true; }
   
   // Checks for param_rate, obs_rate
   MFEM_VERIFY(n_steps % param_rate == 0,
               "param_rate must evenly divide n_steps.");
   MFEM_VERIFY(n_steps % obs_rate == 0,
               "obs_rate must evenly divide n_steps.");
   MFEM_VERIFY(obs_rate % param_rate == 0,
               "obs_rate must be a multiple of param_rate.");
   
   const int param_steps = n_steps/param_rate;
   const int obs_steps   = n_steps/obs_rate;

   const int order_p = order;
   const int order_u = order-1;

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3a. Read-in or create mesh.
   Mesh *mesh;
   if (strlen(mesh_file) != 0)
   {
      mesh = new Mesh(mesh_file, 1, 1);
   }
   else
   {
      int nx = 1; double sx = 1.0;
      int ny = 1; double sy = 1.0;
      int nz = 1; double sz = 1.0;
      if (problem >= 100)
      {
         sx = 8.0; nx = 16; int xfac =  8;
         sy = 8.0; ny = 16; int yfac = 16;
         sz = 4.0; nz = 16;
         sx *= xfac; nx *= xfac;
         sy *= yfac; ny *= yfac;
         cout << endl << "Initial mesh element dimensions:" << endl;
         cout << "  sx/nx = " << (sx/nx)*(1000/Cascadia::l0) << " km" << endl;
         cout << "  sy/ny = " << (sy/ny)*(1000/Cascadia::l0) << " km" << endl;
         cout << "  sz/nz = " << (sz/nz)*(1000/Cascadia::l0) << " km" << endl;
      }
      mesh = new Mesh(Mesh::MakeCartesian3D(
                      nx, ny, nz, Element::HEXAHEDRON,
                      sx, sy, sz));
   }
   const int dim = mesh->Dimension();
   
   // 3b. Get mesh bounding box.
   Vector min_mesh(dim);
   Vector max_mesh(dim);
   mesh->GetBoundingBox(min_mesh, max_mesh, 0);
   WaveSolution::xmin = min_mesh(0);
   WaveSolution::xmax = max_mesh(0);
   WaveSolution::ymin = min_mesh(1);
   WaveSolution::ymax = max_mesh(1);
   WaveSolution::zmin = min_mesh(2);
   WaveSolution::zmax = max_mesh(2);
   
   // 4. Uniformly h-refine the mesh.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }
   
   cout << endl << "++++++++++++++++++++++++++++++"
        << endl << "Printing 3D Mesh Info" << endl;
   mesh->PrintCharacteristics();
   
   // 5. Define the ODE solver used for time integration.
   ODESolver *ode_solver = nullptr;
   switch (ode_solver_type)
   {
      // Explicit methods
      case 1: ode_solver = new ForwardEulerSolver; break;
      // RK2: a=1/2 (midpoint method); a=1 (Heun's method); a=2/3 (min truncation error)
      case 2: ode_solver = new RK2Solver(1.0); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;
      // Implicit (L-stable) methods
      case 11: ode_solver = new BackwardEulerSolver; break;
      default:
         MFEM_ABORT("Undefined ODE solver type: " << ode_solver_type);
   }
   
   // 6a. Define the finite element spaces for velocity u and pressure p.
   FiniteElementCollection *l2_coll(new L2_FECollection(order_u, dim));
   FiniteElementCollection *h1_coll(new H1_FECollection(order_p, dim));

   FiniteElementSpace *DG_space = new FiniteElementSpace(mesh, l2_coll, dim, Ordering::byNODES);
   FiniteElementSpace *CG_space = new FiniteElementSpace(mesh, h1_coll);
   
   cout << "DG_space->GetVSize() = " << DG_space->GetVSize() << endl;
   cout << "CG_space->GetVSize() = " << CG_space->GetVSize() << endl;
   
   cout << "DG_space->GetTrueVSize() = " << DG_space->GetTrueVSize() << endl;
   cout << "CG_space->GetTrueVSize() = " << CG_space->GetTrueVSize() << endl;
   cout << "++++++++++++++++++++++++++++++" << endl;

   // 6b. Mark boundary attributes (essential, natural BCs).
   Array<int> ess_tdof_list;
   MFEM_VERIFY(mesh->bdr_attributes.Size() == 6,
               "Expect 6 boundary attributes for this problem.");
   Array<int> ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 0; // 0: natural BCs ; 1: essential BCs
   CG_space->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   MFEM_VERIFY(ess_tdof_list.Size() == 0,
               "Using only natural BCs for this weak formulation.");
   
   Array<int> nat_bdr(6); nat_bdr = 1;
   Array<int> surfac_bdr(6); surfac_bdr = 0;
   Array<int> bottom_bdr(6); bottom_bdr = 0;
   Array<int> absorb_bdr(6); absorb_bdr = 0;
   for (int i=0; i < 6; ++i)
   {
      switch (mesh->bdr_attributes[i])
      {
         case 1:
            bottom_bdr[i] = 1;
            break;
         case 2: case 3: case 4: case 5:
            absorb_bdr[i] = 1;
            break;
         case 6:
            surfac_bdr[i] = 1;
            break;
         default:
            MFEM_ABORT("main : unexpected bdr attribute " << mesh->bdr_attributes[i]);
      }
   }
   
   // TEST: Boundary Markers
   bool test_boundary_markers = false;
   if (test_boundary_markers)
   {
      GridFunction xyz(CG_space); xyz = 0;
      ConstantCoefficient thr(3.0);
      xyz.ProjectBdrCoefficient(thr, absorb_bdr);
      ConstantCoefficient two(2.0);
      xyz.ProjectBdrCoefficient(two, surfac_bdr);
      ConstantCoefficient one(1.0);
      xyz.ProjectBdrCoefficient(one, bottom_bdr);
      ParaViewDataCollection paraview_dc_bm("ForwardSeqBM", mesh);
      string pv = WaveIO::output_dir + "/ParaView";
      paraview_dc_bm.SetPrefixPath(pv);
      paraview_dc_bm.SetLevelsOfDetail(1);
      paraview_dc_bm.SetDataFormat(VTKFormat::BINARY);
      paraview_dc_bm.RegisterField("pressure",&xyz);
      paraview_dc_bm.Save();
   }
   // END TEST

   // 6c. Create submesh to define the parameter space on the sea bottom
   //     First, extract SubMesh from Mesh via boundary attributes
   SubMesh *bottom_mesh = new SubMesh(SubMesh::CreateFromBoundary(*mesh, bottom_bdr));
   MFEM_VERIFY(bottom_mesh->Dimension() == dim-1,
               "Bottom boundary SubMesh dim must be one lower than Mesh dim.");
   cout << endl << "++++++++++++++++++++++++++++++"
        << endl << "Printing 2D Bottom_Mesh Info" << endl;
   bottom_mesh->PrintCharacteristics();
   
   //     Create a FE collection and FE space on the sea bottom SubMesh,
   //     using linear continuous elements for the parameter space
   const int order_m = 1;
   FiniteElementCollection *param_h1_coll(new H1_FECollection(order_m, dim-1));
   FiniteElementSpace *param_space = new FiniteElementSpace(bottom_mesh, param_h1_coll);
   cout << "param_space->GetVSize()     = " << param_space->GetVSize()     << endl;
   cout << "param_space->GetTrueVSize() = " << param_space->GetTrueVSize() << endl;
   cout << "++++++++++++++++++++++++++++++" << endl;
   const int n_param = param_space->GetVSize();
   
   //     Define transfer operators between H1 GridFunctions on Mesh and SubMesh
   WaveTransferMap wave_map(*CG_space, *param_space);
   
   GridFunction m_gf(param_space);
   GridFunction CG_gf(CG_space);
   
   // TEST: TransferMap: SubMesh<->Mesh (Parameter<->State)
   bool test_transfer_maps = false;
   if (test_transfer_maps)
   {
      m_gf = 0; CG_gf = 0;
      FunctionCoefficient lin(WaveSolution::LinearFunction);
      CG_gf.ProjectBdrCoefficient(lin, bottom_bdr);
      wave_map.StateToParam(CG_gf, m_gf);
      
      ParaViewDataCollection paraview_dc_sub("ForwardSeq_SubMesh", bottom_mesh);
      string pv = WaveIO::output_dir + "/ParaView";
      paraview_dc_sub.SetPrefixPath(pv);
      paraview_dc_sub.SetLevelsOfDetail(1);
      paraview_dc_sub.SetDataFormat(VTKFormat::BINARY);
      paraview_dc_sub.RegisterField("parameter",&m_gf);
      paraview_dc_sub.Save();
      
      m_gf = 0; CG_gf = 0;
      m_gf.ProjectCoefficient(lin);
      wave_map.ParamToState(m_gf, CG_gf);
      
      ParaViewDataCollection paraview_dc_mesh("ForwardSeq_Mesh", mesh);
      paraview_dc_mesh.SetPrefixPath(pv);
      paraview_dc_mesh.SetLevelsOfDetail(1);
      paraview_dc_mesh.SetDataFormat(VTKFormat::BINARY);
      paraview_dc_mesh.RegisterField("pressure",&CG_gf);
      paraview_dc_mesh.Save();
   }
   // END TEST
   
   // 6d. Create mappings between state space and data space:
   //     Map from pressure variable to pointwise observational data, and its transpose.

   MFEM_VERIFY(nx_obs > 1 && ny_obs > 1, "nx_obs, ny_obs");
   //     Set number of sensor locations
   const double xx_pts_min = 0.25*WaveSolution::xmax;
   const double xx_pts_max = 0.75*WaveSolution::xmax;
   const double xx_dist    = (xx_pts_max-xx_pts_min) / (nx_obs-1);
   const double yy_pts_min = 0.25*WaveSolution::ymax;
   const double yy_pts_max = 0.75*WaveSolution::ymax;
   const double yy_dist    = (yy_pts_max-yy_pts_min) / (ny_obs-1);
   const int n_obs = nx_obs*ny_obs; // number of observation points (sensors)
   
   //     Specify physical coordinates for each sensor
   double sensor_pts[n_obs][dim];
   for (int j = 0; j < ny_obs; j++)
   {
      double yy = yy_pts_min + j*yy_dist;
      for (int i = 0; i < nx_obs; i++)
      {
         double xx = xx_pts_min + i*xx_dist;
         sensor_pts[j*nx_obs+i][0] = xx;
         sensor_pts[j*nx_obs+i][1] = yy;
         sensor_pts[j*nx_obs+i][2] = 0;
      }
   }
   
   //     Create dense matrix with coordinates of observation points
   DenseMatrix pts_mat(&(sensor_pts[0][0]), dim, n_obs);
   cout << endl << "Observation points:" << endl;
   pts_mat.PrintT(cout, dim);
   
   //     Create Observation Operator (defines mappings between state and data space)
   WaveObservationOp wave_obs(*param_space, wave_map, pts_mat);

   // TEST: Mapping between: State space <-> Data space (pointwise observations)
   bool test_data_maps = false;
   if (test_data_maps)
   {
      GridFunction pressure_gf(CG_space);
      GridFunction adjoint_gf(CG_space);
      FunctionCoefficient lin(WaveSolution::LinearFunction);
      pressure_gf.ProjectCoefficient(lin);

      Vector obs_data(n_obs);
      
      wave_obs.StateToObservation(pressure_gf, obs_data);
      // obs_data.Print(cout, 1);
      wave_obs.ObservationToState(obs_data, adjoint_gf);
      
      Array<int> sensor_elem_ids;
      Array<IntegrationPoint> sensor_ips;
      int n_found = mesh->FindPoints(pts_mat, sensor_elem_ids, sensor_ips);
      MFEM_VERIFY(n_found == n_obs,
                  "test_data_maps: Found " << n_found << " / " << n_obs << " points.");

      double val1, val2;
      for (int i = 0; i < n_obs; i++)
      {
         val1 = pressure_gf.GetValue(sensor_elem_ids[i], sensor_ips[i]);
         val2 = adjoint_gf.GetValue(sensor_elem_ids[i], sensor_ips[i]);
         MFEM_VERIFY(fabs(val1-val2) < 1.0e-10,
                     "test_data_maps: val1 = " << val1 << ", val2 = " << val2);
      }
   }
   // END TEST

   // 7. Initialize wave operators, BCs, forcing terms, and analytical solution (if any)
   WaveSolution::Init(problem);
   
   int height = DG_space->GetTrueVSize() + CG_space->GetTrueVSize();
   
   WaveOperator *wave_fwd = nullptr;
   if (fwd)
   {
      cout << endl << "Creating WaveOperator (forward)..." << endl;
      wave_fwd = new WaveOperator(*DG_space, *CG_space, *param_space, wave_map, wave_obs,
         ess_bdr, nat_bdr, surfac_bdr, bottom_bdr, absorb_bdr, lump,
         height, WaveSolution::c1, WaveSolution::c2, WaveSolution::c3,
         param_rate, param_steps, obs_rate, obs_steps, false);
   }
   
   WaveOperator *wave_adj = nullptr;
   if (adj)
   {
      cout << endl << "Creating WaveOperator (adjoint)..." << endl;
      wave_adj = new WaveOperator(*DG_space, *CG_space, *param_space, wave_map, wave_obs,
         ess_bdr, nat_bdr, surfac_bdr, bottom_bdr, absorb_bdr, lump,
         height, WaveSolution::c1, WaveSolution::c2, WaveSolution::c3,
         param_rate, param_steps, obs_rate, obs_steps, true);
   }
   
   // 8. Find CFL-save time-step dt.
   
   //      Set CFL max parameter.
   const double cfl_max = 0.25;
   
   //      Determine the minimum element size.
   double h_min = mesh->GetElementSize(0, 1);
   for (int i = 1; i < mesh->GetNE(); i++)
   {
      h_min = min(mesh->GetElementSize(i, 1), h_min);
   }
   
   //      Determine maximal wave speed.
   const double c_wave = sqrt(WaveSolution::c1 * WaveSolution::c2);
   
   //      Determine CFL estimate for maximal time step.
   const double max_dt = cfl_max*(h_min/(3*c_wave))/(2*order+1);
   
   //      Print CFL estimate
   cout << "CFL estimate:" << endl
        << "    cfl_max = " << cfl_max << endl
        << "    h_min   = " << h_min*(Cascadia::l0/1000) << " km" << endl
        << "    c_wave  = " << c_wave*(Cascadia::u0/1000) << " km/s" << endl
        << "    max_dt  = " << max_dt*(Cascadia::t0) << " s" << endl << endl;
   
   //      Determine actual time-step size dt
   const double dt = WaveSolution::dt;
   cout << "Time-stepping parameters:" << endl
        << "    t_final = " << t_final*(Cascadia::t0) << " s" << endl
        << "    n_steps = " << n_steps << endl
        << "    dt      = " << dt*(Cascadia::t0) << " s" << endl
        << " param_rate = " << param_rate << endl
        << "   obs_rate = " << obs_rate << endl << endl;
   
   if (dt > max_dt)
   {
      MFEM_WARNING("CFL estimate warning: dt > max_dt");
   }
   
   // 9. Create vis object
   string wave_dc_name = "WaveSolution";
   WaveVis wave_vis(mesh, visualization, vis_steps, order, wave_dc_name);
   
   // 10a. Create p2o operator
   WaveParamToObs *wave_p2o = nullptr;
   if (wave_fwd || wave_adj)
   {
      wave_p2o = new WaveParamToObs(wave_fwd, wave_adj, wave_obs, wave_map, wave_vis,
                                    *ode_solver, n_steps, dt);
   }
   
   // 10b. Create I/O object
   WaveIO wave_io(param_space, wave_obs, n_steps, dt,
                  n_param, n_obs, param_rate, obs_rate, hdf);
   
   // TEST: Create loads for all time-steps using GridFunctionCoefficients
   //       -> only used for unknown solutions
   GridFunction **params = nullptr;
   bool test_store_load = true;
   if (test_store_load && WaveSolution::IsUnknown() && fwd==1)
   {
      cout << "Storing p2o fwd load (parameter) via GridFunctionCoefficients." << endl << endl;
      params = wave_p2o->ParamToGF(WaveSolution::mParameter);
   }
   else if (fwd==1)
   {
      cout << "Defining p2o fwd load (parameter) via time-dependent function." << endl << endl;
   }
   // END TEST
   
   // Read parameter from file
   bool param_from_file = false;
   if (param_from_file)
   {
      if (params)
      {
         for (int k = 0; k < param_steps; k++) { delete params[k]; params[k] = nullptr; }
         delete[] params; params = nullptr;
      }
      string filename(output_dir);
      filename += (hdf) ? "/param_vec.h5" : "/param_vec.txt";
      params = wave_io.ParamFromFile(filename);
   }
   
   // Write parameter to ParaView output
   if (params && visualization)
   {
      cout << "Writing parameter field to ParaView output..." << endl;
      string param_dc_name = "Parameter";
      WaveVis param_vis(bottom_mesh, visualization, vis_steps, order_m, param_dc_name);
      
      ParaViewDataCollection *param_dc = param_vis.ParaviewDC();
      GridFunction *m_gf = nullptr;
      double t = 0;
      for (int k = 0; k < param_steps; k++)
      {
         if (k % 5 == 0)
         {
            m_gf = params[k];
            param_dc->RegisterField("parameter", m_gf);
            t = k*param_rate*dt;
            param_dc->SetCycle(k*param_rate);
            param_dc->SetTime(t * Cascadia::t0);
            param_dc->Save();
         }
      }
   }
   
   // 11a. Specify parameter (load) and call p2o map
   Vector **obs = nullptr;
   if (fwd == 1) // do one forward solve
   {
      if (WaveSolution::IsKnown())
      {
         // For known solution, parameters are already defined
         // by TD function WaveSolution::mParameter
         wave_p2o->GetObs(obs);
      }
      else
      {
         // For unknown solution, define parameters (input) for p2o map
         if (params)
         {
            // load defined by GridFunctionCoefficient for each time step
            wave_p2o->Mult(params, obs);
         }
         else
         {
            // load defined by a time-dependent function
            wave_p2o->Mult(WaveSolution::mParameter, obs);
         }
      }
      // 11b. Write observations to file
      if (observations)
      {
         wave_io.FwdToFile(obs);
      }
      
      bool test_misfit = false; // unit test
      if (test_misfit)
      {
//         double misfit = wave_p2o->EvalMisfit(obs, obs);
         double misfit = EvalMisfit(obs, obs, obs_steps);
         MFEM_VERIFY(abs(misfit) < 1.0e-14, "misfit > 0");
         
         Vector **noisy_obs = new Vector*[obs_steps];
         for (int k = 0; k < obs_steps; k++)
         {
            Vector &tmp = *(obs[k]);
            noisy_obs[k] = new Vector(tmp);
         }
         double rel_noise = 0.01;
//         double noise_var = wave_p2o->AddNoise(noisy_obs, rel_noise);
//         misfit = wave_p2o->EvalMisfit(obs, noisy_obs);
         double noise_var = AddNoise(noisy_obs, rel_noise, obs_steps);
         misfit = EvalMisfit(obs, noisy_obs, obs_steps);
         
         cout << endl << "rel noise = " << rel_noise << endl;
         cout         << "noise var = " << noise_var << endl;
         cout         << "   misfit = " << misfit << endl << endl;
         
         for (int k = 0; k < obs_steps; k++) { delete noisy_obs[k]; noisy_obs[k] = nullptr; }
         delete[] noisy_obs; noisy_obs = nullptr;
      }
      
      bool test_obs_io = false; // unit test
      if (test_obs_io)
      {
         string filename(output_dir);
         filename += (hdf) ? "/obs_vec.h5" : "/obs_vec.txt";
         
         // Write observations to file
         double rel_noise = 0.0;
         double noise_cov = 1.0;
         wave_io.ObsToFile(filename, obs, rel_noise, noise_cov);
         
         // Read observations from file
         Vector **obs_in = wave_io.ObsFromFile(filename);
         
         for (int k = 0; k < obs_steps; k++)
         {
            Vector &o_out = *obs[k];
            Vector &o_in = *obs_in[k];
            //cout << endl << "obs step k = " << k << endl;
            for (int i = 0; i < n_obs; i++)
            {
               //cout << "o_out[" << i << "] = " << o_out[i] << endl;
               //cout << "o_in [" << i << "] = " << o_in[i] << endl;
               double diff = abs(o_out[i]-o_in[i]);
               MFEM_VERIFY(diff < 1.0e-12 || diff < 1.0e-5*abs(o_out[i]),
                           "o_out != o_in");
            }
         }
         
         for (int k = 0; k < obs_steps; k++) { delete obs_in[k]; obs_in[k] = nullptr; }
         delete[] obs_in; obs_in = nullptr;
      }
   }
   else if (fwd == 2)
   {
      cout << "Doing repeated forward solves to create and export p2o map." << endl;

      // Observations are defined every m-th step of the parameter input
      // so the Block Toeplitz structure occurs at the level of m subblocks
      const int m = param_steps / obs_steps; // SWITCH: write compact p2o map (first column block)
//      const int m = param_steps; // write full p2o map (all column blocks)

      // Create data (load) for adjoint problem
      params = new GridFunction*[param_steps];
      for (int k = 0; k < param_steps; k++)
      {
         params[k] = new GridFunction(param_space);
         *(params[k]) = 0.0;
      }

      // Activate one nodal dof (one parameter) at a time
      // to obtain one column of p2o map
      chrono.Clear();
      for (int k = 0; k < m; k++) // do m := param_steps to write all columns (not compact version)
      {
         MFEM_VERIFY(params[k]->Size() == n_param, "Check vec size.");
         for (int j = 0; j < n_param; j++)
         {
            cout << endl << "+++ Forward solve number " << (k*n_param)+(j+1) << " / " << m*n_param << endl << endl;
            // Activate j-th parameter at k-th time (forward load),
            // yielding the (k*n_param+j)-th column of the p2o map
            (*(params[k]))[j] = 1.0;
            wave_p2o->Mult(params, obs);
            wave_io.FwdToFile(obs);
            wave_fwd->ResetLoad();
            (*(params[k]))[j] = 0.0;
         }
      }
      cout << endl << m*n_param << " forward solves took " << chrono.RealTime() << " seconds." << endl;
   }
   
   Vector **data = nullptr;
   GridFunction **adj_gf = nullptr;
   if (adj == 1) {
      // 12a. Specify data (adjoint load) and call transpose p2o map
      if (fwd == 1)
      {
         data = obs; // use observations from fwd solve as input to adjoint
      }
      else
      {
         data = new Vector*[obs_steps];
         for (int k = 0; k < obs_steps; k++)
         {
            data[k] = new Vector(n_obs);
            *(data[k]) = 0.0;
            (*(data[k]))[0] = 1.0; // TODO: define data (adjoint load)
         }
      }
   
      // load defined by GridFunctionCoefficient from data for each time step
      wave_p2o->MultTranspose(data, adj_gf);
      
      // 12b. Write adjoint output to file
      if (observations)
      {
         wave_io.AdjToFile(adj_gf);
      }
   }
   else if (adj == 2)
   {
      cout << "Doing repeated adjoint solves to create and export adjoint p2o map." << endl;

      // Parameters are defined at every observation time
      const int m = 1; // SWITCH: write compact adjoint p2o map (last column block)
//      const int m = obs_steps; // write full adjoint p2o map (all column blocks)

      // Create data (load) for adjoint problem
      data = new Vector*[obs_steps];
      for (int k = 0; k < obs_steps; k++)
      {
         data[k] = new Vector(n_obs);
         *(data[k]) = 0.0;
      }

      // Activate one nodal dof (one sensor) at a time
      // to obtain one column of adjoint p2o map
      chrono.Clear();
      int first_vec = 0;
      int last_vec = n_obs;
      if (adjvec >= 0)
      {
         first_vec = adjvec;
         last_vec = adjvec+1;
         if (last_vec > n_obs)
         {
            MFEM_WARNING("adjvec >= n_obs! Setting to n_obs-1.");
            last_vec = n_obs;
         }
      }
      for (int k = obs_steps-m; k < obs_steps; k++) // do m := obs_steps to write all columns (not compact version)
      {
         for (int j = first_vec; j < last_vec; j++)
         {
            cout << endl << "+++ Adjoint solve number " << (k-(obs_steps-m))*n_obs+(j+1) << " / " << m*n_obs << endl << endl;
            // Compacted version: activate j-th sensor at terminal time (adjoint load),
            // yielding the j-th column of the last column block of adjoint p2o map
            (*(data[k]))[j] = 1.0; // TODO: testing
            wave_p2o->MultTranspose(data, adj_gf);
            wave_io.AdjToFile(adj_gf, j);
            wave_adj->ResetLoad();
            (*(data[k]))[j] = 0.0;
         }
      }
      cout << endl << m*(last_vec-first_vec) << " adjoint solves took " << chrono.RealTime() << " seconds." << endl;
   }
   
   // 13. Prior
   int prior_type = prior%10;
   MFEM_VERIFY(prior_type <= 1, "Bi-Laplacian Prior not yet implemented.");
   if (prior)
   {
      height = n_param * param_steps;
      WavePrior wave_prior(*param_space,
                           height, prior_type,
                           dt, param_rate, param_steps,
                           alpha1, alpha2, alpha3);
      
      // Create global prior CSR matrix and write to file
      if (prior > 10)
      {
         wave_prior.PriorToFile(hdf);
      }
      
      // Compute regularization part of cost functional
      Vector param_vec(height);
      if (params)
      {
         // Use previously defined parameter GridFunctions
         int offset = 0;
         for (int k = 0; k < param_steps; k++)
         {
            param_vec.SetVector(*(params[k]), offset);
            offset += n_param;
         }
      }
      else
      {
         // Create parameter vector
         FunctionCoefficient param_coeff(WaveSolution::mParameter);
         GridFunction param_gf(param_space);
         double t = 0;
         int offset = 0;
         for (int k = 0; k < param_steps; k++)
         {
            t = k*param_rate*dt;
            param_coeff.SetTime(t);
            param_gf.ProjectCoefficient(param_coeff);
            param_vec.SetVector(param_gf, offset);
            offset += n_param;
         }
      }
      Vector reg_vec(height);
      wave_prior.Mult(param_vec, reg_vec);
      
      double reg_val = InnerProduct(param_vec, reg_vec);
      
      wave_prior.MultReg1(param_vec, reg_vec);
      double r1 = InnerProduct(param_vec, reg_vec);
      
      wave_prior.MultReg2(param_vec, reg_vec);
      double r2 = InnerProduct(param_vec, reg_vec);
      
      wave_prior.MultReg3(param_vec, reg_vec);
      double r3 = InnerProduct(param_vec, reg_vec);
      
      double reg_sum = r1 + r2 + r3;
      
      cout << endl << "Regularization part of cost functional:" << endl;
      cout << "  alpha1 = " << alpha1  << endl;
      cout << "  alpha2 = " << alpha2  << endl;
      cout << "  alpha3 = " << alpha3  << endl;
      cout << " reg_val = " << reg_val << endl;
      cout << "      r1 = " <<      r1 << endl;
      cout << "      r2 = " <<      r2 << endl;
      cout << "      r3 = " <<      r3 << endl;
      cout << " reg_sum = " << reg_sum << endl;
      
      // TESTING PARAM I/O
      bool test_param_io = false; // unit test
      if (test_param_io)
      {
         // Make sure params has been initialized as GridFunction**
         if (!params)
         {
            params = new GridFunction*[param_steps];
            for (int k = 0; k < param_steps; k++)
            {
               params[k] = new GridFunction(param_space);
               *(params[k]) = 0.0;
            }
            int offset = 0;
            for (int k = 0; k < param_steps; k++)
            {
               params[k]->MakeRef(param_vec, offset, n_param);
               offset += n_param;
            }
         }
         
         WaveIO::CreateDirectory(output_dir, nullptr, 0);
         string filename(output_dir);
         filename += (hdf) ? "/param_vec.h5" : "/param_vec.txt";
         
         // Write params to file
         wave_io.ParamToFile(filename, params);
         
         // Read params from file
         GridFunction **params_in = wave_io.ParamFromFile(filename);
         
         for (int k = 0; k < param_steps; k++)
         {
            Vector &p_out = *params[k];
            Vector &p_in = *params_in[k];
            //cout << endl << "param step k = " << k << endl;
            for (int i = 0; i < n_param; i++)
            {
               //cout << "params_out[" << i << "] = " << p_out[i] << endl;
               //cout << "params_in [" << i << "] = " << p_in[i] << endl;
               double diff = abs(p_out[i]-p_in[i]);
               MFEM_VERIFY(diff < 1.0e-12 || diff < 1.0e-5*abs(p_out[i]),
                           "p_out != p_in");
            }
         }
         
         for (int k = 0; k < param_steps; k++) { delete params_in[k]; params_in[k] = nullptr; }
         delete[] params_in; params_in = nullptr;
      }
   }
   
   // 14. Create a noisy data vector from observations
   bool noise_obs = true;
   if (noise_obs)
   {
      if (!obs)
      {
         // Read observations from file
         string filename(output_dir);
         filename += (hdf) ? "/obs_vec.h5" : "/obs_vec.txt";
         obs = wave_io.ObsFromFile(filename);
      }
//      double misfit = wave_p2o->EvalMisfit(obs, obs);
      double misfit = EvalMisfit(obs, obs, obs_steps);
      MFEM_VERIFY(abs(misfit) < 1.0e-14, "misfit > 0");
      
      Vector **noisy_obs = new Vector*[obs_steps];
      for (int k = 0; k < obs_steps; k++)
      {
         Vector &tmp = *(obs[k]);
         noisy_obs[k] = new Vector(tmp);
      }
      double rel_noise = 0.01;
//      double noise_var = wave_p2o->AddNoise(noisy_obs, rel_noise);
//      misfit = wave_p2o->EvalMisfit(obs, noisy_obs);
      double noise_var = AddNoise(noisy_obs, rel_noise, obs_steps);
      misfit = EvalMisfit(obs, noisy_obs, obs_steps);
      
      cout << endl << "rel noise = " << rel_noise << endl;
      cout         << "noise var = " << noise_var << endl;
      cout         << "   misfit = " << misfit << endl << endl;
      
      // Write observations to file
      string filename(output_dir);
      filename += (hdf) ? "/noisy_obs_vec.h5" : "/noisy_obs_vec.txt";
//      wave_io.ObsToFile(filename, noisy_obs, rel_noise, noise_var);
      wave_io.ObsToFile(filename, noisy_obs, rel_noise, 1.0);
      
      for (int k = 0; k < obs_steps; k++) { delete noisy_obs[k]; noisy_obs[k] = nullptr; }
      delete[] noisy_obs; noisy_obs = nullptr;
   }

   // 15. Free the used memory.
   cout << endl << "cascadia_seq: freeing memory" << endl;
   if (params)
   {
      for (int k = 0; k < param_steps; k++) { delete params[k]; params[k] = nullptr; }
      delete[] params; params = nullptr;
   }
   if (obs)
   {
      for (int k = 0; k < obs_steps; k++) { delete obs[k]; obs[k] = nullptr; }
      if (data==obs) { data = nullptr; }
      delete[] obs; obs = nullptr;
   }
   if (data)
   {
      for (int k = 0; k < obs_steps; k++) { delete data[k]; data[k] = nullptr; }
      delete[] data; data = nullptr;
   }
   if (adj_gf)
   {
      for (int k = 0; k < param_steps; k++) { delete adj_gf[k]; adj_gf[k] = nullptr; }
      delete[] adj_gf; adj_gf = nullptr;
   }
   
   delete wave_p2o;
   delete wave_fwd;
   delete wave_adj;
   delete ode_solver;
   delete param_space;
   delete param_h1_coll;
   delete bottom_mesh;
   delete DG_space;
   delete CG_space;
   delete l2_coll;
   delete h1_coll;
   delete mesh;

   cout << "cascadia_seq: all done." << endl << endl;
   return 0;
}


/// Evaluate misfit part of cost functional (between model observations and data)
double EvalMisfit(Vector** &obs, Vector** &data, int obs_steps)
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
double AddNoise(Vector** &obs, double rel_noise, int obs_steps)
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
      for (int i = 0; i < obs_vec.Size(); i++)
      {
         obs_vec[i] += distr(rnd_gen);
      }
   }
   
   return noise_std_dev * noise_std_dev;
}
