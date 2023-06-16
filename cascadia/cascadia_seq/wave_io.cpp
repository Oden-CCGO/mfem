
#include "cascadia.hpp"
#include "hdf5.h"

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

namespace mfem
{

using namespace std;

string WaveIO::output_dir = ".";
bool WaveIO::adj_reverse_order = false;

WaveIO::WaveIO(FiniteElementSpace* param_space_, WaveObservationOp &wave_obs_,
               int n_steps_, double dt_, int n_param_, int n_obs_,
               int param_rate_, int obs_rate_, bool binary_)
   : param_space(param_space_), wave_obs(wave_obs_),
     n_steps(n_steps_), dt(dt_), n_param(n_param_), n_obs(n_obs_),
     param_rate(param_rate_), obs_rate(obs_rate_), binary(binary_)
{
   param_steps = n_steps/param_rate;
   obs_steps = n_steps/obs_rate;
   
   count_fwd_text = 0;
   count_fwd_binary = 0;
   
   init_meta_fwd = false;
   init_meta_adj = false;
   
   memcpy = false;
}

void WaveIO::MetaToFile(bool adj)
{
   if (adj)
   {
      if (init_meta_adj) { return; }
      init_meta_adj = true;
   }
   else
   {
      if (init_meta_fwd) { return; }
      init_meta_fwd = true;
   }

   string prefix, suffix;
   
   rel_path = output_dir + "/p2o/";
   if (binary) { rel_path += "binary/"; suffix = ".h5"; }
   else        { rel_path += "text/"; suffix = ".txt"; }
   
   string filename = rel_path;
   if (adj)
   {
      filename += "meta_adj";
      prefix_adj = "adj/vec_";
      prefix = prefix_adj;
      CreateDirectory(rel_path+"adj", nullptr, 0);
   }
   else
   {
      filename += "meta_fwd";
      prefix_fwd = "fwd/vec_";
      prefix = prefix_fwd;
      CreateDirectory(rel_path+"fwd", nullptr, 0);
   }
   
   // Note: the output here defines the block structure of the
   //       block Toeplitz system used in the FFT matvec algorithm
   MFEM_VERIFY(param_steps % obs_steps == 0,
               "Block Toeplitz assumes param_rate evenly divides obs_rate.");
   // Observations are defined every m-th step of the parameter input
   int m = param_steps / obs_steps;
   int params_per_block = n_param * m;
   
   cout << "MetaToFile: writing to " << filename << endl;
   ofstream meta_file;
   meta_file.open(filename);
   if (meta_file.is_open())
   {
      meta_file << n_obs << endl;            // <int> number of observations per block
      meta_file << params_per_block << endl; // <int> number of parameters per block
      meta_file << obs_steps << endl;        // <int> number of blocks
      meta_file << prefix << endl;        // <string> prefix for vector files
      meta_file << suffix << endl;        // <string> suffix for vector files (.txt or .h5)
      if (adj)
      {
         meta_file << adj_reverse_order << endl; // <bool> adj_vec written in block-reverse order
      }
      else
      {
         const int nr_sensors = wave_obs.GetNrSensors();
         const DenseMatrix *sensors_pts = wave_obs.GetSensorCoords();
         for (int i = 0; i < nr_sensors; i++)
         {
            meta_file << (*sensors_pts)(0,i) << endl;
            meta_file << (*sensors_pts)(1,i) << endl;
         }
      }
      meta_file.close();
   }
   else
   {
      MFEM_WARNING("Could not open meta_file.");
   }
   cout << "MetaToFile: done." << endl;
}

void WaveIO::FwdToFile(Vector **obs)
{
   ostringstream filename_oss;
   if (binary)
   {
      MetaToFile(false);
      filename_oss << rel_path << prefix_fwd
                   << setfill('0') << setw(6) << count_fwd_binary << ".h5";
      count_fwd_binary++;
   }
   else
   {
      MetaToFile(false);
      filename_oss << rel_path << prefix_fwd
                   << setfill('0') << setw(6) << count_fwd_text << ".txt";
      count_fwd_text++;
   }

   cout << "FwdToFile: writing to " << filename_oss.str() << endl;
   ObsToFile(filename_oss.str(), obs);
   cout << "FwdToFile: done." << endl;
}

void WaveIO::AdjToFile(GridFunction **adj, int adjvec)
{
   ostringstream filename_oss;
   if (binary)
   {
      MetaToFile(true);
      filename_oss << rel_path << prefix_adj
                   << setfill('0') << setw(6) << adjvec << ".h5";
   }
   else
   {
      MetaToFile(true);
      filename_oss << rel_path << prefix_adj
                   << setfill('0') << setw(6) << adjvec << ".txt";
   }

   cout << "AdjToFile: writing to " << filename_oss.str() << endl;
   ParamToFile(filename_oss.str(), adj);
   cout << "AdjToFile: done." << endl;
}

void WaveIO::ObsToFile(const std::string &filename, Vector **obs,
                       double rel_noise, double noise_cov)
{
   if (binary)
   {
      // Create file
      hid_t file_id;
      herr_t status;
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // Create dataspace
      hsize_t rank = 1;
      hsize_t ddim = obs_steps*n_obs;
      hid_t dspace_id;
      dspace_id = H5Screate_simple(rank, &ddim, NULL);

      // Create dataset
      hid_t dset_id;
      dset_id = H5Dcreate(file_id, "vec", H5T_NATIVE_DOUBLE, dspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Attach attributes to datasets
      hid_t attr_dspace_id = H5Screate(H5S_SCALAR);
      hid_t attr_id;
      
      // - Attribute 1: obs_steps
      attr_id = H5Acreate(dset_id, "obs_steps", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &obs_steps);
      status = H5Aclose(attr_id);
      
      // - Attribute 2: dt
      attr_id = H5Acreate(dset_id, "dt", H5T_NATIVE_DOUBLE,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &dt);
      status = H5Aclose(attr_id);
      
      // - Attribute 3: obs_rate
      attr_id = H5Acreate(dset_id, "obs_rate", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &obs_rate);
      status = H5Aclose(attr_id);
      
      // - Attribute 4: n_obs
      attr_id = H5Acreate(dset_id, "n_obs", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &n_obs);
      status = H5Aclose(attr_id);
      
      // - Attribute 5: rel_noise
      attr_id = H5Acreate(dset_id, "rel_noise", H5T_NATIVE_DOUBLE,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &rel_noise);
      
      // - Attribute 6: noise_cov
      attr_id = H5Acreate(dset_id, "noise_cov", H5T_NATIVE_DOUBLE,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &noise_cov);
      
      // - Attribute 7: reindex
      int reindex = 0;
      attr_id = H5Acreate(dset_id, "reindex", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &reindex);
      status = H5Aclose(attr_id);
      
      status = H5Sclose(attr_dspace_id);

      // Write dataset
      if (memcpy)
      {
         // Copying sub-vectors into global vector; then write to dataset
         double *dset_data = new double[obs_steps*n_obs];
         for (int k = 0; k < obs_steps; k++)
         {
            Vector &tmp = *(obs[k]);
            MFEM_VERIFY(tmp.Size() == n_obs,
                        "Size of obs vector does not match.");
            for (int i = 0; i < n_obs; i++)
            {
               dset_data[k*n_obs+i] = tmp[i];
            }
         }
         status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
      }
      else
      {
         // Writing sub-vectors directly to the dataset one-by-one
         hsize_t offset = 0;
         hsize_t count = n_obs;
         hsize_t stride = 1;
         hsize_t block = 1;
         
         // Select memory dataspace (same in each iteration)
         hid_t memspace_id = H5Screate_simple(rank, &count, NULL);
         
         for (int k = 0; k < obs_steps; k++)
         {
            MFEM_VERIFY(obs[k]->Size() == n_obs,
                        "Size of obs vector does not match.");
            const double *dset_data = obs[k]->HostRead();
            
            // Select file dataspace
            status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset, &stride, &count, &block);
            hssize_t numel = H5Sget_select_npoints(dspace_id);
            MFEM_VERIFY(numel == n_obs, "Inconsistent number of elements");
            
            // Write sub-vector to dataset
            status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, dset_data);
            
            offset += count;
         }
      }

      // Close dataset
      status = H5Dclose(dset_id);

      // Close dataspace
      status = H5Sclose(dspace_id);

      // Close file
      status = H5Fclose(file_id);
      if (status < 0) { MFEM_WARNING("H5Fclose < 0"); }
   }
   else
   {
      ofstream obs_file;
      obs_file.open(filename);
      if (obs_file.is_open())
      {
         for (int k = 0; k < obs_steps; k++)
         {
            Vector &tmp = *(obs[k]);
            MFEM_VERIFY(tmp.Size() == n_obs,
                        "Size of obs vector does not match.");
            for (int i = 0; i < n_obs; i++)
            {
               obs_file << setprecision(10) << scientific << tmp[i] << endl;
            }
         }
         obs_file.close();
      }
      else
      {
         MFEM_WARNING("Could not open obs_file.");
      }
   }
   
   cout << "ObsToFile: done." << endl;
}

Vector** WaveIO::ObsFromFile(const std::string &filename)
{
   Vector **obs = new Vector*[obs_steps];
   if (binary)
   {
      // Open existing file
      hid_t file_id;
      herr_t status;
      file_id = H5Fopen(filename.c_str(), H5P_DEFAULT, H5P_DEFAULT);

      // Open existing dataset
      hid_t dset_id;
      dset_id = H5Dopen(file_id, "/vec", H5P_DEFAULT);

      // Read dataset
      if (memcpy)
      {
         // Read global vector from dataset; then, copy into sub-vectors
         double *dset_data = new double[obs_steps*n_obs];
         
         status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
         
         for (int k = 0; k < obs_steps; k++)
         {
            obs[k] = new Vector(n_obs);
            Vector &tmp = *(obs[k]);
            
            for (int i = 0; i < n_obs; i++)
            {
               tmp[i] = dset_data[k*n_obs+i];
            }
         }
      }
      else
      {
         // Reading sub-vectors directly from the dataset one-by-one
         hsize_t offset = 0;
         hsize_t count = n_obs;
         hsize_t stride = 1;
         hsize_t block = 1;
         
         // Create dataspace
         hsize_t rank = 1;
         hsize_t ddim = obs_steps*n_obs;
         hid_t dspace_id;
         dspace_id = H5Screate_simple(rank, &ddim, NULL);
         
         // Select memory dataspace (same in each iteration)
         hid_t memspace_id = H5Screate_simple(rank, &count, NULL);
         
         for (int k = 0; k < obs_steps; k++)
         {
            obs[k] = new Vector(n_obs);
            Vector &tmp = *(obs[k]);
            
            // Select file dataspace
            status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset, &stride, &count, &block);
            hssize_t numel = H5Sget_select_npoints(dspace_id);
            MFEM_VERIFY(numel == n_obs, "Inconsistent number of elements");
            
            // Read sub-vector from dataset
            double *dset_data = tmp.HostWrite();
            status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, dset_data);
            
            offset += count;
         }
         
         // Close dataspace
         status = H5Sclose(dspace_id);
      }

      // Close dataset
      status = H5Dclose(dset_id);

      // Close file
      status = H5Fclose(file_id);
      if (status < 0) { MFEM_WARNING("H5Fclose < 0"); }
   }
   else
   {
      ifstream obs_file;
      obs_file.open(filename);
      if (obs_file.is_open())
      {
         for (int k = 0; k < obs_steps; k++)
         {
            obs[k] = new Vector(n_obs);
            Vector &tmp = *(obs[k]);
            
            for (int i = 0; i < n_obs; i++)
            {
               obs_file >> tmp[i];
            }
         }
      }
      else
      {
         MFEM_WARNING("Could not open obs_file.");
      }
      obs_file.close();
   }
   
   cout << "ObsFromFile: done." << endl;
   
   return obs;
}

void WaveIO::ParamToFile(const std::string &filename, GridFunction **param)
{
   if (binary)
   {
      // Create file
      hid_t file_id;
      herr_t status;
      file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // Create dataspace
      hsize_t rank = 1;
      hsize_t ddim = param_steps*n_param;
      hid_t dspace_id;
      dspace_id = H5Screate_simple(rank, &ddim, NULL);

      // Create dataset
      hid_t dset_id;
      dset_id = H5Dcreate(file_id, "vec", H5T_NATIVE_DOUBLE, dspace_id,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Attach attributes to datasets
      hid_t attr_dspace_id = H5Screate(H5S_SCALAR);
      hid_t attr_id;
      
      // - Attribute 1: param_steps
      attr_id = H5Acreate(dset_id, "param_steps", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &param_steps);
      status = H5Aclose(attr_id);
      
      // - Attribute 2: n_param
      attr_id = H5Acreate(dset_id, "n_param", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &n_param);
      status = H5Aclose(attr_id);
      
      // - Attribute 3: reindex
      int reindex = 0;
      attr_id = H5Acreate(dset_id, "reindex", H5T_NATIVE_INT,
                          attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &reindex);
      status = H5Aclose(attr_id);
      
      status = H5Sclose(attr_dspace_id);
      
      // Write dataset
      if (memcpy)
      {
         // Copying sub-vectors into global vector; then write to dataset
         double *dset_data = new double[param_steps*n_param];
         
         for (int k = 0; k < param_steps; k++)
         {
            Vector &tmp = *(param[k]);
            MFEM_VERIFY(tmp.Size() == n_param,
                        "Size of param vector does not match.");
            for (int i = 0; i < n_param; i++)
            {
               dset_data[k*n_param+i] = tmp[i];
            }
         }
         status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
      }
      else
      {
         // Writing sub-vectors directly to the dataset one-by-one
         hsize_t offset = 0;
         hsize_t count = n_param;
         hsize_t stride = 1;
         hsize_t block = 1;
         
         // Select memory dataspace (same in each iteration)
         hid_t memspace_id = H5Screate_simple(rank, &count, NULL);
         
         for (int k = 0; k < param_steps; k++)
         {
            MFEM_VERIFY(param[k]->Size() == n_param,
                        "Size of param vector does not match.");
            const double *dset_data = param[k]->HostRead();
            
            // Select file dataspace
            status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset, &stride, &count, &block);
            hssize_t numel = H5Sget_select_npoints(dspace_id);
            MFEM_VERIFY(numel == n_param, "Inconsistent number of elements");
            
            // Write sub-vector to dataset
            status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, dset_data);
            
            offset += count;
         }
      }

      // Close dataset
      status = H5Dclose(dset_id);

      // Close dataspace
      status = H5Sclose(dspace_id);

      // Close file
      status = H5Fclose(file_id);
      if (status < 0) { MFEM_WARNING("H5Fclose < 0"); }
   }
   else
   {
      ofstream param_file;
      param_file.open(filename);
      if (param_file.is_open())
      {
         for (int k = 0; k < param_steps; k++)
         {
            GridFunction &tmp = *(param[k]);
            MFEM_VERIFY(tmp.Size() == n_param,
                        "Size of param vector does not match.");
            for (int i = 0; i < n_param; i++)
            {
               param_file << setprecision(10) << scientific << tmp[i] << endl;
            }
         }
         param_file.close();
      }
      else
      {
         MFEM_WARNING("Could not open param_file.");
      }
   }
   
   cout << "ParamToFile: done." << endl;
}

GridFunction** WaveIO::ParamFromFile(const std::string &filename)
{
   GridFunction** params = new GridFunction*[param_steps];
   if (binary)
   {
      // Open existing file
      hid_t file_id;
      herr_t status;
      file_id = H5Fopen(filename.c_str(), H5P_DEFAULT, H5P_DEFAULT);

      // Open existing dataset
      hid_t dset_id;
      dset_id = H5Dopen(file_id, "/vec", H5P_DEFAULT);

      // Read dataset
      if (memcpy)
      {
         // Read global vector from dataset; then, copy into sub-vectors
         double *dset_data = new double[param_steps*n_param];
         
         status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
         
         for (int k = 0; k < param_steps; k++)
         {
            params[k] = new GridFunction(param_space);
            GridFunction &tmp = *(params[k]);
            
            for (int i = 0; i < n_param; i++)
            {
               tmp[i] = dset_data[k*n_param+i];
            }
         }
      }
      else
      {
         // Reading sub-vectors directly from the dataset one-by-one
         hsize_t offset = 0;
         hsize_t count = n_param;
         hsize_t stride = 1;
         hsize_t block = 1;
         
         // Create dataspace
         hsize_t rank = 1;
         hsize_t ddim = param_steps*n_param;
         hid_t dspace_id;
         dspace_id = H5Screate_simple(rank, &ddim, NULL);
         
         // Select memory dataspace (same in each iteration)
         hid_t memspace_id = H5Screate_simple(rank, &count, NULL);

         for (int k = 0; k < param_steps; k++)
         {
            params[k] = new GridFunction(param_space);
            GridFunction &tmp = *(params[k]);

            // Select file dataspace
            status = H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset, &stride, &count, &block);
            hssize_t numel = H5Sget_select_npoints(dspace_id);
            MFEM_VERIFY(numel == n_param, "Inconsistent number of elements");

            // Read sub-vector from dataset
            double *dset_data = tmp.HostWrite();
            status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dspace_id, H5P_DEFAULT, dset_data);

            offset += count;
         }
         
         // Close dataspace
         status = H5Sclose(dspace_id);
      }

      // Close dataset
      status = H5Dclose(dset_id);

      // Close file
      status = H5Fclose(file_id);
      if (status < 0) { MFEM_WARNING("H5Fclose < 0"); }
   }
   else
   {
      ifstream param_file;
      param_file.open(filename);
      if (param_file.is_open())
      {
         for (int k = 0; k < param_steps; k++)
         {
            params[k] = new GridFunction(param_space);
            GridFunction &tmp = *(params[k]);
            
            for (int i = 0; i < n_param; i++)
            {
               param_file >> tmp[i];
            }
         }
      }
      else
      {
         MFEM_WARNING("Could not open param_file.");
      }
      param_file.close();
   }
   
   cout << "ParamFromFile: done." << endl;
   
   return params;
}

/// static method (copied from DataCollection::create_directory)
int WaveIO::CreateDirectory(const std::string &dir_name,
                            const Mesh *mesh, int myid)
{
   // create directories recursively
   const char path_delim = '/';
   std::string::size_type pos = 0;
   int err_flag;
#ifdef MFEM_USE_MPI
   const ParMesh *pmesh = dynamic_cast<const ParMesh*>(mesh);
#endif

   do
   {
      pos = dir_name.find(path_delim, pos+1);
      std::string subdir = dir_name.substr(0, pos);

#ifndef MFEM_USE_MPI
      err_flag = mkdir(subdir.c_str(), 0777);
      err_flag = (err_flag && (errno != EEXIST)) ? 1 : 0;
#else
      if (myid == 0 || pmesh == NULL)
      {
         err_flag = mkdir(subdir.c_str(), 0777);
         err_flag = (err_flag && (errno != EEXIST)) ? 1 : 0;
      }
#endif
   }
   while ( pos != std::string::npos );

#ifdef MFEM_USE_MPI
   if (pmesh)
   {
      MPI_Bcast(&err_flag, 1, MPI_INT, 0, pmesh->GetComm());
   }
#endif

   return err_flag;
}

} // close namespace mfem
