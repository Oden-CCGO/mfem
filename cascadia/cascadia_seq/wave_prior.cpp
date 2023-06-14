
#include "cascadia.hpp"
#include "hdf5.h"

namespace mfem
{

using namespace std;

WavePrior::WavePrior(FiniteElementSpace &fes_m_,
                     int height_, int type_,
                     double time_step_, int param_rate_, int param_steps_,
                     double alpha1_, double alpha2_, double alpha3_)
   : Operator(height_, 0.0),
     fes_m(fes_m_),
     M(nullptr), m_var(nullptr), K(nullptr), k_var(nullptr),
     D(nullptr), C(nullptr), L(nullptr), B(nullptr), block_offsets(param_steps_+1),
     alpha1(nullptr), alpha2(nullptr), alpha3(nullptr),
     type(type_), time_step(time_step_), dt(time_step_*param_rate_),
     param_rate(param_rate_), param_steps(param_steps_)
{
   StopWatch chrono;
   chrono.Start();
   
   // Preliminary
   MFEM_VERIFY(type==1, "Only Laplacian prior is implemented.");
   
   MFEM_VERIFY(alpha1_>=0.0, "alpha1 < 0.");
   MFEM_VERIFY(alpha2_>=0.0, "alpha2 < 0.");
   MFEM_VERIFY(alpha3_>=0.0, "alpha3 < 0.");
   
   alpha1 = new ConstantCoefficient(alpha1_);
   alpha2 = new ConstantCoefficient(alpha2_);
   alpha3 = new ConstantCoefficient(alpha3_);
   
   // 1. Assemble FE mass matrix: alpha1(m,v)
   m_var = new BilinearForm(&fes_m);
   m_var->AddDomainIntegrator(new MassIntegrator(*alpha1));
   m_var->Assemble();
   m_var->Finalize();
   M = &(m_var->SpMat());
//   M->PrintInfo(mfem::out);
   
   cout << endl << "Timer: Mass       matrix assembly (prior): "
        << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();
   
   // 2. Assemble FE stiffness matrix: alpha2(grad m, grad v)
   k_var = new BilinearForm(&fes_m);
   k_var->AddDomainIntegrator(new DiffusionIntegrator(*alpha2));
   k_var->Assemble();
   k_var->Finalize();
   K = &(k_var->SpMat());
//   K->PrintInfo(mfem::out);
   
   if (M->NumNonZeroElems() != K->NumNonZeroElems())
   {
      MFEM_WARNING("  M->NumNonZeroElems() = " << M->NumNonZeroElems() <<
                   ", K->NumNonZeroElems() = " << K->NumNonZeroElems());
   }
   
   cout << "Timer: Stiffness  matrix assembly (prior): "
        << chrono.RealTime() << " seconds." << endl;
   chrono.Clear();

   // 3. Assemble space-time prior (reg matrix)
   //    R is Block-Diagonal if and only if alpha3 = 0
   //      each block on the diagonal is the same
   //      R = | M+K  0   0   0   |
   //          | 0    M+K 0   0   |
   //          | 0    0   M+K 0   |
   //          | 0    0   0   M+K |
   
   //    If alpha3 > 0, then R is Block-Tri-Diagonal
   //      R = | C   L   0   0 |
   //          | L   D   L   0 |
   //          | 0   L   D   L |
   //          | 0   0   L   C |
   //      where D := M+K+C
   
   // 3a. create block offset array
   const int nr_params = fes_m.GetVSize();
   int offset = 0;
   for (int i = 0; i <= param_steps; i++)
   {
      block_offsets[i] = offset;
      offset += nr_params;
   }
   
   // 3b. create abstract block matrix and set blocks
   B = new BlockMatrix(block_offsets);
   
   if (alpha3_ > 0)
   {
      C = new SparseMatrix(*M);
      *C *= 2*(alpha3_/alpha1_)/(dt*dt);
   
      L = new SparseMatrix(*M);
      *L *= -(alpha3_/alpha1_)/(dt*dt);
      
      D = new SparseMatrix(*M);
      *D += *K;
      *D += *C;
      
      for (int i = 1; i < param_steps-1; i++)
      {
         B->SetBlock(i, i  , D); // diag block
         B->SetBlock(i, i-1, L); // left block
         B->SetBlock(i, i+1, L); // right block
      }
      // Set first row block (behaves as if a zero Dir.BC is to the left)
      B->SetBlock(0, 0, C);
      B->SetBlock(0, 1, L);
      
      // Set last row block (behaves as if a zero Dir.BC is to the right)
      int i = param_steps-1;
      B->SetBlock(i, i  , C);
      B->SetBlock(i, i-1, L);
   }
   else
   {
      D = new SparseMatrix(*M);
      *D += *K;
      
      for (int i = 0; i < param_steps; i++)
      {
         // set diagonal blocks
         B->SetBlock(i, i, D);
      }
   }

   cout << "Timer: Space-time matrix assembly (prior): "
        << chrono.RealTime() << " seconds." << endl;
   
   // Specify whether PriorToFile re-indexes CSR matrix before
   // writing to file, as is needed by the FFT matvec code
   reindex = true;
}

void WavePrior::Mult(const Vector &x, Vector &y) const
{
   B->Mult(x, y);
}

void WavePrior::MultReg1(const Vector &x, Vector &y) const
{
   if (!(alpha1->constant > 0))
   {
      MFEM_WARNING("! (alpha1 > 0)");
      y = 0.0;
      return;
   }
   
   BlockMatrix A(block_offsets);
   for (int i = 1; i < param_steps-1; i++)
   {
      A.SetBlock(i, i, M);
   }
   
   if (!(alpha3->constant > 0))
   {
      A.SetBlock(0, 0, M);
      int i = param_steps-1;
      A.SetBlock(i, i, M);
   }
   
   A.Mult(x, y);
}

void WavePrior::MultReg2(const Vector &x, Vector &y) const
{
   if (!(alpha2->constant > 0))
   {
      MFEM_WARNING("! (alpha2 > 0)");
      y = 0.0;
      return;
   }
   
   BlockMatrix A(block_offsets);
   for (int i = 1; i < param_steps-1; i++)
   {
      A.SetBlock(i, i, K);
   }
   
   if (!(alpha3->constant > 0))
   {
      A.SetBlock(0, 0, K);
      int i = param_steps-1;
      A.SetBlock(i, i, K);
   }
   
   A.Mult(x, y);
}

void WavePrior::MultReg3(const Vector &x, Vector &y) const
{
   if (!(alpha3->constant > 0))
   {
      MFEM_WARNING("! (alpha3 > 0)");
      y = 0.0;
      return;
   }
   
   BlockMatrix A(block_offsets);
   for (int i = 1; i < param_steps-1; i++)
   {
      A.SetBlock(i, i  , C); // diag block
      A.SetBlock(i, i-1, L); // left block
      A.SetBlock(i, i+1, L); // right block
   }
   // Set first row block (behaves as if a zero Dir.BC is to the left)
   A.SetBlock(0, 0, C);
   A.SetBlock(0, 1, L);
   
   // Set last row block (behaves as if a zero Dir.BC is to the right)
   int i = param_steps-1;
   A.SetBlock(i, i  , C);
   A.SetBlock(i, i-1, L);
   
   A.Mult(x, y);
}

SparseMatrix* WavePrior::ReindexCSR(const SparseMatrix *R) const
{
   const int size = R->Height();
   const int nnz = R->NumNonZeroElems();
   const int *R_I = R->GetI();
   const int *R_J = R->GetJ();
   const double *R_A = R->GetData();
   
   const int n_param = size/param_steps;
   MFEM_VERIFY(size%param_steps==0,
               "param_steps does not evenly divide matrix size.");
   
   // The original matrix R has the dof ordering:
   // time steps (outer) -> space dofs (inner)
   // number of blocks   = param_steps x param_steps
   // size of each block = n_param x n_param
   
   // The new matrix has the dof ordering:
   // space dofs (outer) -> time steps (inner)
   // number of blocks   = n_param x n_param
   // size of each block = param_steps x param_steps
   
   int *I = new int[size+1]; I[0] = 0;
   int *J = new int[nnz];
   double *A = new double[nnz];
   
   for (int i = 0; i < size; i++)
   {
      int r_i = (i % param_steps) * n_param + i / param_steps;
      // start test
      int ii = (r_i % n_param) * param_steps + r_i / n_param;
      MFEM_VERIFY(i == ii, "i != ii"); // TODO: remove
      // end test
      int row_nnz = R_I[r_i+1] - R_I[r_i];
      I[i+1] = I[i] + row_nnz;
      
      for (int k = 0; k < row_nnz; k++)
      {
         int r_l = R_I[r_i] + k;
         int r_j = R_J[r_l];
         int l = I[i] + k;
         J[l] = (r_j % n_param) * param_steps + r_j / n_param;
         A[l] = R_A[r_l];
      }
   }
      
   return new SparseMatrix(I, J, A, size, size);
}

void WavePrior::PriorToFile(bool binary)
{
   StopWatch chrono;
   chrono.Start();

   // Get global sparse matrix from an abstract BlockMatrix
   // note: copies memory internally
   chrono.Clear();
   SparseMatrix *R = B->CreateMonolithic();
   cout << "Timer: Creating csr matrix (prior)       : " << chrono.RealTime() << " seconds." << endl;
   
   // Re-index csr matrix (ordering of dofs)
   if (reindex)
   {
      chrono.Clear();
      SparseMatrix *S = ReindexCSR(R);
      delete R;
      R = S;
      cout << "Timer: Re-indexing csr matrix (prior)    : " << chrono.RealTime() << " seconds." << endl;
   }
   
   // Sort column indices (may be needed by FFT matvec code)
   if (!(R->ColumnsAreSorted()))
   {
      chrono.Clear();
      R->SortColumnIndices();
      cout << "Timer: Sorting csr matrix columns (prior): " << chrono.RealTime() << " seconds." << endl;
   }
   
   if (R->IsSymmetric() > 0.0) { MFEM_WARNING("R Matrix is not sym!"); }
   
   // Create directory
   string prefix, suffix;
   
   string rel_path = WaveIO::output_dir + "/prior/";
   string filename = "csr";
   if (binary) { rel_path += "binary"; suffix = ".h5"; }
   else        { rel_path += "text"; suffix = ".txt"; }
   
   WaveIO::CreateDirectory(rel_path, nullptr, 0);
   
   // Write meta file with info about prior
   string meta_filename = rel_path + "/meta_prior";
   cout << endl << "PriorToFile: writing to " << meta_filename << endl;
   ofstream meta_file;
   meta_file.open(meta_filename);
   meta_file << param_steps << endl;          // <int> number of blocks
   meta_file << R->Height() << endl;          // <int> size of the matrix
   meta_file << R->NumNonZeroElems() << endl; // <int> number of non-zero entries
   meta_file << type << endl;                 // type of prior (1: Laplacian; 2 Bi-Laplacian)
   meta_file << alpha1->constant << endl;     // <double> (regularization parameter for |m|)
   meta_file << alpha2->constant << endl;     // <double> (regularization parameter for |grad m|)
   meta_file << alpha3->constant << endl;     // <double> (regularization parameter for |dm/dt|)
   meta_file << filename << endl;             // <string> filename for csr matrix file
   meta_file << suffix << endl;               // <string> suffix for csr matrix file (.txt or .h5)
   meta_file << reindex << endl;              // <bool> specifies whether csr matrix was re-indexed
   meta_file.close();
   cout << "PriorToFile: done." << endl << endl;
   
   R->PrintInfo(mfem::out); // TODO: remove verbose
   
   // Write prior (csr matrix) to file
   chrono.Clear();
   
   string filename_csr = rel_path+"/"+filename+suffix;
   cout << "PriorToFile: writing to " << filename_csr << endl;
   
   if (binary)
   {
      // Create file
      hid_t file_id;
      herr_t status;
      file_id = H5Fcreate(filename_csr.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      // Create dataspaces
      hsize_t rank = 1;
      hsize_t I_dim = R->Height()+1;
      hsize_t J_dim = R->NumNonZeroElems();
      hsize_t A_dim = R->NumNonZeroElems();
      hid_t I_dspace_id = H5Screate_simple(rank, &I_dim, NULL);
      hid_t J_dspace_id = H5Screate_simple(rank, &J_dim, NULL);
      hid_t A_dspace_id = H5Screate_simple(rank, &A_dim, NULL);

      // Create datasets
      hid_t I_dset_id;
      const char* I_dset_name = "indptr";
      I_dset_id = H5Dcreate(file_id, I_dset_name, H5T_NATIVE_INT, I_dspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      hid_t J_dset_id;
      const char* J_dset_name = "indices";
      J_dset_id = H5Dcreate(file_id, J_dset_name, H5T_NATIVE_INT, J_dspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      hid_t A_dset_id;
      const char* A_dset_name = "data";
      A_dset_id = H5Dcreate(file_id, A_dset_name, H5T_NATIVE_DOUBLE, A_dspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write datasets
      const int *I = R->GetI();
      const int *J = R->GetJ();
      const double *A = R->GetData();
      status = H5Dwrite(I_dset_id, H5T_NATIVE_INT   , H5S_ALL, H5S_ALL, H5P_DEFAULT, I);
      status = H5Dwrite(J_dset_id, H5T_NATIVE_INT   , H5S_ALL, H5S_ALL, H5P_DEFAULT, J);
      status = H5Dwrite(A_dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, A);

      // Attach attributes to datasets
      const int size = R->Height();
      hid_t attr_dspace_id = H5Screate(H5S_SCALAR);
      hid_t attr_id = H5Acreate(I_dset_id, "size", H5T_NATIVE_INT,
                                attr_dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, H5T_NATIVE_INT, &size);
      
      status = H5Aclose(attr_id);
      status = H5Sclose(attr_dspace_id);

      // Close datasets
      status = H5Dclose(I_dset_id);
      status = H5Dclose(J_dset_id);
      status = H5Dclose(A_dset_id);

      // Close dataspaces
      status = H5Sclose(I_dspace_id);
      status = H5Sclose(J_dspace_id);
      status = H5Sclose(A_dspace_id);

      // Close file
      status = H5Fclose(file_id);
      if (status < 0) { MFEM_WARNING("H5Fclose < 0"); }
   }
   else
   {
      ofstream csr_file;
      csr_file.open(filename_csr);
      R->PrintCSR2(csr_file);
      csr_file.close();
   }
   
   cout << "PriorToFile: done." << endl << endl;

   cout << "Timer: Space-time csr output (prior)     : " << chrono.RealTime() << " seconds." << endl;
   
   delete R;
}

WavePrior::~WavePrior()
{
   delete m_var;
   delete k_var;
   
   delete C;
   delete L;
   delete D;
   delete B;
   
   delete alpha1;
   delete alpha2;
   delete alpha3;
}

}
