#ifndef CASCADIA_SUM_OPERATOR
#define CASCADIA_SUM_OPERATOR

#include "mfem.hpp"

namespace mfem
{

/// Sum Operator: x -> A(x) + b B(x).
class SumOperator : public mfem::Operator
{
private:
   const mfem::Operator &A_;
   const mfem::Operator &B_;
   double b_;
   mutable mfem::Vector z_;

public:
   /// Create an operator which is the sum of A and B.
   explicit SumOperator(const mfem::Operator *A, const mfem::Operator *B, double b)
      : mfem::Operator(A->Height(), A->Width()), A_(*A), B_(*B), b_(b) {
      MFEM_VERIFY(A->Height()==B->Height(),
                 "SumOperator: A->Height()!=B->Height()");
      MFEM_VERIFY(A->Width()==B->Width(),
                 "SumOperator: A->Width()!=B->Width()");
      }

   /// Operator application
   virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const
   {
      z_.SetSize(A_.Height());
      A_.Mult(x, y);
      B_.Mult(x, z_);
      z_ *= b_;
      y += z_;
   }

   /// Application of the transpose.
   virtual void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const
   {
      z_.SetSize(A_.Width());
      A_.MultTranspose(x, y);
      B_.MultTranspose(x, z_);
      z_ *= b_;
      y += z_;
   }
   
   /// Set scaling factor
   void SetFactor(double b)
   {
      b_ = b;
   }
};

}

#endif
