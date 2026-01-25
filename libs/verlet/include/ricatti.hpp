// =============================================================================
//  FILE: ricatti.hpp  -  Riccati equation integrator
// =============================================================================
//  Purpose :
//    * ricatti_coef  - holds coefficients of a (quadratic) vector field:
//          f(x)_i = A_i + sum_j B_{i,j} x_j + sum_{j,k} C_{i,j,k} x_j x_k
//      plus helper to symmetrise C in the last two indices.
//    * ricatti_core  - performs one kvasi-implicit step for state x using
//      a linear system solve based on a locally assembled matrix.
//
//  Notes :
//    - Uses math::lu_naive + math::fb_naive (no pivoting). Numeric stability
//      relies on diagonal dominance / well-conditioned system matrix.
//    - The step builds an n*n matrix "mat" and a vector "vec", then solves
//      mat . vec = rhs in place, and swaps vec into x.
// =============================================================================

#pragma once
#include <vector>
#include <ling.hpp>

namespace verlet 
{

// -------------------------------------------------------------------------
//  ricatti_coef
// -------------------------------------------------------------------------
//  Coefficients describing a quadratic (Riccati-like) vector field.
//
//  Stored data:
//    n   - dimension
//    A   - constant term (size n)
//    B   - linear term matrix (size n*n)
//    C   - quadratic term tensor (size n*n*n), interpreted as:
//            C_{i,j,k} multiplies x_k in the contribution for component i,j
// -------------------------------------------------------------------------

    template<class U>
    class ricatti_coef 
    {
        private:
            
            size_t n;

        public:

            std::vector<U> A; 
            std::vector<std::vector<U>> B;
            std::vector<std::vector<std::vector<U>>> C;

            size_t dim() const 
            { 
                return n; 
            }   

            void init(size_t dim) 
            {
                n = dim;
                A.resize(dim, U(0));
                B.resize(dim, std::vector<U>(dim, U(0)));
                C.resize(dim, std::vector<std::vector<U>>(dim, std::vector<U>(dim, U(0))));
            }   

            void symmetrize_C()
            {
                for (size_t i = 0; i < n; ++i)
                    for (size_t j = 0; j < n; ++j)
                        for (size_t k = 0; k < n; ++k)
                            C[i][j][k] = (C[i][j][k] + C[i][k][j]) / U(2);
            }  
    }; 
    
// -------------------------------------------------------------------------
//  ricatti_core
// -------------------------------------------------------------------------
//  Core stepping routine operating on state vector x (size n).
//
//  Internal workspaces:
//    vec - system vector workspace (size n)
//    mat - system matrix workspace (size n*n), overwritten by LU factors
//
//  External dependency:
//    coef - pointer to coefficient set (must be valid during step()).
//
//  Limitations:
//    - math::lu_naive has no pivoting -> requires numerically safe matrix.
//    - mat and vec are overwritten in-place during step().
// -------------------------------------------------------------------------

    template<class U>
    class ricatti_core
    {
        private:

            std::vector<U> vec; 
            std::vector<std::vector<U>> mat;

        protected:
        
            ricatti_coef<U>* coef;
    
        public:

            std::vector<U> x;

            size_t dim() const 
            { 
                return coef->dim(); 
            }

            void init(ricatti_coef<U>* coefficients) 
            {
                coef = coefficients;
                x.resize(dim(), U(0));
                vec.resize(dim(), U(0));
                mat.resize(dim(), std::vector<U>(dim(), U(0)));   
            }
 
            void step(const U h)
            {
                for (size_t i = 0; i < dim(); ++i) {
                    auto& mat_i = mat[i];
                    const auto& B_i = coef->B[i];
                    const auto& C_i = coef->C[i];
                    for (size_t j = 0; j < dim(); ++j) {
                        U tmp = B_i[j] / 2;
                        const auto& C_ij = C_i[j];
                        for (size_t k = 0; k < dim(); ++k)
                            tmp += C_ij[k] * x[k];
                        mat_i[j] = -h * tmp;
                    }
                    mat_i[i] += U(1);
                }

                for (size_t i = 0; i < dim(); ++i) {
                    U tmp = coef->A[i];
                    const auto& B_i = coef->B[i];
                    for (size_t j = 0; j < dim(); ++j)
                        tmp += B_i[j] * x[j];
                    vec[i] = x[i] + h * tmp / 2;
                }

                math::lu_naive(dim(), mat);
                math::fb_naive(dim(), mat, vec);                

                x.swap(vec);
            }
    };
}