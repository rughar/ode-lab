// =============================================================================
//  FILE: ricatti.hpp  -  Riccati equation integrator
// =============================================================================
//
//  Purpose
//  -------
//  This file defines the class ricatti_core<U>, which provides a symetric
//  integrator for systems of ordinary differential equations of Riccati type:
//
//      \dot{x}_i = A_i
//                  + sum_j     B_{i,j} x_j
//                  + sum_{j,k} C_{i,j,k} x_j x_k
//
//  The dimension of the system is specified at construction time.
//
//
//  How to use
//  ----------
//  1) Derive a class from ricatti_core<U>
//  2) Pass the system dimension to the base-class constructor
//  3) Override following method to define the vector field:
//
//        void set_coef();
//
//  4) Inside this method, assign coefficients using the macros:
//
//        A_COEF(i)        = value;
//        B_COEF(i, j)     = value;
//        C_COEF(i, j, k)  = value;
//
//  5) Set the current state in the public vector `u`
//  6) Advance the solution by calling step(h)
//
//
//  Notes
//  -----
//  * The state vector `u` is updated in place by each call to step(h).
//
// =============================================================================

#pragma once
#include <vector>
#include <functional>
#include <ling.hpp>

namespace verlet 
{
    template<class U>
    class ricatti_core
    {
        private:

            size_t n;
        
        protected:

            U h_half;
            std::vector<U> mat, vec;    

            struct ACoefProxy 
            {
                ricatti_core& self;
                size_t i;
                void operator=(U value) 
                {
                    self.u[i] += 2 * self.h_half * value;
                }
            };

            struct BCoefProxy 
            {
                ricatti_core& self;
                size_t i, j;
                void operator=(U value) 
                {
                    self.mat[self.n * i + j] -= value * self.h_half;
                    self.u[i] += self.h_half * value * self.vec[j];
                }
            };

            struct CCoefProxy 
            {
                ricatti_core& self;
                size_t i, j, k;
                void operator=(U value) 
                {
                    U tmp = value * self.h_half;
                    self.mat[self.n * i + j] -= tmp * self.vec[k];
                    self.mat[self.n * i + k] -= tmp * self.vec[j];
                }
            };

            ACoefProxy a_coef(const size_t i)
            {
                return ACoefProxy(*this, i);
            };

            BCoefProxy b_coef(const size_t i, const size_t j)
            {
                return BCoefProxy(*this, i, j);
            };

            CCoefProxy c_coef(const size_t i, const size_t j, const size_t k)
            {
                return CCoefProxy(*this, i, j, k);
            };

        public:

            std::vector<U> u;

            virtual void set_coef() {}

            ricatti_core(const size_t size) : n(size)
            {
                vec.resize(n, U(0));
                mat.resize(n * n, U(0));  
            }

            size_t dim() const 
            { 
                return n; 
            }

            void step(const U h)
            {
                h_half = h / 2;
                
                vec = u;
                std::fill(mat.begin(), mat.end(), U(0));
            
                set_coef();

                math::lu_naive(n, mat);    
                math::fb_naive(n, mat, u);     
            }
    };
}