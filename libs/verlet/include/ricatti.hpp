#pragma once
#include <vector>
#include <ling.hpp>

namespace verlet 
{
    template<class U>
    struct ricatti_data 
    {
        size_t n;
        std::vector<U> vec; 
        std::vector<std::vector<U>> mat;
        void init(size_t dim) 
        {
            n = dim;
            vec.resize(n, U(0));
            mat.resize(n, std::vector<U>(n, U(0)));
        }   
    };
    

    template<class U>
    class ricatti_core
    {
        private:
            size_t n;
        protected:
            std::vector<U> A;
            std::vector<std::vector<U>> B;
            std::vector<std::vector<std::vector<U>>> C;
        public:

            virtual void set_A()
            {
                std::fill(A.begin(), A.end(), U(0));
            }

            virtual void set_B()
            {
                for (auto& row : B)
                    std::fill(row.begin(), row.end(), U(0));
            }
            
            virtual void set_C()
            {
                for (auto& mat : C)
                    for (auto& row : mat)
                        std::fill(row.begin(), row.end(), U(0));
            }           

            void symmetrize_C()
            {
                for (size_t i = 0; i < n; ++i)
                    for (size_t j = 0; j < n; ++j)
                        for (size_t k = 0; k < n; ++k)
                            C[i][j][k] = (C[i][j][k] + C[i][k][j]) / U(2);
            }               
            
            void init(size_t dim) 
            {
                n = dim;
                A.resize(n, U(0));
                B.resize(n, std::vector<U>(n, U(0)));
                C.resize(n, std::vector<std::vector<U>>(n, std::vector<U>(n, U(0))));
            }

            template<class V>
            void step(const U h, V& x, ricatti_data<U>& data) const
            {
                for (size_t i = 0; i < n; ++i) {
                    auto& mat_i = data.mat[i];
                    const auto& B_i = B[i];
                    const auto& C_i = C[i];
                    for (size_t j = 0; j < n; ++j) {
                        U tmp = B_i[j] / 2;
                        const auto& C_ij = C_i[j];
                        for (size_t k = 0; k < n; ++k)
                            tmp += C_ij[k] * x[k];
                        mat_i[j] = -h * tmp;
                    }
                    mat_i[i] += U(1);
                }

                for (size_t i = 0; i < n; ++i) {
                    U tmp = A[i];
                    const auto& B_i = B[i];
                    for (size_t j = 0; j < n; ++j)
                        tmp += B_i[j] * x[j];
                    data.vec[i] = x[i] + h * tmp / 2;
                }

                math::lu_naive(n, data.mat);
                math::fb_naive(n, data.mat, data.vec);                

                for (size_t i = 0; i < n; ++i)
                    x[i] = data.vec[i];
            }
    };
}