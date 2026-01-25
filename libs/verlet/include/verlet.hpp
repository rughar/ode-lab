#pragma once
#include <ricatti.hpp>

namespace verlet 
{
    template<class U>
    class verlet_core : public ricatti_core<U>, protected ricatti_coef<U>
    {
        private:

            void drift(const U h)
            {
                for (size_t i = 0; dim(); ++i) 
                    ricatti_core<U>::x[i] += h * u[i];   
            }

            void kick(const U h)
            {
                ricatti_core<U>::step(h);
            }

        public:

            std::vector<U> u;

            size_t dim() const 
            { 
                return ricatti_coef<U>::dim(); 
            }   

            void init(const size_t n) 
            {
                ricatti_coef<U>::init(n);
                ricatti_core<U>::init(this);
                u.resize(n, U(0));
            }

            void step(const U h)
            {
                drift(h/2);
                kick(h);        
                drift(h/2);
            }
    };
       
}