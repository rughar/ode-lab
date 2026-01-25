#pragma once
#include <ricatti.hpp>

namespace verlet 
{
    template<class U>
    class verlet_core : public ricatti_core<U>
    {
        private:

            void drift(const U h)
            {
                for (size_t i = 0; i < ricatti_core<U>::dim(); ++i) 
                    x[i] += h * ricatti_core<U>::u[i];   
            }

            void kick(const U h)
            {
                update_coefs();
                ricatti_core<U>::step(h);
            }

        public:

            std::vector<U> x;
 
            void init(const size_t n) 
            {
                x.resize(n, U(0));
            }

            virtual void update_coefs() = 0;

            void step(const U h)
            {
                drift(h/2);
                kick(h);        
                drift(h/2);
            }
    };
       
}