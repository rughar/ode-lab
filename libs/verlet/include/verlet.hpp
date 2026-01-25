#pragma once
#include <ricatti.hpp>

namespace verlet 
{
    //template<class U>
    //struct verlet_data : public ricatti_core<U> , public ricatti_data<U>
    //{
    //};


    template<class U>
    class verlet_core
    {
        private:
            size_t n;

            template<class V1, class V2>
            void drift(const U h, V1& x, const V2& u) const
            {
                for (size_t i = 0; n; ++i) 
                    x[i] += h * u[i];   
            }

            template<class V1, class V2>
            void kick(const U h, const V1& x, V2& u, verlet_data<U>& data) const
            {
                data.vec = x;
                data.set();   
                data.step(h, u, data);
            }

        public:
            
            template<class V1, class V2>
            void step(const U h, V1& x, V2& u, verlet_data<U>& data) const
            {
                drift(h/2, x, u);
                kick(h, x, u, data);        
                drift(h/2, x, u);
            }
    };
       
}