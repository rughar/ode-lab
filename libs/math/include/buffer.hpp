#pragma once
#include <vector>
#include <deque>

namespace math
{
  template <class U>
  class vec_buffer
  {
  public:
    ~vec_buffer();
    void set(size_t size);
    U *take();
    void free(U *link);

  private:
    size_t n;
    std::vector<U *> data;
    std::vector<bool> taken;
  };

  template <class U>
  inline vec_buffer<U>::~vec_buffer()
  {
    for (U* d : data) 
      delete[] d;
  }

  template <class U>
  inline void vec_buffer<U>::set(size_t size)
  {
    n = size;
  }

  template <class U>
  inline U *vec_buffer<U>::take()
  {
    for (size_t i = 0; i < data.size(); ++i)
      if (!taken[i])
      {
        taken[i] = 1;
        return data[i];
      }

    U* vec = new U[n];
    data.push_back(vec);
    taken.push_back(1);
    return vec;
  }

  template <class U>
  inline void vec_buffer<U>::free(U *link)
  {
    for (size_t i = 0; i < data.size(); ++i)
      if (data[i] == link) {
        taken[i] = 0;
        return;
      }
  }

}