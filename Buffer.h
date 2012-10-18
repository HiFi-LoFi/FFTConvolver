// ==================================================================================
// Copyright (c) 2012 HiFi-LoFi
//
// This is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ==================================================================================

#ifndef _FFTCONVOLVER_BUFFER_H
#define _FFTCONVOLVER_BUFFER_H

#include "Configuration.h"

#include <cassert>
#include <cstddef>
#include <cstring>


namespace fftconvolver
{

template<typename T>
class Buffer
{
public:  
  explicit Buffer(size_t initialSize = 0) :
    _data(0),
    _size(0)
  {
    resize(initialSize);
  }
  
  virtual ~Buffer()
  {
    clear();
  }
    
  void clear()
  {
    internal::DeallocateBuffer(_data);
    _data = 0;
    _size = 0;
  }
  
  void resize(size_t size)
  {
    if (_size != size)
    {
      clear();
      
      if (size > 0)
      {
        _data = internal::AllocateBuffer<T>(size);
        _size = size;
      }
    }
    setZero();
  }
  
  size_t size() const
  {
    return _size;
  }
  
  void setZero()
  {
    memset(_data, 0, _size * sizeof(T));
  }
  
  void copyFrom(const Buffer<T>& other)
  {
    assert(_size == other._size);
    memcpy(_data, other._data, _size * sizeof(T));
  }
  
  operator T*()
  {
    return _data;
  }
  
  operator const T*() const
  {
    return _data;
  }
  
  T* data()
  {
    return _data;
  }
  
  const T* data() const
  {
    return _data;
  }
  
private:
  T* _data;
  size_t _size;
  
  // Prevent uncontrolled usage
  Buffer(const Buffer&);
  Buffer& operator=(const Buffer&);
};
  
} // End of namespace

#endif // Header guard
