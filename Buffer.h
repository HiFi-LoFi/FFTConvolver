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

#include <algorithm>
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
        assert(!_data && _size == 0);
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
    ::memset(_data, 0, _size * sizeof(T));
  }
  
  void copyFrom(const Buffer<T>& other)
  {
    assert(_size == other._size);
    if (this != &other)
    {
      ::memcpy(_data, other._data, _size * sizeof(T));
    }
  }

  T& operator[](size_t index)
  {
    assert(_data && index < _size);
    return _data[index];
  }

  const T& operator[](size_t index) const
  {
    assert(_data && index < _size);
    return _data[index];
  }

  operator bool() const
  {
    return (_data != 0 && _size > 0);
  }

  T* data()
  {
    return _data;
  }

  const T* data() const
  {
    return _data;
  }

  static void Swap(Buffer<T>& a, Buffer<T>& b)
  {
    std::swap(a._data, b._data);
    std::swap(a._size, b._size);
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
