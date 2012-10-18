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

#ifndef Convolution_SplitComplex_h
#define Convolution_SplitComplex_h

#include "Sample.h"


namespace fftconvolver
{

class SplitComplex
{
public:
  explicit SplitComplex(size_t initialSize = 0) : 
    _size(0),
    _re(),
    _im()
  {
    resize(initialSize);
  }
  
  ~SplitComplex()
  {
    clear();
  }
  
  void clear()
  {
    _re.clear();
    _im.clear();
    _size = 0;
  }
  
  void resize(size_t newSize)
  {
    _re.resize(newSize);
    _im.resize(newSize);
    _size = newSize;
  }
  
  void setZero()
  {
    _re.setZero();
    _im.setZero();
  }
  
  void copyFrom(const SplitComplex& other)
  {
    _re.copyFrom(other._re);
    _im.copyFrom(other._im);
  }
  
  Sample* re()
  {
    return _re;
  }
  
  const Sample* re() const
  {
    return _re;
  }
  
  Sample* im()
  {
    return _im;
  }
  
  const Sample* im() const
  {
    return _im;
  }
  
  size_t size() const
  {
    return _size;
  }
  
private:
  size_t _size;
  SampleBuffer _re;
  SampleBuffer _im;
  
  // Prevent uncontrolled usage
  SplitComplex(const SplitComplex&);
  SplitComplex& operator=(const SplitComplex&);
};
  
} // End of namespace

#endif // Header guard
