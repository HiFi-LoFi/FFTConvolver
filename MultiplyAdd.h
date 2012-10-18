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

#ifndef _FFTCONVOLVER_PREMULTIPLYADD_H
#define _FFTCONVOLVER_PREMULTIPLYADD_H

#include "Configuration.h"
#include "SplitComplex.h"

#include <cstddef>
#include <vector>


namespace fftconvolver
{ 
  
void MultiplyAdd(Sample* FFTCONVOLVER_RESTRICT re, 
                 Sample* FFTCONVOLVER_RESTRICT im,
                 const Sample* FFTCONVOLVER_RESTRICT reA,
                 const Sample* FFTCONVOLVER_RESTRICT imA,
                 const Sample* FFTCONVOLVER_RESTRICT reB,
                 const Sample* FFTCONVOLVER_RESTRICT imB,
                 const size_t len);


// ======================================================
  
  
class MultiplyAddEngine
{
public:
  MultiplyAddEngine();
  virtual ~MultiplyAddEngine();
    
  virtual void init(const std::vector<SplitComplex*>& ir);
  
  virtual void setAudio(size_t index, const SplitComplex& audio);
  
  struct Pair
  {
    size_t indexIr;
    size_t indexAudio;
  };
  
  virtual void multiplyAdd(const std::vector<Pair>& pairs);

  virtual const SplitComplex& getResult();
  
private:
  void clear();
  
  std::vector<SplitComplex*> _ir;
  std::vector<SplitComplex*> _audio;
  SplitComplex _result;
  
  // Prevent uncontrolled usage
  MultiplyAddEngine(const MultiplyAddEngine&);
  MultiplyAddEngine& operator=(const MultiplyAddEngine&);
};

} // End of namespace fftconvolver

#endif // Header guard
