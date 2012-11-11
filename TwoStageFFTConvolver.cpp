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

#include "TwoStageFFTConvolver.h"

#include "Buffer.h"

#include <algorithm>
#include <cassert>


namespace fftconvolver
{

namespace internal
{

template<typename T>
T NextPowerOf2(const T& val)
{
  T nextPowerOf2 = 1;
  while (nextPowerOf2 < val)
  {
    nextPowerOf2 *= 2;
  }
  return nextPowerOf2;
}

} // End of namespace internal

  

TwoStageFFTConvolver::TwoStageFFTConvolver() :
  _headBlockSize(0),
  _tailBlockSize(0),
  _headConvolver(),
  _tailConvolver(),
  _tailInput(),
  _tailOutput(),
  _tailInputFill(0),
  _tailPos(0)
{
}

  
TwoStageFFTConvolver::~TwoStageFFTConvolver()
{
  TwoStageFFTConvolver::clear();
}

  
void TwoStageFFTConvolver::clear()
{
  _headBlockSize = 0;
  _tailBlockSize = 0;
  _headConvolver.clear();
  _tailConvolver.clear();
  _tailInput.clear();
  _tailOutput.clear();;
  _tailInputFill = 0;
  _tailPos = 0;
}

  
bool TwoStageFFTConvolver::init(size_t headBlockSize,
                                size_t tailBlockSize,
                                const Sample* ir,
                                size_t irLen)
{
  clear();

  if (headBlockSize == 0 || tailBlockSize == 0)
  {
    return false;
  }
  
  if (headBlockSize > tailBlockSize)
  {
    assert(false);
    std::swap(headBlockSize, tailBlockSize);
  }
  
  // Ignore zeros at the end of the impulse response because they only waste computation time
  while (irLen > 0 && ::fabs(ir[irLen-1]) < 0.000001f)
  {
    --irLen;
  }

  if (irLen == 0)
  {
    return true;
  }
  
  _headBlockSize = internal::NextPowerOf2(headBlockSize);
  _tailBlockSize = internal::NextPowerOf2(tailBlockSize);

  const size_t headIrLen = std::min(irLen, _tailBlockSize);
  _headConvolver.init(_headBlockSize, ir, headIrLen, new MultiplyAddEngine());

  if (irLen > _tailBlockSize)
  {
    _tailConvolver.init(_tailBlockSize, ir+_tailBlockSize, irLen-_tailBlockSize, new MultiplyAddEngine());
    _tailInput.resize(_tailBlockSize);
    _tailOutput.resize(_tailBlockSize);
    _tailInputFill = 0;
    _tailPos = 0;
  }
  
  return true;
}


void TwoStageFFTConvolver::process(const Sample* input, Sample* output, size_t len)
{
  // Convolve head
  _headConvolver.process(input, output, len);

  // Convolve tail
  if (_tailOutput.size() > 0 && _tailInput.size() > 0)
  {
    size_t processed = 0;
    while (processed < len)
    {
      const size_t remaining = len - processed;
      const size_t processing = std::min(remaining, _tailBlockSize-_tailInputFill);
      assert(_tailInputFill + processing <= _tailBlockSize);

      // Sum head and tail
      waitForTailConvolution();
      const size_t sumBegin = processed;
      const size_t sumEnd = processed + processing;
      assert(sumEnd <= len);
      for (size_t i=sumBegin; i<sumEnd; ++i)
      {
        assert(_tailPos < _tailBlockSize);
        output[i] += _tailOutput[_tailPos];
        ++_tailPos;
      }
    
      // Append to tail input
      ::memcpy(_tailInput+_tailInputFill, input+processed, processing * sizeof(Sample));
      _tailInputFill += processing;
    
      // Schedule tail convolution if necessary
      assert(_tailInputFill <= _tailBlockSize);
      if (_tailInputFill == _tailBlockSize)
      {
        scheduleTailConvolution();
        _tailInputFill = 0;
        _tailPos = 0;
      }
      processed += processing;
    }
  }
}


void TwoStageFFTConvolver::scheduleTailConvolution()
{
  _tailConvolver.process(_tailInput, _tailOutput, _tailBlockSize);
}


void TwoStageFFTConvolver::waitForTailConvolution()
{
}
    
} // End of namespace fftconvolver
