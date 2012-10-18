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

#ifndef _FFTCONVOLVER_FFTCONVOLVER_H
#define _FFTCONVOLVER_FFTCONVOLVER_H

#include "AudioFFT.h"
#include "Buffer.h"
#include "MultiplyAdd.h"
#include "Sample.h"
#include "SplitComplex.h"

#include <cstddef>

#include <vector>


namespace fftconvolver
{ 

class FFTConvolver
{  
public:
  FFTConvolver();  
  virtual ~FFTConvolver();
  
  virtual void clear();
  virtual bool init(size_t blockSize,
                    const Sample* ir,
                    size_t irLen,
                    MultiplyAddEngine* muliplyAddEngine);
  virtual void process(const Sample* input, Sample* output, size_t len);
  
protected:
  virtual void multiplyAdd(SplitComplex& result, const SplitComplex& a, const SplitComplex& b) const;
  
private:  
  size_t _blockSize;
  size_t _segSize;
  size_t _segCount;
  size_t _fftComplexSize;
  std::vector<SplitComplex*> _segments;
  std::vector<SplitComplex*> _segmentsIR;
  SampleBuffer _fftBuffer;
  audiofft::AudioFFT _fft;
  SplitComplex _preMultiplied;
  SplitComplex _conv;
  SampleBuffer _overlap;
  size_t _current;
  SampleBuffer _inputBuffer;
  size_t _inputBufferFill;
  MultiplyAddEngine* _multiplyAddEngine;
  std::vector<MultiplyAddEngine::Pair> _preMultiplyAddPairs;
};
  
} // End of namespace fftconvolver

#endif // Header guard
