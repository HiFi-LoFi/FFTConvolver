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

#include "FFTConvolver.h"

#include "Buffer.h"
#include "Configuration.h"

#include <cassert>
#include <cmath>


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
  
  
template<typename T>
void Sum(T* FFTCONVOLVER_RESTRICT result, const T* FFTCONVOLVER_RESTRICT a, const T* FFTCONVOLVER_RESTRICT b, size_t len)
{
  const size_t end8 = 8 * (len / 8);
  for (size_t i=0; i<end8; i+=8)
  {
    result[i+0] = a[i+0] + b[i+0];
    result[i+1] = a[i+1] + b[i+1];
    result[i+2] = a[i+2] + b[i+2];
    result[i+3] = a[i+3] + b[i+3];
    result[i+4] = a[i+4] + b[i+4];
    result[i+5] = a[i+5] + b[i+5];
    result[i+6] = a[i+6] + b[i+6];
    result[i+7] = a[i+7] + b[i+7];
  }
  for (size_t i=end8; i<len; ++i)
  {
    result[i] = a[i] + b[i];
  }
}


template<typename T>
void CopyAndPad(Buffer<T>& dest, const T* src, size_t srcSize)
{
  assert(dest.size() >= srcSize);
  ::memcpy(dest, src, srcSize * sizeof(T));
  ::memset(dest + srcSize, 0, (dest.size()-srcSize) * sizeof(T)); 
}
  
} // End of namespace internal

  

FFTConvolver::FFTConvolver() :
  _blockSize(0),
  _segSize(0),
  _segCount(0),
  _fftComplexSize(0),
  _segments(),
  _segmentsIR(),
  _fftBuffer(),
  _fft(),
  _preMultiplied(),
  _conv(),
  _overlap(),
  _current(0),
  _inputBuffer(),
  _inputBufferFill(0),
  _multiplyAddEngine(0),
  _preMultiplyAddPairs()
{
}

  
FFTConvolver::~FFTConvolver()
{
  FFTConvolver::clear();
}

  
void FFTConvolver::clear()
{
  delete _multiplyAddEngine;
  _multiplyAddEngine = 0;
  
  for (size_t i=0; i<_segCount; ++i)
  {
    delete _segments[i];
    delete _segmentsIR[i];
  }
  
  _blockSize = 0;
  _segSize = 0;
  _segCount = 0;
  _fftComplexSize = 0;
  _segments.clear();
  _segmentsIR.clear();
  _fftBuffer.clear();
  _fft.init(0);
  _preMultiplied.clear();
  _conv.clear();
  _overlap.clear();
  _current = 0;
  _inputBuffer.clear();
  _inputBufferFill = 0;
  _preMultiplyAddPairs.clear();
}

  
bool FFTConvolver::init(size_t blockSize, const Sample* ir, size_t irLen, MultiplyAddEngine* muliplyAddEngine)
{
  clear();

  if (blockSize == 0 || !muliplyAddEngine)
  {
    return false;
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
  
  _blockSize = internal::NextPowerOf2(blockSize);
  _segSize = 2 * _blockSize;
  _segCount = static_cast<size_t>(::ceil(static_cast<float>(irLen) / static_cast<float>(_blockSize)));
  _fftComplexSize = audiofft::AudioFFT::ComplexSize(_segSize);
  
  // FFT
  _fft.init(_segSize);
  _fftBuffer.resize(_segSize);
  
  // Prepare segments
  for (size_t i=0; i<_segCount; ++i)
  {
    _segments.push_back(new SplitComplex(_fftComplexSize));    
  }
  
  // Prepare IR
  for (size_t i=0; i<_segCount; ++i)
  {
    SplitComplex* segment = new SplitComplex(_fftComplexSize);
    const size_t remaining = irLen - (i * _blockSize);
    const size_t sizeCopy = (remaining >= _blockSize) ? _blockSize : remaining;
    internal::CopyAndPad(_fftBuffer, &ir[i*_blockSize], sizeCopy);
    _fft.fft(_fftBuffer, segment->re(), segment->im());
    _segmentsIR.push_back(segment);
  }
  
  // Prepare convolution buffers  
  _preMultiplied.resize(_fftComplexSize);
  _conv.resize(_fftComplexSize);
  _overlap.resize(_blockSize);
  
  // Prepare input buffer
  _inputBuffer.resize(_blockSize);
  _inputBufferFill = 0;

  // Reset current position
  _current = 0;
  
  // Prepare pre-multiply-add
  _multiplyAddEngine = muliplyAddEngine;
  _multiplyAddEngine->init(_segmentsIR);
  if (_segCount >= 2)
  {
    _preMultiplyAddPairs.resize(_segCount-2);
  }
  
  return true;
}


void FFTConvolver::process(const Sample* input, Sample* output, size_t len)
{
  if (_segCount == 0)
  {
    ::memset(output, 0, len * sizeof(Sample));
    return;
  }

  size_t processed = 0;
  while (processed < len)
  {
    const bool inputBufferWasEmpty = (_inputBufferFill == 0);
    const size_t processing = std::min(len-processed, _blockSize-_inputBufferFill);
    const size_t inputBufferPos = _inputBufferFill;
    ::memcpy(_inputBuffer+inputBufferPos, input+processed, processing * sizeof(Sample));
    
    // Forward FFT
    internal::CopyAndPad(_fftBuffer, &_inputBuffer[0], _blockSize); 
    _fft.fft(_fftBuffer, _segments[_current]->re(), _segments[_current]->im());
    
    // Complex multiplication
    if (inputBufferWasEmpty)
    {
      if (_segCount < 2)
      {
        _preMultiplied.setZero();
      }
      if (_segCount >= 2)
      {
        _preMultiplied.copyFrom(_multiplyAddEngine->getResult());
        multiplyAdd(_preMultiplied, *_segments[(_current + 1) % _segCount], *_segmentsIR[1]);
      }
    }
    _conv.copyFrom(_preMultiplied);
    multiplyAdd(_conv, *_segments[_current], *_segmentsIR[0]);
    
    // Backward FFT
    _fft.ifft(_fftBuffer, _conv.re(), _conv.im());
    
    // Add overlap
    internal::Sum(output+processed, _fftBuffer+inputBufferPos, _overlap+inputBufferPos, processing);
    
    // Input buffer full => Next block
    _inputBufferFill += processing;
    if (_inputBufferFill == _blockSize)
    {
      // Input buffer is empty again now
      _inputBuffer.setZero();
      _inputBufferFill = 0;
      
      // Save the overlap
      ::memcpy(_overlap, _fftBuffer + _blockSize, _blockSize * sizeof(Sample));
      
      // Update current segment
      _multiplyAddEngine->setAudio(_current, *_segments[_current]);
      _current = (_current > 0) ? (_current - 1) : (_segCount - 1);
      
      // Trigger pre-calculation
      if (_segCount > 2)
      {
        for (size_t i=2; i<_segCount; ++i)
        {
          _preMultiplyAddPairs[i-2].indexIr = i;
          _preMultiplyAddPairs[i-2].indexAudio = (_current + i) % _segCount;
        }
        _multiplyAddEngine->multiplyAdd(_preMultiplyAddPairs);
      }
    }
    
    processed += processing;
  }
}
  
  
void FFTConvolver::multiplyAdd(SplitComplex& result, const SplitComplex& a, const SplitComplex& b) const
{
  assert(result.size() == a.size());
  assert(result.size() == b.size());
  MultiplyAdd(result.re(), result.im(), a.re(), a.im(), b.re(), b.im(), result.size());
}
  
} // End of namespace fftconvolver
