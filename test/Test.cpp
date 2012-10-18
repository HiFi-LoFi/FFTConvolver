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

#include <algorithm>
#include <cstring>
#include <vector>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../FFTConvolver.h"
#include "../MultiplyAdd.h"


template<typename T>
void SimpleConvolve(const T* input, size_t inLen, const T* ir, size_t irLen, T* output)
{
  if (irLen > inLen)
  {
    SimpleConvolve(ir, irLen, input, inLen, output);
    return;
  }
  
  ::memset(output, 0, (inLen+irLen-1) * sizeof(T));
  
  for (size_t n=0; n<irLen; ++n)
  {
    for (size_t m=0; m<=n; ++m)
    {
      output[n] += ir[m] * input[n-m];
    }
  }
  
  for (size_t n=irLen; n<inLen; ++n)
  {
    for (size_t m=0; m<irLen; ++m)
    {
      output[n] += ir[m] * input[n-m];
    }
  }
  
  for (size_t n=inLen; n<inLen+irLen-1; ++n)
  {
    for (size_t m=n-inLen+1; m<irLen; ++m)
    {
      output[n] += ir[m] * input[n-m];
    }
  }
}


static bool TestConvolver(size_t inputSize,
                          size_t irSize,
                          size_t blockSizeMin,
                          size_t blockSizeMax,
                          size_t blockSizeConvolver,
                          bool refCheck,
                          fftconvolver::MultiplyAddEngine* multiplyAddEngine)
{
  // Prepare input and IR
  std::vector<fftconvolver::Sample> in(inputSize);
  for (size_t i=0; i<inputSize; ++i)
  {
    in[i] = static_cast<float>(::cos(static_cast<float>(i) / 123.0));
  }
  
  std::vector<fftconvolver::Sample> ir(irSize);
  for (size_t i=0; i<irSize; ++i)
  {
    ir[i] = static_cast<float>(::sin(static_cast<float>(i) / 79.0));
  }
  
  // Simple convolver
  std::vector<fftconvolver::Sample> outSimple(in.size() + ir.size() - 1, fftconvolver::Sample(0.0));
  if (refCheck)
  {
    SimpleConvolve(&in[0], in.size(), &ir[0], ir.size(), &outSimple[0]);
  }
  
  // Orgami convolver
  std::vector<fftconvolver::Sample> out(in.size() + ir.size() - 1, fftconvolver::Sample(0.0));
  {
    fftconvolver::FFTConvolver convolver;
    convolver.init(blockSizeConvolver, &ir[0], ir.size(), multiplyAddEngine);
    std::vector<fftconvolver::Sample> inBuf(blockSizeMax);
    size_t processedOut = 0;
    size_t processedIn = 0;
    while (processedOut < out.size())
    {
      const size_t blockSize = blockSizeMin + (static_cast<size_t>(rand()) % (1+(blockSizeMax-blockSizeMin))); 
      
      const size_t remainingOut = out.size() - processedOut;
      const size_t remainingIn = in.size() - processedIn;
      
      const size_t processingOut = std::min(remainingOut, blockSize);
      const size_t processingIn = std::min(remainingIn, blockSize);
      
      memset(&inBuf[0], 0, inBuf.size() * sizeof(fftconvolver::Sample));
      if (processingIn > 0)
      {
        memcpy(&inBuf[0], &in[processedIn], processingIn * sizeof(fftconvolver::Sample));
      }
      
      convolver.process(&inBuf[0], &out[processedOut], processingOut);
      
      processedOut += processingOut;
      processedIn += processingIn;
    }
  }
  
  if (refCheck)
  {
    size_t diffSamples = 0;
    for (size_t i=0; i<outSimple.size(); ++i)
    {
      const fftconvolver::Sample a = out[i];
      const fftconvolver::Sample b = outSimple[i];
      
      if (::fabs(a-b) > 0.05)
      {
        ++diffSamples;
      }
    }
    printf("Correctness Test (input %zu, IR %zu, blocksize %zu-%zu) => %s\n", inputSize, irSize, blockSizeMin, blockSizeMax, (diffSamples == 0) ? "[OK]" : "[FAILED]");
    return (diffSamples == 0);
  }
  else
  {
    printf("Performance Test (input %zu, IR %zu, blocksize %zu-%zu) => Completed\n", inputSize, irSize, blockSizeMin, blockSizeMax);
    return true;
  }
}



int main()
{ 
  // Correctness
#define TEST_CORRECTNESS
#ifdef TEST_CORRECTNESS
  TestConvolver(1, 1, 1, 1, 1, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(2, 2, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(3, 3, 3, 3, 3, true, new fftconvolver::MultiplyAddEngine());
  
  TestConvolver(3, 2, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(4, 2, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(4, 3, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(9, 4, 3, 3, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(171, 7, 5, 5, 5, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(1979, 17, 7, 7, 5, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100, 10, 3, 5, 5, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(123, 45, 12, 34, 34, true, new fftconvolver::MultiplyAddEngine());
  
  TestConvolver(2, 3, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(2, 4, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(3, 4, 2, 2, 2, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(4, 9, 3, 3, 3, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(7, 171, 5, 5, 5, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(17, 1979, 7, 7, 7, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(10, 100, 3, 5, 5, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(45, 123, 12, 34, 34, true, new fftconvolver::MultiplyAddEngine());
  
  TestConvolver(100000, 1234, 100,  128,  128, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 1234, 100,  256,  256, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 1234, 100,  512,  512, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 1234, 100, 1024, 1024, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 1234, 100, 2048, 2048, true, new fftconvolver::MultiplyAddEngine());
  
  TestConvolver(100000, 4321, 100,  128,  128, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 4321, 100,  256,  256, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 4321, 100,  512,  512, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 4321, 100, 1024, 1024, true, new fftconvolver::MultiplyAddEngine());
  TestConvolver(100000, 4321, 100, 2048, 2048, true, new fftconvolver::MultiplyAddEngine());
#endif
  
//#define TEST_PERFORMANCE
#ifdef TEST_PERFORMANCE
  // Performance
  TestConvolver(60*44100, 5*44100, 50, 100, 1024, false, new fftconvolver::MultiplyAddEngine());
#endif
  
  return 0;
}
