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

#include "MultiplyAdd.h"

#include "Buffer.h"

#include <cassert>
#include <cmath>


#if defined (FFTCONVOLVER_USE_SSE)
  #include <xmmintrin.h>
#endif


namespace fftconvolver
{

void MultiplyAdd(Sample* FFTCONVOLVER_RESTRICT re, 
                 Sample* FFTCONVOLVER_RESTRICT im,
                 const Sample* FFTCONVOLVER_RESTRICT reA,
                 const Sample* FFTCONVOLVER_RESTRICT imA,
                 const Sample* FFTCONVOLVER_RESTRICT reB,
                 const Sample* FFTCONVOLVER_RESTRICT imB,
                 const size_t len)
{
#if defined(FFTCONVOLVER_USE_SSE)
  const size_t end4 = 4 * (len / 4);
  for (size_t i=0; i<end4; i+=4)
  {
    __m128 mul1;
    __m128 mul2;
    __m128 ra = _mm_load_ps(&reA[i]);
    __m128 rb = _mm_load_ps(&reB[i]);
    __m128 ia = _mm_load_ps(&imA[i]);
    __m128 ib = _mm_load_ps(&imB[i]);
    __m128 real = _mm_load_ps(&re[i]);
    __m128 imag = _mm_load_ps(&im[i]);
    
    mul1 = _mm_mul_ps(ra, rb);
    mul2 = _mm_mul_ps(ia, ib);
    real = _mm_add_ps(real, mul1);
    real = _mm_sub_ps(real, mul2);
    _mm_store_ps(&re[i], real);
    
    mul1 = _mm_mul_ps(ra, ib);
    mul2 = _mm_mul_ps(ia, rb);
    imag = _mm_add_ps(imag, mul1);
    imag = _mm_add_ps(imag, mul2);
    _mm_store_ps(&im[i], imag);
  }
  for (size_t i=end4; i<len; ++i)
  {
    re[i] += reA[i] * reB[i] - imA[i] * imB[i];
    im[i] += reA[i] * imB[i] + imA[i] * reB[i];
  }
#else
  const size_t end4 = 4 * (len / 4);
  for (size_t i=0; i<end4; i+=4)
  {
    re[i+0] += reA[i+0] * reB[i+0] - imA[i+0] * imB[i+0];
    re[i+1] += reA[i+1] * reB[i+1] - imA[i+1] * imB[i+1];
    re[i+2] += reA[i+2] * reB[i+2] - imA[i+2] * imB[i+2];
    re[i+3] += reA[i+3] * reB[i+3] - imA[i+3] * imB[i+3];
    im[i+0] += reA[i+0] * imB[i+0] + imA[i+0] * reB[i+0];
    im[i+1] += reA[i+1] * imB[i+1] + imA[i+1] * reB[i+1];
    im[i+2] += reA[i+2] * imB[i+2] + imA[i+2] * reB[i+2];
    im[i+3] += reA[i+3] * imB[i+3] + imA[i+3] * reB[i+3];
  }
  for (size_t i=end4; i<len; ++i)
  {
    re[i] += reA[i] * reB[i] - imA[i] * imB[i];
    im[i] += reA[i] * imB[i] + imA[i] * reB[i];
  }
#endif
}


void MultiplyAdd(SplitComplex& result, const SplitComplex& a, const SplitComplex& b)
{
  assert(result.size() == a.size());
  assert(result.size() == b.size());
  MultiplyAdd(result.re(), result.im(), a.re(), a.im(), b.re(), b.im(), result.size());
}

} // End of namespace fftconvolver
