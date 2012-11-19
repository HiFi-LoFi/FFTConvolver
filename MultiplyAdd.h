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

void MultiplyAdd(SplitComplex& result, const SplitComplex& a, const SplitComplex& b);

} // End of namespace fftconvolver

#endif // Header guard
