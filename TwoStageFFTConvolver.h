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

#ifndef _FFTCONVOLVER_TWOSTAGEFFTCONVOLVER_H
#define _FFTCONVOLVER_TWOSTAGEFFTCONVOLVER_H

#include "FFTConvolver.h"
#include "Sample.h"


namespace fftconvolver
{ 

class TwoStageFFTConvolver
{  
public:
  TwoStageFFTConvolver();  
  virtual ~TwoStageFFTConvolver();
  
  bool init(size_t headBlockSize, size_t tailBlockSize, const Sample* ir, size_t irLen);  
  void process(const Sample* input, Sample* output, size_t len);
  void reset();
  
protected:
  virtual void startBackgroundProcessing();
  virtual void waitForBackgroundProcessing();

  void doBackgroundProcessing();

private:
  size_t _headBlockSize;
  size_t _tailBlockSize;
  FFTConvolver _headConvolver;
  FFTConvolver _tailConvolver0;
  SampleBuffer _tailOutput0;
  SampleBuffer _tailPrecalculated0;
  FFTConvolver _tailConvolver;
  SampleBuffer _tailOutput;
  SampleBuffer _tailPrecalculated;
  SampleBuffer _tailInput;
  size_t _tailInputFill;
  size_t _precalculatedPos;
  SampleBuffer _backgroundProcessingInput;

  // Prevent uncontrolled usage
  TwoStageFFTConvolver(const TwoStageFFTConvolver&);
  TwoStageFFTConvolver& operator=(const TwoStageFFTConvolver&);
};
  
} // End of namespace fftconvolver

#endif // Header guard
