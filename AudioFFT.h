// ==================================================================================
// Copyright (c) 2013 HiFi-LoFi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// ==================================================================================

#ifndef _AUDIOFFT_H
#define _AUDIOFFT_H


/**
* AudioFFT provides real-to-complex/complex-to-real FFT routines.
* 
* Features:
*
* - Real-complex FFT and complex-real inverse FFT for power-of-2-sized real data.
*
* - Uniform interface to different FFT implementations (currently Ooura, FFTW3 and Apple Accelerate).
*
* - Complex data is handled in "split-complex" format, i.e. there are separate
*   arrays for the real and imaginary parts which can be useful for SIMD optimizations
*   (split-complex arrays have to be of length (size/2+1) representing bins from DC
*   to Nyquist frequency).
*
* - Output is "ready to use" (all scaling etc. is already handled internally).
*
* - No allocations/deallocations after the initialization which makes it usable
*   for real-time audio applications (that's what I wrote it for and using it).
*
*
* How to use it in your project:
* 
* - Add the .h and .cpp file to your project - that's all.
*
* - To get extra speed, you can link FFTW3 to your project and define
*   AUDIOFFT_FFTW3 (however, please check whether your project suits the
*   according license).
*
* - To get the best speed on Apple platforms, you can link the Apple
*   Accelerate framework to your project and define
*   AUDIOFFT_APPLE_ACCELERATE  (however, please check whether your
*   project suits the according license).
*
*
* Remarks:
*
* - AudioFFT is not intended to be the fastest FFT, but to be a fast-enough
*   FFT suitable for most audio applications.
*
* - AudioFFT uses the quite liberal MIT license.
*
*
* Example usage:
* @code
* #include "AudioFFT.h"
*
* void Example()
* {
*   const size_t fftSize = 1024; // Needs to be power of 2!
*
*   std::vector<float> input(fftSize, 0.0f);
*   std::vector<float> re(fftaudio::AudioFFT::ComplexSize(fftSize)); 
*   std::vector<float> im(fftaudio::AudioFFT::ComplexSize(fftSize)); 
*   std::vector<float> output(fftSize);
*
*   audiofft::AudioFFT fft;
*   fft.init(1024);
*   fft.fft(input.data(), re.data(), im.data());
*   fft.ifft(output.data(), re.data(), im.data());
* }
* @endcode
*/


#include <cstddef>


#if defined(AUDIOFFT_APPLE_ACCELERATE)
  #define AUDIOFFT_APPLE_ACCELERATE_USED
  #include <Accelerate/Accelerate.h>
  #include <vector>
#elif defined (AUDIOFFT_FFTW3)
  #define AUDIOFFT_FFTW3_USED
  #include <fftw3.h>
#else
  #ifndef AUDIOFFT_OOURA
    #define AUDIOFFT_OOURA
  #endif
  #define AUDIOFFT_OOURA_USED
  #include <vector>
#endif


namespace audiofft
{

/**
* @class AudioFFTBase
* @brief Interface class for FFT implementations
*/
class AudioFFTBase
{
public:
  /**
  * @brief Constructor
  */
  AudioFFTBase()
  {
  }

  /**
  * @brief Destructor
  */
  virtual ~AudioFFTBase()
  {
  }

  /**
  * @brief Initializes the FFT object
  * @param size Size of the real input (must be power 2)
  */
  virtual void init(size_t size) = 0;
  
  /**
  * @brief Performs the forward FFT
  * @param data The real input data (has to be of the length as specified in init())
  * @param re The real part of the complex output (has to be of length as returned by ComplexSize())
  * @param im The imaginary part of the complex output (has to be of length as returned by ComplexSize())
  */
  virtual void fft(const float* data, float* re, float* im) = 0;

  /**
  * @brief Performs the inverse FFT
  * @param data The real output data (has to be of the length as specified in init())
  * @param re The real part of the complex input (has to be of length as returned by ComplexSize())
  * @param im The imaginary part of the complex input (has to be of length as returned by ComplexSize())
  */
  virtual void ifft(float* data, const float* re, const float* im) = 0;

  /**
  * @brief Calculates the necessary size of the real/imaginary complex arrays
  * @param size The size of the real data
  * @return The size of the real/imaginary complex arrays
  */
  static size_t ComplexSize(size_t size)
  {
    return (size / 2) + 1;
  }
  
private:
  AudioFFTBase(const AudioFFTBase&);
  AudioFFTBase& operator=(const AudioFFTBase&);
};


// ================================================================


#ifdef AUDIOFFT_APPLE_ACCELERATE_USED


/**
* @class AppleAccelerateFFT
* @brief FFT implementation using the Apple Accelerate framework internally
*/
class AppleAccelerateFFT : public AudioFFTBase
{
public:
  AppleAccelerateFFT();
  virtual ~AppleAccelerateFFT();
  virtual void init(size_t size);
  virtual void fft(const float* data, float* re, float* im);
  virtual void ifft(float* data, const float* re, const float* im);

private:
  size_t _size;
  size_t _powerOf2;
  FFTSetup _fftSetup;
  std::vector<float> _re;
  std::vector<float> _im;
  
  // Prevent uncontrolled usage
  AppleAccelerateFFT(const AppleAccelerateFFT&);
  AppleAccelerateFFT& operator=(const AppleAccelerateFFT&);
};

typedef AppleAccelerateFFT AudioFFT;

#endif // AUDIOFFT_APPLE_ACCELERATE_USED


// ================================================================


#ifdef AUDIOFFT_FFTW3_USED

/**
* @class FFTW3FFT
* @brief FFT implementation using FFTW3 internally (see fftw.org)
*/
class FFTW3FFT : public AudioFFTBase
{
public:
  FFTW3FFT();
  virtual ~FFTW3FFT();
  virtual void init(size_t size);
  virtual void fft(const float* data, float* re, float* im);
  virtual void ifft(float* data, const float* re, const float* im);

private:
  size_t _size;
  size_t _complexSize;
  fftwf_plan _planForward;
  fftwf_plan _planBackward;
  float* _data;
  float* _re;
  float* _im;
  
  // Prevent uncontrolled usage
  FFTW3FFT(const FFTW3FFT&);
  FFTW3FFT& operator=(const FFTW3FFT&);
};

typedef FFTW3FFT AudioFFT;

#endif // AUDIOFFT_FFTW3_USED


// ================================================================


#ifdef AUDIOFFT_OOURA_USED

/**
* @class OouraFFT
* @brief FFT implementation based on the great radix-4 routines by Takuya Ooura
*/
class OouraFFT : public AudioFFTBase
{
public:
  OouraFFT();
  virtual void init(size_t size);
  virtual void fft(const float* data, float* re, float* im);
  virtual void ifft(float* data, const float* re, const float* im);

private:
  size_t _size;
  std::vector<int> _ip;
  std::vector<double> _w;
  std::vector<double> _buffer;

  // The original FFT routines by Takuya Ooura (see http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html)
  void rdft(int n, int isgn, double *a, int *ip, double *w);
  void makewt(int nw, int *ip, double *w);
  void makect(int nc, int *ip, double *c);
  void bitrv2(int n, int *ip, double *a);
  void cftfsub(int n, double *a, double *w);
  void cftbsub(int n, double *a, double *w);
  void rftfsub(int n, double *a, int nc, double *c);
  void rftbsub(int n, double *a, int nc, double *c);
  void cft1st(int n, double *a, double *w);
  void cftmdl(int n, int l, double *a, double *w);
  
  // Prevent uncontrolled usage
  OouraFFT(const OouraFFT&);
  OouraFFT& operator=(const OouraFFT&);
};

typedef OouraFFT AudioFFT;

#endif // AUDIOFFT_OOURA_USED

} // End of namespace

#endif // Header guard
