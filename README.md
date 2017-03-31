FFTConvolver
============

FFTConvolver is a C++ library for highly efficient convolution of
audio data (e.g. for usage in real-time convolution reverbs etc.).

- Partitioned convolution algorithm (using uniform block sizes)
- Optional support for non-uniform block sizes (TwoStageFFTConvolver)
- No external dependencies (FFT already included)
- Optional optimization for SSE

## Example: FFTConvolver ##

    // Some impulse response consisting of 12345 samples
    const float* impulseResponse = ...;
    const size_t impulseResponseSize = 12345;

    // Setup the convolver object (block size: 128 samples)
    // (typically, this is done in some initialization and not in the audio callback)
    fftconvolver::FFTConvolver convolver;
    convolver.init(128, impulseResponse, impulseResponseSize);

    // Now consecutive audio chunks of arbitrary size can be convolved
    // (typically, this is done in the audio callback)
    convolver.process(audioIn, audioOut, audioLen);
    .
    .
    .
    convolver.process(audioIn, audioOut, audioLen);


## Example: TwoStageFFTConvolver ##

    // Some impulse response consisting of 12345 samples
    const float* impulseResponse = ...;
    const size_t impulseResponseSize = 12345;

    // Setup the convolver object
    // (typically, this is done in some initialization and not in the audio callback)
    fftconvolver::TwoStageFFTConvolver convolver;
    const size_t headBlockSize = 64;
    const size_t tailBlockSize = 1024;
    convolver.init(headBlockSize, tailBlockSize, impulseResponse, impulseResponseSize);

    // Now consecutive audio chunks of arbitrary size can be convolved
    // (typically, this is done in the audio callback)
    convolver.process(audioIn, audioOut, audioLen);
    .
    .
    .
    convolver.process(audioIn, audioOut, audioLen);
