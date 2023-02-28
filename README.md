FFTConvolver
============

FFTConvolver is a C++ library for highly efficient convolution of
audio data (e.g. for usage in real-time convolution reverbs etc.).

- Partitioned convolution algorithm (using uniform block sizes)
- Optional support for non-uniform block sizes (TwoStageFFTConvolver)
- No external dependencies (FFT already included)
- Optional optimization for SSE (enabled by defining FFTCONVOLVER_USE_SSE)

# Building

Use [CMake](https://cmake.org) to generate the project, e.g.:

```sh
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

# Or enable SSE optimizations (OFF by default)
cmake .. -DCMAKE_BUILD_TYPE=Release -DFFTCONVOLVER_USE_SSE=ON
```

You can build the project by using the generated build files directly, or by
running the build tool through CMake:

```sh
cmake --build .
```

## Running the Tests

The build produces a `Test` executable and a `FFTConvolver` static library.
Simply run the executable as follows:

```sh
./Test
```
