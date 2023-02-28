project(FFTConvolver CXX)

cmake_minimum_required(VERSION 3.12)


option(FFTCONVOLVER_USE_SSE "Enable SSE optimizations" OFF)

# AudioFFT static library
set(AudioFFT_SOURCES
    AudioFFT.cpp
    AudioFFT.h
)
add_library(AudioFFT STATIC ${AudioFFT_SOURCES})

# FFTConvolver static library
set(FFTConvolver_SOURCES
    FFTConvolver.cpp
    FFTConvolver.h

    TwoStageFFTConvolver.cpp
    TwoStageFFTConvolver.h

    Utilities.cpp
    Utilities.h
)

add_library(FFTConvolver STATIC ${FFTConvolver_SOURCES})
target_link_libraries(FFTConvolver PRIVATE AudioFFT)
if(FFTCONVOLVER_USE_SSE)
    target_compile_definitions(FFTConvolver PRIVATE -DFFTCONVOLVER_USE_SSE)
endif()

# Test executable
set(Test_SOURCES
    test/Test.cpp
)

add_executable(Test ${Test_SOURCES})
target_link_libraries(Test PRIVATE FFTConvolver)
