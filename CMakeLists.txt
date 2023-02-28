project(FFTConvolver CXX)

cmake_minimum_required(VERSION 3.12)


set(SOURCES
    AudioFFT.cpp
    AudioFFT.h

    FFTConvolver.cpp
    FFTConvolver.h

    TwoStageFFTConvolver.cpp
    TwoStageFFTConvolver.h

    Utilities.cpp
    Utilities.h
)

add_library(FFTConvolver STATIC ${SOURCES})

set(Test_SOURCES
    test/Test.cpp
)

add_executable(Test ${Test_SOURCES})
target_link_libraries(Test PRIVATE FFTConvolver)
