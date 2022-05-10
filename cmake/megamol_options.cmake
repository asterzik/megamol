# C++ standard
set(CMAKE_CXX_STANDARD 17)

# Word size
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(BITS 64)
else()
  set(BITS 32)
endif()

# Compiler flags (inspired by OSPRay build)
option(DISABLE_WARNINGS "Disables all compiler warnings" ON)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(MEGAMOL_COMPILER_GCC TRUE)
  include(gcc)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(MEGAMOL_COMPILER_CLANG TRUE)
  include(clang)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(MEGAMOL_COMPILER_MSVC TRUE)
  include(msvc)
else()
  message(FATAL_ERROR
    "Unsupported compiler specified: '${CMAKE_CXX_COMPILER_ID}'")
endif()

# OpenGL
option(ENABLE_GL "Enable GL support" ON)
if (ENABLE_GL)
  add_compile_definitions(WITH_GL)
  find_package(OpenGL REQUIRED)
endif()

# CUDA
option(ENABLE_CUDA "Enable CUDA support" OFF)
if(ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_ARCHITECTURES FALSE)
endif()

# MPI
option(ENABLE_MPI "Enable MPI support" OFF)
set(MPI_GUESS_LIBRARY_NAME "undef" CACHE STRING "Override MPI library name, e.g., MSMPI, MPICH2")
if(ENABLE_MPI)
  if(MPI_GUESS_LIBRARY_NAME STREQUAL "undef")
    message(FATAL_ERROR "you must set MPI_GUESS_LIBRARY_NAME to ovveride automatic finding of unwanted MPI libraries (or empty for default)")
  endif()
  find_package(MPI REQUIRED)
  if(MPI_C_FOUND)
    target_compile_definitions(MPI::MPI_C INTERFACE "-DWITH_MPI")
  endif()
endif()

# GLFW
option(USE_GLFW "Use GLFW" ON)

# Profiling
if (ENABLE_GL)
  option(ENABLE_PROFILING "Enable profiling code" OFF)
  if (ENABLE_PROFILING)
    add_compile_definitions(PROFILING)
  endif()
endif()

# VR Service / mwk-mint, interop, Spout2
if (ENABLE_GL)
  option(ENABLE_VR_SERVICE_UNITY_KOLABBW "Enable KolabBW-Unity-Interop in VR Service" OFF)
  if(ENABLE_VR_SERVICE_UNITY_KOLABBW)
    add_compile_definitions(WITH_VR_SERVICE_UNITY_KOLABBW)
  endif()
endif()

# CUE
if (WIN32)
  option(ENABLE_CUESDK "Enable CUE for highlighting hotkeys on Corsair Keyboards" OFF)
  if (ENABLE_CUESDK)
    add_compile_definitions(CUESDK_ENABLED)
  endif()
endif()

# Threading (XXX: this is a bit wonky due to Ubuntu/clang)
include(CheckFunctionExists)
check_function_exists(pthread_create HAVE_PTHREAD)
if(HAVE_PTHREAD)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads REQUIRED)
endif()

# OpenMP
if(UNIX)
  find_package(OpenMP REQUIRED)
endif()
if(OPENMP_FOUND OR WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
