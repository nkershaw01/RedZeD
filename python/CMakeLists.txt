cmake_minimum_required(VERSION 3.15...3.31)
project(redzed LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Development.Module (not the full Development) is what a extension module needs,
# and is the component that works when there are no static Python libraries around
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core cpp/bindings.cpp)

target_compile_definitions(_core PRIVATE
    REDZED_NO_MAIN      # drop the command line driver in redzed.cpp
    NDEBUG
)

if(MSVC)
    target_compile_options(_core PRIVATE /O2)
    # windows.h defines min and max as macros, which breaks std::max in redzed.cpp
    target_compile_definitions(_core PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
else()
    target_compile_options(_core PRIVATE -O3)
endif()

install(TARGETS _core DESTINATION redzed_tda)
