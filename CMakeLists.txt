cmake_minimum_required(VERSION 3.9)
# Get project name from folder name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
# all .cpp files
set(SOURCES src/main.cpp src/alpha_safe_paths.cpp src/optimal_paths.cpp src/safety_windows.cpp src/draw_subgraph.cpp src/write_json.cpp)
add_executable(${ProjectId} ${SOURCES})
# compiler flags for debugging
if (APPLE)
    set(CMAKE_CXX_FLAGS "-O3 -openmp")
elseif(EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS "-O3 -s WASM=1 -s USE_SDL=2")
    # Uncomment the next line ONLY if you have pthreads and OpenMP support in Emscripten
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_PTHREADS=1 -fopenmp -s PTHREAD_POOL_SIZE=4")
else()
    set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
endif()

# GMP
if(NOT EMSCRIPTEN)
    find_package(GMP REQUIRED)
    include_directories(./ ${GMP_INCLUDE_DIR})
    target_link_libraries(${ProjectId} ${GMP_LIBRARY})
else()
    # Use gmp-wasm for Emscripten builds
    include_directories($ENV{HOME}/opt-wasm/include)
    target_link_libraries(${ProjectId} $ENV{HOME}/opt-wasm/lib/libgmp.a)
endif()

# OpenMP
#find_package(OpenMP REQUIRED)
#target_link_libraries(${ProjectId} OpenMP::OpenMP_CXX)
# runnable output file
set_target_properties(${ProjectId} PROPERTIES OUTPUT_NAME "emerald")

if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WASM=1 -s USE_SDL=2")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sEXPORTED_RUNTIME_METHODS=callMain,FS")
endif()
