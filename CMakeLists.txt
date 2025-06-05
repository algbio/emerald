cmake_minimum_required(VERSION 3.9)
# Get project name from folder name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
# all .cpp files
set(SOURCES src/main.cpp src/alpha_safe_paths.cpp src/optimal_paths.cpp src/safety_windows.cpp src/draw_subgraph.cpp src/write_json.cpp src/wasm_interface.cpp)
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
    
    # Export both runtime methods and the generate_alignment_json function
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sEXPORTED_RUNTIME_METHODS=callMain,FS,UTF8ToString,stringToUTF8")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sEXPORTED_FUNCTIONS=_main,_generate_alignment_json,_malloc,_free")
    
    # Add memory settings for large protein sequences
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s INITIAL_MEMORY=512MB")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s MAXIMUM_MEMORY=2GB")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s ALLOW_MEMORY_GROWTH=1")
    
    # Create JavaScript bindings for easier use from JS
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/emerald_bindings.js "
// JavaScript bindings for Emerald WASM module
Module.onRuntimeInitialized = function() {
  console.log('WASM runtime initialized, binding functions...');
  
  // Create user-friendly wrapper for generate_alignment_json
  Module.generateAlignmentJson = function(refSeq, refDesc, memSeq, memDesc, alpha = 0.75, delta = 0, gapCost = -1, startGap = -11) {
    // Allocate memory for string parameters
    var refSeqBuf = Module._malloc(refSeq.length + 1);
    var refDescBuf = Module._malloc(refDesc.length + 1);
    var memSeqBuf = Module._malloc(memSeq.length + 1);
    var memDescBuf = Module._malloc(memDesc.length + 1);
    
    // Write strings to memory
    Module.stringToUTF8(refSeq, refSeqBuf, refSeq.length + 1);
    Module.stringToUTF8(refDesc, refDescBuf, refDesc.length + 1);
    Module.stringToUTF8(memSeq, memSeqBuf, memSeq.length + 1);
    Module.stringToUTF8(memDesc, memDescBuf, memDesc.length + 1);
    
    try {
      // Call the C++ function
      var resultPtr = Module._generate_alignment_json(refSeqBuf, refDescBuf, memSeqBuf, memDescBuf, alpha, delta, gapCost, startGap);
      
      // Convert result to JavaScript string
      var result = Module.UTF8ToString(resultPtr);
      return result;
    } finally {
      // Free allocated memory
      Module._free(refSeqBuf);
      Module._free(refDescBuf);
      Module._free(memSeqBuf);
      Module._free(memDescBuf);
    }
  };
  
  console.log('Emerald WASM bindings loaded successfully!');
};
")
    
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --post-js ${CMAKE_CURRENT_BINARY_DIR}/emerald_bindings.js")
endif()
