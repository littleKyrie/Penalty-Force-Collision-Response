if (TBB_FOUND)
    return()
endif ()

find_path(TBB_INCLUDE_DIR /tbb.h
        PATHS ${CMAKE_SOURCE_DIR}/extern/oneTBB/include/oneapi
)

# list(APPEND CMAKE_MODULE_PATH "/Users/kyriezhang/Public/code/c++/test-FEM-3d/extern/oneTBB/lib/cmake/TBB")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/extern/oneTBB/lib/cmake/TBB")
include(TBBConfig)

# Provide TBB_FOUND variable
set(TBB_FOUND TRUE)