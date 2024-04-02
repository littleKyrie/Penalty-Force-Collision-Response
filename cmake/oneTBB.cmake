if(TARGET TBB::tbb)
    return()
endif ()

include(FetchContent)
set(FETCHCONTENT_SOURCE_DIR_TBB ${PROJECT_SOURCE_DIR}/external/oneTBB)
FetchContent_Declare(
        oneTBB
        SOURCE_DIR ${FETCHCONTENT_SOURCE_DIR_TBB}
)
FetchContent_MakeAvailable(oneTBB)