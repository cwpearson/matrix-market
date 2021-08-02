# Copyright (C) 2021 Carl Pearson
# This code is released under the GPLv3 license

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native CXX_HAS_MARCH)
check_cxx_compiler_flag(-mcpu=native CXX_HAS_MCPU)

function(mm_tool_properties target)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD 11)
    set_property(TARGET ${target} PROPERTY CXX_EXTENSIONS OFF)
    set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)
endfunction()

function(mm_tool_options target)
    target_compile_options(${target} PRIVATE
        -ffast-math;
        -Wall
        -Wextra
        -Wcast-align;
        -Wdisabled-optimization;
        -Wformat=2;
        -Winit-self;
        -Wmissing-include-dirs;
        -Woverloaded-virtual;
        -Wpointer-arith;
        -Wshadow;
        -Wstrict-aliasing;
        -Wswitch-enum;
        -Wvla;
    )
    # prefer march=native, but fallback to mcpu=native if needed
    if (CXX_HAS_MARCH)
        message(STATUS "-march=native detected")
        target_compile_options(
        ${target}
        PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-march=native>
        )
    elseif (CXX_HAS_MCPU)
        message(STATUS "-mcpu=native -mtune=native instead of -march=native")
        target_compile_options(
        ${target}
        PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-mcpu=native -mtune=native>
        )
    endif()
endfunction()


add_executable(mtx-stats mtx_stats.cpp)
mm_tool_properties(mtx-stats)
mm_tool_options(mtx-stats)
target_include_directories(mtx-stats PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(mtx-stats mm)