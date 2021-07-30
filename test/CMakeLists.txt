add_executable(test-cpu
test.cpp)



target_compile_options(test-cpu PRIVATE
 -Wall
 -Wextra
 -Wcast-align;
 -Wdisabled-optimization;
 -Wformat=2;
 -Winit-self;
 -Wmissing-include-dirs;
 -Woverloaded-virtual;
 -Wpointer-arith;
 -Wpedantic;
 -Wshadow;
 -Wstrict-aliasing;
 -Wswitch-enum;
 -Wvla;
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(test-cpu PRIVATE -Wlogical-op)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif()

target_include_directories(test-cpu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(test-cpu mm)
add_test(NAME test-cpu COMMAND test-cpu "${CMAKE_CURRENT_SOURCE_DIR}/data")