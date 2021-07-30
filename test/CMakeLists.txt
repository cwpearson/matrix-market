add_executable(test-cpu
test.cpp)

target_compile_options(test-cpu PRIVATE
 -Wall
 -Wextra
 -Wcast-align;
 -Wdisabled-optimization;
 -Wformat=2;
 -Winit-self;
 -Wlogical-op'
 -Wmissing-include-dirs;
 -Woverloaded-virtual;
 -Wpointer-arith;
 -Wshadow;
 -Wstrict-aliasing;
 -Wswitch-enum;
 -Wvla;
)

target_include_directories(test-cpu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(test-cpu mm)
add_test(NAME test-cpu COMMAND test-cpu "${CMAKE_CURRENT_SOURCE_DIR}/data")