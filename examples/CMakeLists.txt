add_executable(readme readme.cpp)

target_compile_options(readme PRIVATE
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

target_include_directories(readme PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(readme mm)