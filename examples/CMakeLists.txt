# Copyright (C) 2021 Carl Pearson
# This code is released under the GPLv3 license

add_executable(readme readme.cpp)

target_include_directories(readme PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(readme mm)