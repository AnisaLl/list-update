cmake_minimum_required(VERSION 3.10)
project(list_update_final)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost COMPONENTS system filesystem REQUIRED)

add_executable(list_update_final main.cpp)

target_link_libraries(list_update_final
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        )


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")