cmake_minimum_required(VERSION 3.5)
project(Practicum3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp dopri5.cpp sabuhi.cpp maxim.cpp)
add_executable(Practicum3 ${SOURCE_FILES})