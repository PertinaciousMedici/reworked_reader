cmake_minimum_required(VERSION 4.2)
project(NovoLeitor)

set(CMAKE_CXX_STANDARD 23)

add_executable(NovoLeitor src/main.cpp
        src/instance.h
        src/parser.cpp
        src/parser.h
)
