cmake_minimum_required(VERSION 2.6)
project(DataStructure)

include_directories("include/")
file(GLOB SOURCES "src/*.cpp")

add_executable(a.exe ${SOURCES})
