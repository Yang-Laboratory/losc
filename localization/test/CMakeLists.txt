cmake_minimum_required(VERSION 3.8)

set(project_name localization_test)
project(${project_name} C CXX)

file(GLOB Sources
    "*.cpp"
    "*.cc")
file(GLOB Headers
    "*.hpp"
    "*.h")

add_executable(
    ${project_name}
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name} PUBLIC
    localization
    gtest_main)

add_test(
    NAME ${project_name}
    COMMAND ${project_name})
