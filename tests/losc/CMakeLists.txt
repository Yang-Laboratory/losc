set(PROJECT_NAME losc_test)
project(${PROJECT_NAME} C CXX)

file(GLOB SOURCE
    "*.cpp"
    "*.cc")

add_executable(
    ${PROJECT_NAME}
    ${SOURCE})

target_link_libraries(
    ${PROJECT_NAME} PUBLIC
    ${PROJECT_LOSC}
    gtest_main)

add_test(
    NAME ${PROJECT_NAME}
    COMMAND ${PROJECT_NAME})
