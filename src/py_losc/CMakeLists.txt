set(PROJECT_NAME py_losc_core)
project(${PROJECT_NAME} C CXX)

# sources
file(GLOB SOURCES "*.cpp")

pybind11_add_module(
    ${PROJECT_NAME} SHARED
    ${SOURCES})

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${PROJECT_LOSC})

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(py_losc.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# ==> install py_losc module <==
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})
install(FILES __init__.py py_losc.py DESTINATION lib/${PROJECT_NAME})
