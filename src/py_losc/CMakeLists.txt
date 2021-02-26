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

file(COPY __init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY py_losc.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# ==> install py_losc module <==
install(TARGETS ${PROJECT_NAME} DESTINATION liblosc/py_losc)
install(FILES __init__.py py_losc.py DESTINATION liblosc/py_losc)
