# C++ losc library
add_subdirectory(losc)

# C losc library
add_subdirectory(c_losc)

# python losc library
add_subdirectory(py_losc)

install(DIRECTORY psi4_losc DESTINATION lib/psi4_losc)
