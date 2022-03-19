set(PROJECT_NAME pyscf_losc)

# sources
file(GLOB SOURCES "*.py")
file(COPY ${SOURCES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# ==> install py_losc module <==
install(FILES ${SOURCES} DESTINATION liblosc/pyscf_losc)

