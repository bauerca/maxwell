include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(AnasaziInterface AnasaziInterface.cpp)
add_dependencies(AnasaziInterface maxwell)
target_link_libraries(AnasaziInterface mxwll
  ${MAXWELL_EXTPROJ_LINK_LIBS} blas lapack z m)

add_executable(RCPTest RCPTest.cpp)
add_dependencies(RCPTest maxwell)
target_link_libraries(RCPTest mxwll
  ${MAXWELL_EXTPROJ_LINK_LIBS} blas lapack z m)

add_test(AnasaziInterface AnasaziInterface)
