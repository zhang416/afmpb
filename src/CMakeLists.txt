FILE(GLOB sources *.cc)
add_executable(afmpb ${sources})

include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/contrib/dashmm/include
  ${PROJECT_SOURCE_DIR}/include
)

link_directories(${PROJECT_SOURCE_DIR}/contrib/dashmm/lib)

target_link_libraries(afmpb PUBLIC dashmm ${HPX_LDFLAGS})

