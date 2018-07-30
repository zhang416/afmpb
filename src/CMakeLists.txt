FILE(GLOB sources *.cc)
add_executable(dafmpb ${sources})

include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${dashmm_SOURCE_DIR}/include
)

target_link_libraries(dafmpb PUBLIC dashmm ${HPX_LDFLAGS})

