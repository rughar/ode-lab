add_executable(ode_lab_utest
  utest_run.cpp
)

target_link_libraries(ode_lab_utest PRIVATE
  qode
  math
)

target_include_directories(ode_lab_utest PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)