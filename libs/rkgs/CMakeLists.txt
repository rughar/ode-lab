add_library(rkgs INTERFACE)

target_include_directories(rkgs INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(rkgs INTERFACE math)