add_library(rkgl INTERFACE)

target_include_directories(rkgl INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(rkgl INTERFACE math)