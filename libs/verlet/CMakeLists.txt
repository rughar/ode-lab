add_library(verlet INTERFACE)

target_include_directories(verlet INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(verlet INTERFACE math)