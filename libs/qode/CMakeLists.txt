add_library(qode INTERFACE)

target_include_directories(qode INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(qode INTERFACE math)