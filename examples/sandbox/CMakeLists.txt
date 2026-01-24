file(GLOB_RECURSE SANDBOX_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(sandbox ${SANDBOX_SOURCES})

target_link_libraries(sandbox PRIVATE verlet)