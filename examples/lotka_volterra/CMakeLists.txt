file(GLOB_RECURSE SANDBOX_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_executable(lotka_volterra ${SANDBOX_SOURCES})

target_link_libraries(lotka_volterra PRIVATE verlet)