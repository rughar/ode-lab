include_guard(GLOBAL)

function(_monorepo_attach_headers TARGET_NAME INC_DIR scope)
  set(options)
  set(oneValueArgs)
  set(multiValueArgs HEADER_GLOBS)
  cmake_parse_arguments(AH "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT AH_HEADER_GLOBS)
    set(AH_HEADER_GLOBS "*.hpp")
  endif()

  set(_hdrs "")
  foreach(glob IN LISTS AH_HEADER_GLOBS)
    file(GLOB_RECURSE _tmp CONFIGURE_DEPENDS "${INC_DIR}/${glob}")
    list(APPEND _hdrs ${_tmp})
  endforeach()

  if(_hdrs)
    target_sources(${TARGET_NAME} ${scope} ${_hdrs})
  endif()
endfunction()


function(monorepo_add_library NAME KIND)
  set(options)
  set(oneValueArgs CXX_STANDARD)
  set(multiValueArgs SOURCES HEADER_GLOBS)
  cmake_parse_arguments(MAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT MAL_CXX_STANDARD)
    set(MAL_CXX_STANDARD 20)
  endif()

  set(LIB_DIR "${CMAKE_CURRENT_LIST_DIR}")
  set(INC_DIR "${LIB_DIR}/include")

  if(KIND STREQUAL "HEADER_ONLY")
    add_library(${NAME} INTERFACE)
    add_library(${NAME}::${NAME} ALIAS ${NAME})

    target_compile_features(${NAME} INTERFACE "cxx_std_${MAL_CXX_STANDARD}")
    target_include_directories(${NAME} INTERFACE
      $<BUILD_INTERFACE:${INC_DIR}>
      $<INSTALL_INTERFACE:include>
    )

    _monorepo_attach_headers(${NAME} "${INC_DIR}" INTERFACE
      HEADER_GLOBS ${MAL_HEADER_GLOBS}
    )

  elseif(KIND STREQUAL "STATIC" OR KIND STREQUAL "SHARED" OR KIND STREQUAL "OBJECT")
    if(NOT MAL_SOURCES)
      message(FATAL_ERROR "monorepo_add_library(${NAME} ${KIND}): you must add SOURCES ...")
    endif()

    add_library(${NAME} ${KIND} ${MAL_SOURCES})
    add_library(${NAME}::${NAME} ALIAS ${NAME})

    target_compile_features(${NAME} PUBLIC "cxx_std_${MAL_CXX_STANDARD}")
    target_include_directories(${NAME} PUBLIC
      $<BUILD_INTERFACE:${INC_DIR}>
      $<INSTALL_INTERFACE:include>
    )

    _monorepo_attach_headers(${NAME} "${INC_DIR}" PRIVATE
      HEADER_GLOBS ${MAL_HEADER_GLOBS}
    )

  else()
    message(FATAL_ERROR "monorepo_add_library: unknown KIND='${KIND}'. Use HEADER_ONLY/STATIC/SHARED/OBJECT.")
  endif()

  install(DIRECTORY "${INC_DIR}/" DESTINATION include)

  if(NOT KIND STREQUAL "OBJECT")
    install(TARGETS ${NAME} EXPORT ${NAME}Targets)
    install(EXPORT ${NAME}Targets
      NAMESPACE ${NAME}::
      DESTINATION "lib/cmake/${NAME}"
    )
  endif()
endfunction()