message(STATUS "Looking for ISPC...")

find_program(ISPC_COMMAND ispc)
if("${ISPC_COMMAND}" MATCHES "ISPC_COMMAND-NOTFOUND")
  message(STATUS "ISPC Not Found")
else("${ISPC_COMMAND}" MATCHES "ISPC_COMMAND-NOTFOUND")
  message(STATUS "ISPC Found: ${ISPC_COMMAND}")
endif("${ISPC_COMMAND}" MATCHES "ISPC_COMMAND-NOTFOUND")
set(ISPC_FLAGS "" CACHE STRING "ISPC compile flags" )
set(ISPC_FLAGS_DEBUG "-g -O0" CACHE STRING "ISPC debug compile flags")
set(ISPC_FLAGS_RELEASE "-O2 -DNDEBUG" CACHE STRING "ISPC release compile flags")

function(ispc_compile filename flags obj)
  get_filename_component(base ${filename} NAME_WE)
  set(base_abs ${CMAKE_CURRENT_BINARY_DIR}/${base})
  set(base_include_abs ${CMAKE_CURRENT_BINARY_DIR}/include/ispc/${base})

  set(output1 ${base_abs}.o)
  set(output2 ${base_include_abs}.ispc.h)

  if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(ispc_compile_flags "${ISPC_FLAGS} ${ISPC_FLAGS_DEBUG} ${flags}")
  else("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    set(ispc_compile_flags "${ISPC_FLAGS} ${ISPC_FLAGS_RELEASE} ${flags}")
  endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  string(REPLACE " " ";" ispc_compile_flags_list ${ispc_compile_flags})

  add_custom_command(
    OUTPUT ${output1}
    COMMAND ${ISPC_COMMAND} -Iinclude ${ispc_compile_flags_list} ${filename} -o ${base_abs}.o -h ${base_include_abs}.ispc.h
    DEPENDS ${filename})

  add_custom_command(
    OUTPUT ${output2}
    COMMAND 
    DEPENDS ${output1})

  set_source_files_properties(${output1} PROPERTIES GENERATED TRUE)
  set_source_files_properties(${output2} PROPERTIES GENERATED TRUE)

  set(${obj} ${base_abs}.o PARENT_SCOPE)
endfunction()

