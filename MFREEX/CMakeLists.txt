file(GLOB_RECURSE MFREEX_SRCS    CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE MFREEX_HEADERS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

add_library(MFREEX STATIC ${MFREEX_SRCS} ${MFREEX_HEADERS})
add_library(MFREEX::MFREEX ALIAS MFREEX)

target_include_directories(MFREEX PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MFREEX PUBLIC
  CVLIB
  EFSFLIB
  triangulate
  util
  eig3
  gpc
  thread_pool
  meschach
  mfreex_compat)

if(MFREEX_ENABLE_WARNINGS AND NOT MSVC)
  target_compile_options(MFREEX PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()
if(MFREEX_WARNINGS_AS_ERRORS AND NOT MSVC)
  target_compile_options(MFREEX PRIVATE -Werror)
endif()
