file(GLOB UTIL_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
add_library(util STATIC ${UTIL_SRCS})
target_include_directories(util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(util PUBLIC meschach)
