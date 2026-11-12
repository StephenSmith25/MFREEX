add_library(gpc STATIC gpc.c)
target_include_directories(gpc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gpc PUBLIC mfreex_compat)
