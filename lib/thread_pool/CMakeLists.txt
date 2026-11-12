# C-Thread-Pool (https://github.com/Pithikos/C-Thread-Pool)
add_library(thread_pool STATIC thpool.c)
target_include_directories(thread_pool PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(thread_pool PUBLIC Threads::Threads mfreex_compat)
