cmake_minimum_required(VERSION 3.6)

file(GLOB_RECURSE SRC_LIST "*.h" "*.cpp")

message(STATUS ${SRC_LIST})

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME})
