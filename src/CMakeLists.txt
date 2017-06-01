cmake_minimum_required(VERSION 3.6)

file(GLOB_RECURSE SRC_LIST "*.h" "*.hpp" "*.cpp")

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME})
