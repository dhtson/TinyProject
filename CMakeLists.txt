cmake_minimum_required(VERSION 3.11)
project(TinyProject VERSION 1.0 LANGUAGES CXX)

aux_source_directory(src MAIN_SRC)
add_executable(${PROJECT_NAME} ${MAIN_SRC})

set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /utf-8 /W1)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
endif()

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/data/machine.data $<TARGET_FILE_DIR:${PROJECT_NAME}>
)
