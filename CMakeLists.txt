cmake_minimum_required(VERSION 3.10)
project(terrier)

set(CMAKE_CXX_STANDARD 14)
find_package(Boost COMPONENTS unit_test_framework)
if (Boost_FOUND)
    message("Boost is found")
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(test_terrier test_simulation.cpp simulation.cpp simulation.h)
    target_link_libraries(test_terrier ${Boost_LIBRARIES})
    enable_testing()
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_test(test_terrier test_terrier)
else ()
    message(WARNING "Boost unit test framework not found, building without test suite.
To enable test suite, please install boost")
endif ()
add_executable(terrier main.cpp simulation.cpp simulation.h)
add_executable(examinations examinations_runs.cpp simulation.cpp simulation.h)

add_executable(single_channel single_channel_runs.cpp simulation.cpp simulation.h)
add_executable(double_channel double_channel_runs.cpp simulation.cpp simulation.h)
target_link_libraries(terrier)
target_link_libraries(examinations)
target_link_libraries(single_channel)
target_link_libraries(double_channel)

add_custom_command(TARGET single_channel POST_BUILD COMMAND ${CMAKE_SOURCE_DIR}/populate.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")
add_custom_command(TARGET double_channel POST_BUILD COMMAND ${CMAKE_SOURCE_DIR}/populate.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")
add_custom_command(TARGET terrier POST_BUILD COMMAND ${CMAKE_SOURCE_DIR}/populate.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")
add_custom_command(TARGET examinations POST_BUILD COMMAND ${CMAKE_SOURCE_DIR}/populate.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")


