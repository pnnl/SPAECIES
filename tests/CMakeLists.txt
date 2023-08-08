add_executable(tests test.cpp)
target_link_libraries(tests Catch2::Catch2WithMain)

catch_discover_tests(tests)
