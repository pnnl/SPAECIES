add_executable(tests test.cpp)
target_link_libraries(tests
  Catch2::Catch2WithMain
  spaecies)

catch_discover_tests(tests)
