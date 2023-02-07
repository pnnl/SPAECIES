add_library(spaecies STATIC spaecies.cpp)
target_include_directories(spaecies PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
