add_library(spaecies STATIC
  exceptions.cpp
  dimension.cpp
  domain.cpp
  variable_descriptor.cpp)
target_include_directories(spaecies PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
