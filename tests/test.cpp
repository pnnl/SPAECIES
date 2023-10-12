#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

#include "spaecies.hpp"

using namespace spaecies;

TEST_CASE( "constructing dimensions using a domain", "[domain]" ) {
  Domain domain;

  SECTION( "a pointer to the returned dimension is kept in the domain" ) {
    DimensionPtr col_dim = domain.add_dimension("column", 8);
    // Note that this is checking for pointer equality, i.e. checking that
    // the returned pointer is to the same location in memory.
    REQUIRE( domain.get_dimension("column") == col_dim );
  }

}

TEST_CASE( "a variable's size can be calculated", "[variable_descriptor]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);
  DimensionPtr lev_dim = domain.add_dimension("level", 10);
  DimensionPtr tracer_dim = domain.add_dimension("tracer", 5);

  SECTION("a variable with no dimensions has a size of 1") {
    VariableDescriptor var_desc("T", Float64Type, {}, "K");
    std::size_t var_size = var_desc.size();
    REQUIRE( var_size == 1 );
  }

  SECTION("a variable with a single dimension has the size of that dimension") {
    VariableDescriptor var_desc("T", Float64Type, {col_dim}, "K");
    std::size_t var_size = var_desc.size();
    REQUIRE( var_size == col_dim->size );
  }

  SECTION("a variable with many dimensions has a size equal to their product") {
    VariableDescriptor var_desc("T", Float64Type,
                                {col_dim, lev_dim, tracer_dim}, "K");
    std::size_t var_size = var_desc.size();
    REQUIRE( var_size == (col_dim->size * lev_dim->size * tracer_dim->size) );
  }

}
