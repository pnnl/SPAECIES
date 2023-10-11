#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

#include "spaecies.hpp"

using namespace spaecies;

TEST_CASE( "constructing dimensions using a domain", "[domain]" ) {
  Domain dom;

  SECTION( "a pointer to the returned dimension is kept in the domain" ) {
    DimensionPtr col_dim = dom.add_dimension("column", 8);
    // Note that this is checking for pointer equality, i.e. checking that
    // the returned pointer is to the same location in memory.
    REQUIRE( dom.get_dimension("column") == col_dim );
  }

}

TEST_CASE( "a variable's size can be calculated", "[variable_descriptor]" ) {
  DimensionPtr col_dim_ptr(new Dimension("column", 8));
  Dimension lev_dim("level", 10);
  Dimension tracer_dim("tracer", 5);

  SECTION("a variable with a single dimension has the size of that dimension") {
    VariableDescriptor var_desc("T", Float64Type, {col_dim_ptr}, "K");
    std::size_t var_size = var_desc.size();
    REQUIRE( var_size == col_dim_ptr->size );
  }

}
