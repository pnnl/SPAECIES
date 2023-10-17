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

  SECTION( "multiple dimensions can be queried in order" ) {
    DimensionPtr col_dim = domain.add_dimension("column", 8);
    DimensionPtr lev_dim = domain.add_dimension("level", 10);
    // Deliberately requesting these out of order to ensure any future refactoring
    // will cause the test to fail unless it respects the requested order.
    std::vector<DimensionPtr> req_dims = domain.get_dimensions({"level", "column"});
    REQUIRE( req_dims.size() == 2 );
    REQUIRE( req_dims[0] == lev_dim );
    REQUIRE( req_dims[1] == col_dim );
  }

}

TEST_CASE( "constructing variables using a domain", "[domain]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);

  SECTION( "a pointer to the returned variable is kept in the domain" ) {
    VarDescPtr temp_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    // Note that this is checking for pointer equality, i.e. checking that
    // the returned pointer is to the same location in memory.
    REQUIRE( domain.get_var_desc("T") == temp_desc );
  }

  SECTION( "multiple variable descriptors can be queried in order" ) {
    VarDescPtr temp_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    VarDescPtr q_desc = domain.add_var_desc("q", Float64Type, {col_dim}, "kg/kg");
    // Deliberately requesting these out of order to ensure any future refactoring
    // will cause the test to fail unless it respects the requested order.
    std::vector<VarDescPtr> req_vars = domain.get_var_descs({"q", "T"});
    REQUIRE( req_vars.size() == 2 );
    REQUIRE( req_vars[0] == q_desc );
    REQUIRE( req_vars[1] == temp_desc );
  }

}

TEST_CASE( "a variable's size can be calculated", "[variable_descriptor]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);
  DimensionPtr lev_dim = domain.add_dimension("level", 10);
  DimensionPtr tracer_dim = domain.add_dimension("tracer", 5);

  SECTION("a variable with no dimensions has a size of 1") {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type, {}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == 1 );
  }

  SECTION("a variable with a single dimension has the size of that dimension") {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == col_dim->size );
  }

  SECTION("a variable with many dimensions has a size equal to their product") {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type,
                                              {col_dim, lev_dim, tracer_dim}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == (col_dim->size * lev_dim->size * tracer_dim->size) );
  }

}
