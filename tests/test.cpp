#define CATCH_CONFIG_MAIN
#include "catch2/catch_test_macros.hpp"

#include "spaecies.hpp"

using namespace spaecies;

TEST_CASE( "constructing dimensions using a domain", "[Domain]" ) {
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

TEST_CASE( "constructing variables using a domain", "[Domain]" ) {
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

TEST_CASE( "a variable's size can be calculated", "[VariableDescriptor]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);
  DimensionPtr lev_dim = domain.add_dimension("level", 10);
  DimensionPtr tracer_dim = domain.add_dimension("tracer", 5);

  SECTION( "a variable with no dimensions has a size of 1" ) {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type, {}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == 1 );
  }

  SECTION( "a variable with a single dimension has the size of that dimension" ) {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == col_dim->size );
  }

  SECTION( "a variable with many dimensions has a size equal to their product" ) {
    VarDescPtr var_desc = domain.add_var_desc("T", Float64Type,
                                              {col_dim, lev_dim, tracer_dim}, "K");
    std::size_t var_size = var_desc->size();
    REQUIRE( var_size == (col_dim->size * lev_dim->size * tracer_dim->size) );
  }

}

TEST_CASE( "accessing values of a contiguous variable", "[ContiguousVariable]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);
  DimensionPtr lev_dim = domain.add_dimension("level", 10);
  DimensionPtr tracer_dim = domain.add_dimension("tracer", 5);

  SECTION( "a scalar variable's value can be accessed" ) {
    VarDescPtr boil_desc = domain.add_var_desc("BOILING_POINT", Float64Type, {}, "K");
    double boiling_point = 373.16;
    ContiguousVariable<double> t_boil(boil_desc, &boiling_point);
    REQUIRE( t_boil.size() == boil_desc->size() );
    REQUIRE( t_boil[0] == boiling_point );
    // Pretend we have a reason to change this. "Oops, we had the wrong value, fix it."
    t_boil[0] = 373.15;
    REQUIRE( t_boil[0] == 373.15 );
  }

  SECTION( "a vector variable's values can be accessed" ) {
    VarDescPtr t_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    double t_buff[col_dim->size];
    for (int i; i != col_dim->size; ++i) {
      t_buff[i] = 2. * i;
    }
    ContiguousVariable<double> t_var(t_desc, t_buff);
    REQUIRE( t_var.size() == t_desc->size() );
    for (int i; i != col_dim->size; ++i) {
      REQUIRE( t_var[i] == 2. * i );
      t_var[i] = t_var[i] + 3.;
      REQUIRE( t_var[i] == 2. * i + 3. );
    }
  }

  SECTION( "a multidimensional variable's flattened representation can be accessed" ) {
    VarDescPtr t_desc = domain.add_var_desc("T", Float64Type, {col_dim, lev_dim}, "K");
    std::size_t total_size = col_dim->size * lev_dim->size;
    double t_buff[total_size];
    for (int i; i != total_size; ++i) {
      t_buff[i] = 2. * i;
    }
    ContiguousVariable<double> t_var(t_desc, t_buff);
    REQUIRE( t_var.size() == t_desc->size() );
    for (int i; i != total_size; ++i) {
      REQUIRE( t_var[i] == 2. * i );
      t_var[i] = t_var[i] + 3.;
      REQUIRE( t_var[i] == 2. * i + 3. );
    }
  }

}
