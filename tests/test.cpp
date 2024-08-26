#define CATCH_CONFIG_MAIN
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_all.hpp"
#include "catch2/generators/catch_generators_all.hpp"

#include "spaecies.hpp"

using namespace spaecies;
using Catch::Matchers::ContainsSubstring, Catch::Matchers::MessageMatches;

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

  SECTION( "querying a non-existent dimension raises an error" ) {
    // Ensure that an exception is thrown with the requested dimension
    // in the what string.
    REQUIRE_THROWS_MATCHES( domain.get_dimension("column"), DimensionNotFoundException, MessageMatches(ContainsSubstring("column") && ContainsSubstring("domain")) );
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

  SECTION( "querying a non-existent variable raises an error" ) {
    // Ensure that an exception is thrown with the requested variable
    // in the what string.
    REQUIRE_THROWS_MATCHES( domain.get_var_desc("qr"), VariableNotFoundException, MessageMatches(ContainsSubstring("qr") && ContainsSubstring("domain")) );
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

TEST_CASE( "asking for human-readable type name", "[VariableDescriptor]" ) {

  SECTION( "correct type names for supported types" ) {
    REQUIRE( spaecies_type_name(Float64Type) == "64-bit float" );
    REQUIRE( spaecies_type_name(Float32Type) == "32-bit float" );
    REQUIRE( spaecies_type_name(Int64Type) == "64-bit integer" );
    REQUIRE( spaecies_type_name(Int32Type) == "32-bit integer" );
    REQUIRE( spaecies_type_name(BoolType) == "boolean" );
    REQUIRE( spaecies_type_name(InvalidType) == "[invalid type]" );
  }

}

TEST_CASE( "raising exceptions when buffer and descriptor type don't match", "[TypeMismatchException]" ) {
  Domain domain;
  double var_double;
  float var_float;
  std::int64_t var_int64;
  std::int32_t var_int32;
  bool var_bool;

  // For this and the following sections, we don't require the exception to use an exact string,
  // but we do check for specific substrings communicating particular information.
  SECTION( "a descriptor only works with buffers of the same type" ) {
    auto desc_type = GENERATE( Float64Type, Float32Type, Int64Type, Int32Type, BoolType, InvalidType );
    VarDescPtr desc = domain.add_var_desc("foo_variable", desc_type, {}, "");
    auto match_desc_type = ContainsSubstring(spaecies_type_name(desc_type));
    auto match_desc_name = ContainsSubstring("foo_variable");
    if (desc_type == Float64Type) {
      REQUIRE_NOTHROW( ContiguousVariableView(desc, &var_double) );
    }
    else {
      auto match_var_type = ContainsSubstring(spaecies_type_name(Float64Type));
      auto match =  MessageMatches(match_desc_type && match_desc_name
                                   && match_var_type);
      REQUIRE_THROWS_MATCHES( ContiguousVariableView(desc, &var_double), TypeMismatchException,
                              match );
    }
    if (desc_type == Float32Type) {
      REQUIRE_NOTHROW( ContiguousVariableView(desc, &var_float) );
    }
    else {
      auto match_var_type = ContainsSubstring(spaecies_type_name(Float32Type));
      auto match =  MessageMatches(match_desc_type && match_desc_name
                                   && match_var_type);
      REQUIRE_THROWS_MATCHES( ContiguousVariableView(desc, &var_float), TypeMismatchException,
                              match );
    }
    if (desc_type == Int64Type) {
      REQUIRE_NOTHROW( ContiguousVariableView(desc, &var_int64) );
    }
    else {
      auto match_var_type = ContainsSubstring(spaecies_type_name(Int64Type));
      auto match =  MessageMatches(match_desc_type && match_desc_name
                                   && match_var_type);
      REQUIRE_THROWS_MATCHES( ContiguousVariableView(desc, &var_int64), TypeMismatchException,
                              match );
    }
    if (desc_type == Int32Type) {
      REQUIRE_NOTHROW( ContiguousVariableView(desc, &var_int32) );
    }
    else {
      auto match_var_type = ContainsSubstring(spaecies_type_name(Int32Type));
      auto match =  MessageMatches(match_desc_type && match_desc_name
                                   && match_var_type);
      REQUIRE_THROWS_MATCHES( ContiguousVariableView(desc, &var_int32), TypeMismatchException,
                              match );
    }
    if (desc_type == BoolType) {
      REQUIRE_NOTHROW( ContiguousVariableView(desc, &var_bool) );
    }
    else {
      auto match_var_type = ContainsSubstring(spaecies_type_name(BoolType));
      auto match =  MessageMatches(match_desc_type && match_desc_name
                                   && match_var_type);
      REQUIRE_THROWS_MATCHES( ContiguousVariableView(desc, &var_bool), TypeMismatchException,
                              match );
    }
  }

}

TEST_CASE( "accessing values of a contiguous variable", "[ContiguousVariableView]" ) {
  Domain domain;
  DimensionPtr col_dim = domain.add_dimension("column", 8);
  DimensionPtr lev_dim = domain.add_dimension("level", 10);
  DimensionPtr tracer_dim = domain.add_dimension("tracer", 5);

  SECTION( "a scalar variable's value can be accessed" ) {
    VarDescPtr boil_desc = domain.add_var_desc("BOILING_POINT", Float64Type, {}, "K");
    double boiling_point = 373.16;
    ContiguousVariableView<double> t_boil(boil_desc, &boiling_point);
    REQUIRE( t_boil.size() == boil_desc->size() );
    REQUIRE( t_boil[0] == boiling_point );
    // Pretend we have a reason to change this. "Oops, we had the wrong value, fix it."
    t_boil[0] = 373.15;
    REQUIRE( t_boil[0] == 373.15 );
  }

  SECTION( "make_contiguous_variable_view can also create views" ) {
    VarDescPtr boil_desc = domain.add_var_desc("BOILING_POINT", Float64Type, {}, "K");
    double boiling_point = 373.16;
    ContiguousVariableView<double> t_boil = make_contiguous_variable_view(boil_desc, &boiling_point);
    REQUIRE( t_boil.size() == boil_desc->size() );
    REQUIRE( t_boil[0] == boiling_point );
    // Pretend we have a reason to change this. "Oops, we had the wrong value, fix it."
    t_boil[0] = 373.15;
    REQUIRE( t_boil[0] == 373.15 );
  }

  SECTION( "make_contiguous_variable_view can create const views" ) {
    VarDescPtr boil_desc = domain.add_var_desc("BOILING_POINT", Float64Type, {}, "K");
    double boiling_point = 373.16;
    const double* const_boil_ptr = &boiling_point;
    const ContiguousVariableView<double> t_boil = make_contiguous_variable_view(boil_desc, const_boil_ptr);
    REQUIRE( t_boil.size() == boil_desc->size() );
    REQUIRE( t_boil[0] == boiling_point );
    // External code changing pointed-to memory.
    boiling_point = 373.15;
    REQUIRE( t_boil[0] == 373.15 );
  }

  SECTION( "a vector variable's values can be accessed" ) {
    VarDescPtr t_desc = domain.add_var_desc("T", Float64Type, {col_dim}, "K");
    double t_buff[col_dim->size];
    for (int i = 0; i != col_dim->size; ++i) {
      t_buff[i] = 2. * i;
    }
    ContiguousVariableView<double> t_var(t_desc, t_buff);
    REQUIRE( t_var.size() == t_desc->size() );
    for (int i = 0; i != col_dim->size; ++i) {
      REQUIRE( t_var[i] == 2. * i );
      t_var[i] = t_var[i] + 3.;
      REQUIRE( t_var[i] == 2. * i + 3. );
    }
  }

  SECTION( "a multidimensional variable's flattened representation can be accessed" ) {
    VarDescPtr t_desc = domain.add_var_desc("T", Float64Type, {col_dim, lev_dim}, "K");
    std::size_t total_size = col_dim->size * lev_dim->size;
    double t_buff[total_size];
    for (int i = 0; i != total_size; ++i) {
      t_buff[i] = 2. * i;
    }
    ContiguousVariableView<double> t_var(t_desc, t_buff);
    REQUIRE( t_var.size() == t_desc->size() );
    for (int i = 0; i != total_size; ++i) {
      REQUIRE( t_var[i] == 2. * i );
      t_var[i] = t_var[i] + 3.;
      REQUIRE( t_var[i] == 2. * i + 3. );
    }
  }

}

TEST_CASE( "interacting with data through an array view", "[VariableArrayView]" ) {
  Domain dom;
  DimensionPtr col_dim = dom.add_dimension("column", 8);
  DimensionPtr lev_dim = dom.add_dimension("level", 10);
  VarDescPtr t_desc = dom.add_var_desc("T", Float64Type, {col_dim, lev_dim}, "K");
  VarDescPtr q_desc = dom.add_var_desc("q", Float64Type, {col_dim, lev_dim}, "kg/kg");
  VarDescPtr p_surf_desc = dom.add_var_desc("p_surf", Float64Type, {col_dim}, "Pa");
  std::size_t t_idx = 0;
  std::size_t q_idx = t_desc->size();
  std::size_t p_idx = q_idx + q_desc->size();
  std::size_t total_size = p_idx + p_surf_desc->size();

  SECTION( "variables can be accessed through a variable array view" ) {
    double data[total_size];
    VariableArrayView<double> view({t_desc, q_desc, p_surf_desc}, &data[0]);
    REQUIRE( view.size() == total_size );
    auto t = view.get_variable("T");
    auto q = view.get_variable("q");
    auto p_surf = view.get_variable("p_surf");
    auto t_array = &data[t_idx];
    auto q_array = &data[q_idx];
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != t.size(); ++i) {
      t[i] = 2.*i;
    }
    for (int i = 0; i != q.size(); ++i) {
      q[i] = i + 0.5;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != t.size(); ++i) {
      REQUIRE( t[i] == 2.*i );
      REQUIRE( t_array[i] == 2.*i );
    }
    for (int i = 0; i != q.size(); ++i) {
      REQUIRE( q[i] == i + 0.5 );
      REQUIRE( q_array[i] == i + 0.5 );
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "getting a pointer to view data" ) {
    double data[total_size];
    VariableArrayView<double> view({t_desc, q_desc, p_surf_desc}, &data[0]);
    double *data_ptr = view.data();
    for (int i = 0; i != view.size(); ++i) {
      data_ptr[i] = 2.*i;
    }
    for (int i = 0; i != total_size; ++i) {
      REQUIRE( data[i] == 2.*i );
    }
  }

  SECTION( "getting a const variable from a const VariableArrayView" ) {
    double data[total_size];
    const double* construct_ptr = &data[0];
    const VariableArrayView<double> view = VariableArrayView<double>::make_const({t_desc, q_desc, p_surf_desc}, construct_ptr);
    REQUIRE( view.size() == total_size );
    auto t = view.get_variable("T");
    auto q = view.get_variable("q");
    auto p_surf = view.get_variable("p_surf");
    auto t_array = &data[t_idx];
    auto q_array = &data[q_idx];
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != t.size(); ++i) {
      t_array[i] = 2.*i;
    }
    for (int i = 0; i != q.size(); ++i) {
      q_array[i] = i + 0.5;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf_array[i] = -2. * i;
    }
    for (int i = 0; i != t.size(); ++i) {
      REQUIRE( t[i] == 2.*i );
      REQUIRE( t_array[i] == 2.*i );
    }
    for (int i = 0; i != q.size(); ++i) {
      REQUIRE( q[i] == i + 0.5 );
      REQUIRE( q_array[i] == i + 0.5 );
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "getting a const pointer to view data" ) {
    double data[total_size];
    const double* construct_ptr = &data[0];
    const VariableArrayView<double> view = VariableArrayView<double>::make_const({t_desc, q_desc, p_surf_desc}, construct_ptr);
    const double *data_ptr = view.data();
    for (int i = 0; i != total_size; ++i) {
      data[i] = 2.*i;
    }
    for (int i = 0; i != view.size(); ++i) {
      REQUIRE( data_ptr[i] == 2.*i );
    }
  }

  SECTION( "VariableArrayView can allocate and own its data" ) {
    VariableArrayView<double> view{t_desc, q_desc, p_surf_desc};
    REQUIRE( view.size() == t_desc->size() + q_desc->size() + p_surf_desc->size() );
    auto t = view.get_variable("T");
    auto q = view.get_variable("q");
    auto p_surf = view.get_variable("p_surf");
    for (int i = 0; i != t.size(); ++i) {
      t[i] = 2.*i;
    }
    for (int i = 0; i != q.size(); ++i) {
      q[i] = i + 0.5;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    // The following checks that the writes worked without overlapping.
    for (int i = 0; i != t.size(); ++i) {
      REQUIRE( t[i] == 2.*i );
    }
    for (int i = 0; i != q.size(); ++i) {
      REQUIRE( q[i] == i + 0.5 );
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
    }
  }

  SECTION( "owning views can be constructed with vectors" ) {
    std::vector<VarDescPtr> var_descs = {t_desc, q_desc, p_surf_desc};
    VariableArrayView<double> view(var_descs);
    REQUIRE( view.size() == t_desc->size() + q_desc->size() + p_surf_desc->size() );
    auto t = view.get_variable("T");
    for (int i = 0; i != t.size(); ++i) {
      t[i] = 2.*i;
    }
    for (int i = 0; i != t.size(); ++i) {
      REQUIRE( t[i] == 2.*i );
    }
  }

  SECTION( "data can be copied into a new owning view" ) {
    // Create "original" data.
    VariableArrayView<double> view{t_desc, q_desc, p_surf_desc};
    REQUIRE( view.size() == t_desc->size() + q_desc->size() + p_surf_desc->size() );
    auto t = view.get_variable("T");
    auto q = view.get_variable("q");
    auto p_surf = view.get_variable("p_surf");
    for (int i = 0; i != t.size(); ++i) {
      t[i] = 2.*i;
    }
    for (int i = 0; i != q.size(); ++i) {
      q[i] = i + 0.5;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    // Make new array.
    VariableArrayView <double> new_view = view.deep_copy();
    auto t_new = new_view.get_variable("T");
    auto q_new = new_view.get_variable("q");
    auto p_surf_new = new_view.get_variable("p_surf");
    // Clear original array.
    for (int i = 0; i != t.size(); ++i) {
      t[i] = 0.;
    }
    for (int i = 0; i != q.size(); ++i) {
      q[i] = 0.;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = 0.;
    }
    // The following checks that the data was copied.
    for (int i = 0; i != t_new.size(); ++i) {
      REQUIRE( t_new[i] == 2.*i );
    }
    for (int i = 0; i != q_new.size(); ++i) {
      REQUIRE( q_new[i] == i + 0.5 );
    }
    for (int i = 0; i != p_surf_new.size(); ++i) {
      REQUIRE( p_surf_new[i] == -2. * i );
    }
  }

  SECTION( "getting an error message when a variable is not in the array" ) {
    VariableArrayView<double> view{t_desc, q_desc, p_surf_desc};
    REQUIRE_THROWS_MATCHES( view.get_variable("invalid"), VariableNotFoundException,
                            MessageMatches(ContainsSubstring("variable array") && ContainsSubstring("invalid")) );
  }

}

TEST_CASE( "interacting with the model state", "[State]" ) {
  Domain dom;
  DimensionPtr col_dim = dom.add_dimension("column", 8);
  VarDescPtr p_surf_desc = dom.add_var_desc("p_surf", Float64Type, {col_dim}, "Pa");
  std::size_t p_idx = 0;
  std::size_t total_size = p_idx + p_surf_desc->size();

  SECTION( "construct and access data using the VariableArrayView interface" ) {
    double data[total_size];
    State<double> state({p_surf_desc}, &data[0]);
    REQUIRE( state.size() == total_size );
    auto p_surf = state.get_variable("p_surf");
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "construct const State with make_const" ) {
    double data[total_size];
    const State<double> state = State<double>::make_const({p_surf_desc}, &data[0]);
    REQUIRE( state.size() == total_size );
    auto p_surf = state.get_variable("p_surf");
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "construct owned State" ) {
    State<double> state{p_surf_desc};
    REQUIRE( state.size() == total_size );
    auto p_surf = state.get_variable("p_surf");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
    }
  }

  SECTION( "construct owned State from vector of VarDescPtr" ) {
    std::vector<VarDescPtr> var_descs = {p_surf_desc};
    State<double> state(var_descs);
    REQUIRE( state.size() == total_size );
    auto p_surf = state.get_variable("p_surf");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
    }
  }

  SECTION( "copy State into a new owning object" ) {
    State<double> state{p_surf_desc};
    auto p_surf = state.get_variable("p_surf");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    State<double> new_state = state.deep_copy();
    auto new_p_surf = new_state.get_variable("p_surf");
    // Clear original array to guarantee we aren't just pointing to it.
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = 0.;
    }
    REQUIRE( new_state.size() == total_size );
    for (int i = 0; i != new_p_surf.size(); ++i) {
      REQUIRE( new_p_surf[i] == -2. * i );
    }
  }
}

TEST_CASE( "interacting a tendency on the model state", "[Tendency]" ) {
  Domain dom;
  DimensionPtr col_dim = dom.add_dimension("column", 8);
  VarDescPtr p_surf_desc = dom.add_var_desc("p_surf_tend", Float64Type, {col_dim}, "Pa");
  std::size_t p_idx = 0;
  std::size_t total_size = p_idx + p_surf_desc->size();

  SECTION( "construct and access data using the VariableArrayView interface" ) {
    double data[total_size];
    Tendency<double> tend({p_surf_desc}, &data[0]);
    REQUIRE( tend.size() == total_size );
    auto p_surf = tend.get_variable("p_surf_tend");
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "construct const Tendency with make_const" ) {
    double data[total_size];
    const Tendency<double> tend = Tendency<double>::make_const({p_surf_desc}, &data[0]);
    REQUIRE( tend.size() == total_size );
    auto p_surf = tend.get_variable("p_surf_tend");
    auto p_surf_array = &data[p_idx];
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
      REQUIRE( p_surf_array[i] == -2. * i );
    }
  }

  SECTION( "construct owned Tendency" ) {
    Tendency<double> tend{p_surf_desc};
    REQUIRE( tend.size() == total_size );
    auto p_surf = tend.get_variable("p_surf_tend");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
    }
  }

  SECTION( "construct owned Tendency from vector of VarDescPtr" ) {
    std::vector<VarDescPtr> var_descs = {p_surf_desc};
    Tendency<double> tend(var_descs);
    REQUIRE( tend.size() == total_size );
    auto p_surf = tend.get_variable("p_surf_tend");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    for (int i = 0; i != p_surf.size(); ++i) {
      REQUIRE( p_surf[i] == -2. * i );
    }
  }

  SECTION( "copy Tendency into a new owning object" ) {
    Tendency<double> tend{p_surf_desc};
    auto p_surf = tend.get_variable("p_surf_tend");
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = -2. * i;
    }
    Tendency<double> new_tend = tend.deep_copy();
    auto new_p_surf = new_tend.get_variable("p_surf_tend");
    // Clear original array to guarantee we aren't just pointing to it.
    for (int i = 0; i != p_surf.size(); ++i) {
      p_surf[i] = 0.;
    }
    REQUIRE( new_tend.size() == total_size );
    for (int i = 0; i != new_p_surf.size(); ++i) {
      REQUIRE( new_p_surf[i] == -2. * i );
    }
  }
}

// Utility macro for testing the copy constructors of exceptions.
// These constructors are pretty trivial pieces of code, but it's easy
// to forget to implement them if they are not used in the library itself.
#define check_exception_copy( T, ... ) { T original(__VA_ARGS__); T copy (original); REQUIRE( original.what() == copy.what() ); }

TEST_CASE( "Copy constructing exceptions" , "[exceptions]" ) {

  // Note that comparing the "what()" outputs here is doing a comparison of the
  // returned char pointers, which is what we want to test because the copy
  // constructor for exceptions should point to the same string.

  SECTION( "DimensionNotFoundException can be copied" ) {
    check_exception_copy(DimensionNotFoundException, "foo", "bar");
  }

  SECTION( "TypeMismatchException can be copied" ) {
    check_exception_copy(TypeMismatchException, "foo1", "foo2", "bar");
  }

  SECTION( "UnreachableException can be copied" ) {
    check_exception_copy(UnreachableException, "foobar");
  }

  SECTION( "VariableNotFoundException can be copied" ) {
    check_exception_copy(VariableNotFoundException, "foo", "bar");
  }

}

#undef check_exception_copy
