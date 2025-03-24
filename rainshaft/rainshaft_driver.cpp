#include "spaecies.hpp"
#include "evaporation.hpp"
#include "explicit_integrator.hpp"
#include "fixed_substep_integrator.hpp"
#include "forcing_integrator.hpp"
#include "limiting_integrator.hpp"
#include "operator_splitting_integrator.hpp"
#include "nudging.hpp"
#include "parsing_utilities.hpp"
#include "rainshaft_constants.hpp"
#include "rainshaft_grid.hpp"
#include "rainshaft_initial.hpp"
#include "rainshaft_solution.hpp"
#include "rainshaft_tend_descs.hpp"
#include "rainshaft_ncio.hpp"
#include "rainshaft_sum_process.hpp"
#include "saturation.hpp"
#include "sed_cfl_integrator.hpp"
#include "sedimentation.hpp"
#include "self_collision.hpp"
#include "sequential_split_integrator.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <chrono>
#include <memory>
#include <optional>
#include <string>
#include "imex_integrator.hpp"
#include "mri_integrator.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;

  // Arg parser
  double dt, dt_partition_1, dt_partition_2, final_time, rel_tol;
  int steps_per_output, num_cases;
  std::string input_file, output_file, method_type, initial_condition, initial_condition_file;
  std::size_t order, icase_in;
  bool do_nudging, postprocess, use_lookup;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
    ("i", po::value<std::string>(&input_file), "(optional) input file for command line arguments")
		("order", po::value<std::size_t>(&order)->default_value(2), "order of method")
		("dt", po::value<double>(&dt)->default_value(0.1), "step size. if integration type is set to MRI/splitting/forcing, this argument is the outer time step size")
    ("dt_partition_1", po::value<double>(&dt_partition_1)->default_value(0), "step size of partition 1 for MRI/splitting/forcing methods. negative values indicate ratios, e.g. dt_partition_1=-0.5 corresponds to dt_partition_1 = dt/2")
    ("dt_partition_2", po::value<double>(&dt_partition_2)->default_value(0), "step size of partition 2 for MRI/splitting/forcing methods. negative values indicate ratios, e.g. dt_partition_2=-0.5 corresponds to dt_partition_2 = dt/2")
    ("rel_tol", po::value<double>(&rel_tol)->default_value(1.e-4), "relative tolerance for adaptive stepping and nonlinear solvers")
    ("postprocess", po::value<bool>(&postprocess)->default_value(false), "postprocesses stages and steps to be positive")
    ("use_lookup", po::value<bool>(&use_lookup)->default_value(false), "uses lookup tables to evaluate fall speeds")
		("type", po::value<std::string>(&method_type)->default_value("explicit"), "type of integrator (e.g. explicit, implicit, imex, mri, original)")
    ("steps", po::value<int>(&steps_per_output)->default_value(-1), "frequency of saved output, e.g. write output to netCDF every steps_per_output timesteps. -1 to save only the first and last steps")
    ("ic_file", po::value<std::string>(&initial_condition_file), "type of initial condiiton (e.g. 'adiabatic' or the filename of E3SM data)")
    ("num_cases", po::value<int>(&num_cases)->default_value(1), "number of E3SM cases to load if initial_condition is set to filename. -1 for all cases.")
    ("case_idx", po::value<std::size_t>(&icase_in), "(optional) specific E3SM case to load if initial_condition is set to filename.")
    ("final_time", po::value<double>(&final_time)->default_value(1800.0), "stopping time for integration")
    ("nudging", po::value<bool>(&do_nudging)->default_value(true), "boolean flag for nudging")
    ("filename", po::value<std::string>(&output_file)->default_value("rainshaft.nc"), "savefile name")
  ;

  // Load from command line to check for input file
  po::variables_map vm{};
  po::store(parse_command_line(argc, argv, desc), vm);
  po::notify(vm);  

  // Load from settings file if available
  if (vm.count("i")) {
    std::cout << "Setting parameters from input file " << vm["i"].as<std::string>() << "..." << std::endl;
    parse_input_file(vm["i"].as<std::string>());

    std::ifstream settings_file("settings_filtered.ini");
    vm = po::variables_map();
    po::store(po::parse_config_file(settings_file, desc), vm);
    po::notify(vm); 
    std::filesystem::remove("settings_filtered.ini");
  } else {
    std::cout << "Setting parameters from command line arguments..." << std::endl;
  }
  std::cout << "---------------------------------------------------" << std::endl;
  
  // Setup dependencies
  option_dependency<std::string>(vm, "case_idx", "ic_file");    // specifying a particular case_idx requires input E3SM data file

  // Print input to program_options
  print_variables_map(vm);

  // Help message
	if (vm.count("help"))
	{
		std::cout << desc << "\n";
		return 1;
	}

  // If no IC file specified, default to adiabatic IC
  if (vm.count("ic_file")) {
    initial_condition = initial_condition_file;
  } else {
    initial_condition = "adiabatic";
    std::cout << "No IC input file specified. Defaulting to adiabatic case." << std::endl;
  }

  // negative partition time steps corresponds to scalings of dt
  if (dt_partition_1 < 0.0) {
    dt_partition_1 = abs(dt_partition_1) * dt;
  }
  if (dt_partition_2 < 0.0) {
    dt_partition_2 = abs(dt_partition_2) * dt;
  }
 
  // Set up model constants.
  // SPS: Choose rho_top in a more principled way?
  RainshaftConstants constants{3.14159265358979323846,
                               287.04, 1.00464e3, 461.50, 997., 2.501e6,
                               0.62197, 1.e-14, 9.80616, 1.e-5, 5.e-3,
                               0.988919555598356, 1.e3, 1.e-4};
  // Approximate model top in meters.
  // (The grid maker will actually use the next higher-altitude E3SM level.)
  double model_top = 2.e3;
  // Surface pressure and temperature in Pa and K, respectively.
  double srf_pres = 1.e5, srf_temp = 293.15;
  // Lapse rate of initial condition in K/m.
  double lapse_rate = 6.5e-3;
  // Relative humidity of initial condition.
  double rel_hum_init = 0.7;
  // Time scale over which to nudge t and q back to initial condition in seconds.
  double nudge_time_scale = 15. * 60.;
  // Time of simulation start.
  double initial_time = 0.;
  // Saturation vapor pressure formula.
  SaturationFormulae sat_form(constants);

  // Set up max grid size
  std::size_t max_levs;
  int max_cases;
  RainshaftGrid default_grid = make_e3sm_like_grid(constants, model_top, srf_pres,
                                                   srf_temp, lapse_rate);
  std::optional<NetcdfReader> reader;
  if (initial_condition == "adiabatic") {
    max_cases = 1;
    num_cases = 1;
    max_levs = default_grid.nlev;
  } else {
    reader = NetcdfReader(initial_condition);
    std::tuple<std::size_t, std::size_t> cases_levs = reader->read_num_cases_and_max_levs();
    max_cases = std::get<0>(cases_levs);

    if (vm.count("num_cases")) {
      if (max_cases < num_cases) {
        throw std::logic_error("Requested more cases than exist on file.");
      }

      if (num_cases == -1) {
        num_cases = max_cases;
      }
    }

    if (vm.count("case_idx")) {
      if (vm.count("num_cases") && (vm["num_cases"].as<int>() != 1)) {
        num_cases = 1;
        std::cout << "Requested a specific case_idx but also specified num_cases != 1. Defaulting to num_cases=1." << std::endl;
      }

      if (max_cases < icase_in) {
        throw std::logic_error("Requested non-existent case_idx in file.");
      }
    }
    // Remove cloud base from levels.
    max_levs = std::get<1>(cases_levs) - 1;
  }
  NetcdfWriter writer(output_file, num_cases, max_levs);

  // Setup list of cases to run
  std::vector<std::size_t> cases_to_run(num_cases);
  if (vm.count("case_idx")) {
    cases_to_run.at(0) = icase_in;
  } else {
    std::iota(cases_to_run.begin(), cases_to_run.end(), 0);
  }

  // Physics types that won't vary between cases (declare here to avoid calculating the
  // lookup tables in the case loop).
  // Sedimentation process.
  Sedimentation sed(constants, use_lookup, false);
  // Self-collision processes.
  SelfCollision self_coll;
  // Evaporation process.
  std::optional<double> dt_for_evap = std::nullopt;
  if (method_type == "original") {
    dt_for_evap = dt;
  }
  Evaporation evap(constants, sat_form, use_lookup, false, dt_for_evap);

  // for (std::size_t icase = 0; icase != num_cases; ++icase) {
  for (const std::size_t &icase : cases_to_run) {
    // Grid for this case.
    RainshaftGrid grid = (initial_condition == "adiabatic")
      ? default_grid : reader->read_grid(icase);
    std::size_t nlev = grid.nlev;
    // Set boundary condition if needed.
    if (initial_condition != "adiabatic") {
      reader->read_boundary_conditions(icase, constants);
    }

    spaecies::Domain dom;
    spaecies::DimensionPtr lev_dim = dom.add_dimension("level", nlev);
    spaecies::VarDescPtr t_desc = dom.add_var_desc("T", spaecies::Float64Type, {lev_dim}, "K");
    spaecies::VarDescPtr q_desc = dom.add_var_desc("q", spaecies::Float64Type, {lev_dim}, "kg/kg");
    spaecies::VarDescPtr nr_desc = dom.add_var_desc("nr", spaecies::Float64Type, {lev_dim}, "1/kg");
    spaecies::VarDescPtr qr_desc = dom.add_var_desc("qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
    VarDescList state_descs = {t_desc, q_desc, nr_desc, qr_desc};
    VarDescList tend_descs = tend_descs_from_state_descs(dom, state_descs);
    State initial_state(state_descs);

    // Set up initial condition.
    if (initial_condition == "adiabatic") {
      bool converged = warm_adiabatic_initial_condition(constants, grid, sat_form, srf_temp,
                                                        lapse_rate, rel_hum_init, initial_state);
      if (!converged) {
        std::cerr << "Initial condition iteration failed to converge." << std::endl;
        return 1;
      }
    } else {
      reader->read_initial_conditions(icase, initial_state);
    }
    RainshaftDerivedVars initial_dvars = RainshaftDerivedVars(constants, grid, initial_state);

    // Nudging to initial condition.
    std::shared_ptr<Nudging> nudge;
    // Sum of all processes.
    std::vector<const RainshaftProcess *> partition_2_process_vec{}, all_process_vec{};
    if (do_nudging) {
      // SPS: Need some kind of span-like interface to avoid having to
      // do this copy.
      VarMut t = initial_state.get_variable("T");
      VarMut q = initial_state.get_variable("q");
      std::vector<double> t_vec, q_vec;
      for (std::size_t i = 0; i != nlev; ++i) {
        t_vec.push_back(t[i]);
        q_vec.push_back(q[i]);
      }
      nudge = std::make_shared<Nudging>(nudge_time_scale, t_vec, q_vec);
      partition_2_process_vec = {&evap, &*nudge, &self_coll};
      all_process_vec = {&sed, &*nudge, &self_coll, &evap};
    } else {
      partition_2_process_vec = {&evap, &self_coll};
      all_process_vec = {&sed, &self_coll, &evap};
    }
    SumProcess partition_2_processes(partition_2_process_vec);
    // Sum of local processes.
    SumProcess partition_1_processes = SumProcess{{&sed}};

    SumProcess all_processes = SumProcess{all_process_vec};

    // List of integrators that need to remain allocated for use by original scheme.
    std::vector<std::shared_ptr<RainshaftIntegrator>> backing_integrators;

    const auto intg = [&]() -> std::unique_ptr<RainshaftIntegrator> {
      if (method_type == "explicit") {
        return std::make_unique<ExplicitIntegrator>(constants, grid, &all_processes, state_descs, tend_descs, dt, order, rel_tol, postprocess, steps_per_output);
      } else if (method_type == "implicit") {
        return std::make_unique<IMEXIntegrator>(constants, grid, &all_processes, nullptr, state_descs, tend_descs, dt, order, rel_tol, postprocess, steps_per_output);  
      } else if (method_type == "imex") {
        return std::make_unique<IMEXIntegrator>(constants, grid, &partition_1_processes, &partition_2_processes, state_descs, tend_descs, dt, order, rel_tol, postprocess, steps_per_output);
      } else if (method_type == "mri") {
        return std::make_unique<MRIIntegrator>(constants, grid, &partition_1_processes, &partition_2_processes, nullptr, state_descs, tend_descs, dt_partition_1, dt, order, rel_tol, postprocess, steps_per_output);
      } else if (method_type == "splitting") {
        return std::make_unique<OperatorSplittingIntegrator>(constants, grid, &partition_1_processes, &partition_2_processes, state_descs, tend_descs, dt, dt_partition_1, dt_partition_2, order, rel_tol, postprocess, steps_per_output);
      } else if (method_type == "forcing") {
        return std::make_unique<ForcingIntegrator>(constants, grid, &partition_1_processes, &partition_2_processes, state_descs, tend_descs, dt, dt_partition_1, dt_partition_2, postprocess, steps_per_output);
      } else if (method_type == "original") {
        std::shared_ptr<ExplicitIntegrator> local_intg = std::make_shared<ExplicitIntegrator>(constants, grid, &partition_2_processes, state_descs, tend_descs, dt, 1, rel_tol);
        backing_integrators.emplace_back(local_intg);
        std::shared_ptr<LimitingIntegrator> local_lim_intg = std::make_shared<LimitingIntegrator>(constants, *local_intg);
        backing_integrators.emplace_back(local_lim_intg);
        std::shared_ptr<SedCflIntegrator> sed_intg = std::make_shared<SedCflIntegrator>(&constants, &grid, tend_descs, &sed);
        backing_integrators.emplace_back(sed_intg);
        std::shared_ptr<SequentialSplitIntegrator> step_intg = std::make_shared<SequentialSplitIntegrator>(std::vector<const RainshaftIntegrator*>({&*local_lim_intg, &*sed_intg}));
        backing_integrators.emplace_back(step_intg);
        return std::make_unique<FixedSubstepIntegrator>(&*step_intg, dt);
      } else {
        std::cout << method_type << std::endl;
        throw std::logic_error("Invalid time integration method type");
      }
    }();

    auto before_sol = high_resolution_clock::now();
    RainshaftSolution solution = intg->integrate(initial_time, final_time, initial_state);
    auto after_sol = high_resolution_clock::now();

    // Time taken for solution.
    duration<double, std::milli> walltime_ms = after_sol - before_sol;
    std::cout << "Case " << icase << ", Time: " << walltime_ms.count() << std::endl;
    // Write out grid and all states.
    std::size_t icase_writer;
    if (vm.count("case_idx")) {
      icase_writer = 0;
    } else {
      icase_writer = icase;
    }
    writer.write_grid(grid, icase_writer);
    writer.write_states(solution.states, icase_writer);
    std::vector<RainshaftDerivedVars> solution_dvars;
    solution_dvars.reserve(solution.states.size());
    for (StateConst state : solution.states) {
      solution_dvars.emplace_back(constants, grid, state);
    }
    writer.write_derived_vars(solution_dvars, icase_writer);
    writer.write_num_rhs_evals(solution.num_rhs_evals, icase_writer);
    writer.write_walltime_ms(walltime_ms.count(), icase_writer);

    // write settings as metadata
    writer.write_metadata(order, dt, dt_partition_1, dt_partition_2, rel_tol, postprocess,
      use_lookup, method_type, steps_per_output, initial_condition_file,
      num_cases, icase_in, final_time, do_nudging);
  }
  return 0;
}
