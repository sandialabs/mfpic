#include "libmfpic/MeshUtilities.hpp"
#include <libmfpic/BuildOutputParametersFromYaml.hpp>
#include <libmfpic/BuildParticleBoundariesFromYaml.hpp>
#include <libmfpic/BuildParticlesFromYaml.hpp>
#include <libmfpic/BuildSpeciesMapFromYaml.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/RunSimulation.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/TimeSteppingFactory.hpp>
#include <libmfpic/VerletTimeIntegrator.hpp>

#include <mfem/mfem.hpp>

#include <yaml-cpp/yaml.h>


namespace mfpic {

void runSimulation(int argc, char* argv[]) {
  std::string input_deck_filename = "mfpic.yaml";

  mfem::OptionsParser options_parser(argc, argv);
  options_parser.AddOption(&input_deck_filename, "-i", "--input-deck", "Input deck to read from.");
  options_parser.ParseCheck(std::cout);

  const YAML::Node main = YAML::LoadFile(input_deck_filename);

  const MeshParameters mesh_parameters = buildMeshParametersFromYAML(main["Mesh"]);
  auto mesh = std::make_shared<mfem::Mesh>(buildMesh(mesh_parameters));

  const YAML::Node fields = main["Fields"];
  const int electrostatic_basis_order = fields["Basis Order"].as<int>();
  Discretization electrostatic_discretization(mesh.get(), electrostatic_basis_order);

  const int mesh_dimension = mesh->Dimension();
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    fields["Boundary Conditions"], mesh_dimension);
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs = buildDirichletBoundaryConditions(
    boundary_attribute_to_dirichlet_value, electrostatic_discretization);

  ElectrostaticFieldOperations electrostatic_field_operations(electrostatic_discretization, std::move(dirichlet_bcs));
  ElectrostaticFieldState electrostatic_field_state(electrostatic_discretization);

  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(main["Species"]);

  const auto [particle_boundary_factories, default_particle_boundary_factory] = buildParticleBoundariesFromYaml(
    main["Particles"],
    mesh_dimension
  );
  ParticleOperations particle_operations(
    electrostatic_discretization,
    particle_boundary_factories,
    default_particle_boundary_factory
  );

  std::default_random_engine generator;
  ParticleContainer particle_container = buildParticlesFromYaml(
    main["Particles"]["Initial Conditions"],
    species_map,
    generator,
    mesh
  );
  dumpParticles(particle_container, 0.0);

  std::vector<LowFidelityState> low_fidelity_states;
  std::vector<std::unique_ptr<LowFidelityOperations>> low_fidelity_operations;

  int dg_euler_order = 0;
  if (main["Euler Fluids"]) {
    const YAML::Node euler_fluids = main["Euler Fluids"];
    if (euler_fluids["Basis Order"]) {
      dg_euler_order = euler_fluids["Basis Order"].as<int>();
    }
  }
  Discretization dg_euler_discretization(mesh.get(), dg_euler_order, FETypes::DG, euler::ConservativeVariables::NUM_VARS);

  if (main["Euler Fluids"]) {
    const YAML::Node euler_fluids = main["Euler Fluids"];
    std::vector<std::unique_ptr<SourceParameters>> list_of_parameters = buildListOfSourceParametersFromYAML(
      euler_fluids["Initial Conditions"], species_map);
    LowFidelityState dg_euler_state = buildEulerState(dg_euler_discretization, list_of_parameters);
    low_fidelity_states.push_back(dg_euler_state);

    std::vector<Species> species_list = dg_euler_state.getSpeciesList();
    std::vector<std::unique_ptr<DGGhostBC>> dg_euler_bcs;
    std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
      dg_euler_discretization,
      electrostatic_discretization,
      species_list,
      dg_euler_bcs);
    low_fidelity_operations.push_back(std::move(dg_euler_operations));
  }

  OutputParameters output_parameters;
  if (main["Output"].IsDefined())
    output_parameters = buildOutputParametersFromYAML(main["Output"]);

  MeshDataWriter mesh_data_writer(output_parameters.mesh_output_folder_name, *mesh);

  IntegratedCharge integrated_charge = particle_operations.assembleCharge(particle_container);

  electrostatic_field_operations.fieldSolve(electrostatic_field_state, integrated_charge);
  mesh_data_writer.output(electrostatic_field_state, low_fidelity_states, 0, 0.);

  std::ofstream csv_file("output.csv");
  csv_file << std::setprecision(std::numeric_limits<double>::digits);
  csv_file << "# Time_Step Time Field_Energy" << std::endl;
  csv_file << 0 << " " << 0.0 << " " << electrostatic_field_operations.fieldEnergy(electrostatic_field_state) << std::endl;

  TimeSteppingParameters time_stepping_parameters = buildTimeSteppingParametersFromYAML(main["Time Stepping"]);
  VerletTimeIntegrator verlet_time_integrator(electrostatic_discretization);
  const double smallest_cell_lengthscale = getSmallestCellLengthscale(*mesh);
  for (int i_timestep = 1; i_timestep <= time_stepping_parameters.number_of_timesteps; ++i_timestep) {
    const double timestep_size = time_stepping_parameters.timestep_size;
    const double begin_time = (i_timestep - 1) * timestep_size;
    const double end_time = i_timestep * timestep_size;

    std::cout << "Time Step: " << i_timestep << "    Time: " << begin_time << std::endl;

    verlet_time_integrator.advanceTimestep(
      low_fidelity_states,
      low_fidelity_operations,
      particle_container,
      particle_operations,
      electrostatic_field_state,
      electrostatic_field_operations,
      timestep_size
    );

    double cfl = 0.;

    for (size_t i_lf_model = 0; i_lf_model < low_fidelity_operations.size(); ++i_lf_model)
      cfl = fmax(cfl, low_fidelity_operations[i_lf_model]->estimateCFL(timestep_size, smallest_cell_lengthscale));

    std::cout << "    Maximum CFL: " << cfl << std::endl;

    if (main["Particles"]["Sources"].IsDefined()) {
      particle_container.addParticles(buildParticlesFromYaml(
        main["Particles"]["Sources"],
        species_map,
        generator,
        mesh
      ));
    }

    if (i_timestep % output_parameters.output_stride == 0) {
      dumpParticles(particle_container, end_time, output_parameters.particle_dump_filename);
      mesh_data_writer.output(electrostatic_field_state, low_fidelity_states, i_timestep, end_time);
      csv_file << i_timestep << " " << end_time << " " << electrostatic_field_operations.fieldEnergy(electrostatic_field_state) << std::endl;
    }
  }

}

}
