#include <libmfpic/BuildElectrostaticFieldOperationsFromYaml.hpp>
#include <libmfpic/BuildOutputParametersFromYaml.hpp>
#include <libmfpic/BuildParticleBoundariesFromYaml.hpp>
#include <libmfpic/BuildParticlesFromYaml.hpp>
#include <libmfpic/BuildSpeciesMapFromYaml.hpp>
#include <libmfpic/DGGhostBC.hpp>
#include <libmfpic/DGEulerBoundaryConditionsFactory.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/DGEulerOperationsFactory.hpp>
#include <libmfpic/DirichletBoundaryConditions.hpp>
#include <libmfpic/DirichletBoundaryConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/ElectrostaticFieldOperations.hpp>
#include <libmfpic/ElectrostaticFieldState.hpp>
#include <libmfpic/Euler.hpp>
#include <libmfpic/LowFidelityOperations.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/MeshDataWriter.hpp>
#include <libmfpic/MeshFactory.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/RunSimulation.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/TextDataWriter.hpp>
#include <libmfpic/TimeIntegrator.hpp>
#include <libmfpic/TimeSteppingFactory.hpp>

#include <memory>
#include <mfem/mfem.hpp>

#include <yaml-cpp/yaml.h>


namespace mfpic {

void runSimulation(int argc, char* argv[]) {
  std::string input_deck_filename = "mfpic.yaml";

  mfem::OptionsParser options_parser(argc, argv);
  options_parser.AddOption(&input_deck_filename, "-i", "--input-deck", "Input deck to read from.");
  options_parser.ParseCheck(std::cout);

  YAML::Node main = YAML::LoadFile(input_deck_filename);

  const MeshParameters mesh_parameters = buildMeshParametersFromYAML(main["Mesh"]);
  auto mesh = std::make_shared<mfem::Mesh>(buildMesh(mesh_parameters));

  YAML::Node fields = main["Fields"];
  int electrostatic_basis_order = 1;
  if (fields["Basis Order"]) {
    electrostatic_basis_order = fields["Basis Order"].as<int>();
  }
  Discretization electrostatic_discretization(mesh.get(), electrostatic_basis_order);

  const int mesh_dimension = mesh->Dimension();
  std::unordered_map<int, double> boundary_attribute_to_dirichlet_value = buildBoundaryAttributeToDirichletValueFromYAML(
    fields["Boundary Conditions"], mesh_dimension);
  std::unique_ptr<DirichletBoundaryConditions> dirichlet_bcs = buildDirichletBoundaryConditions(
    boundary_attribute_to_dirichlet_value, electrostatic_discretization);

  std::unique_ptr<ElectrostaticFieldOperations> electrostatic_field_operations = buildElectrostaticFieldOperationsFromYaml(
    fields,
    electrostatic_discretization,
    std::move(dirichlet_bcs)
  );
  ElectrostaticFieldState particle_electrostatic_field_state(electrostatic_discretization);
  ElectrostaticFieldState variance_reduced_particle_electrostatic_field_state(electrostatic_discretization);

  std::unordered_map<std::string, Species> species_map = buildSpeciesMapFromYaml(main["Species"]);

  const auto [particle_boundary_factories, default_particle_boundary_factory] = buildParticleBoundariesFromYaml(
    main["Particles"],
    mesh_dimension
  );
  ParticleOperations particle_operations(
    electrostatic_discretization,
    particle_boundary_factories,
    default_particle_boundary_factory,
    species_map
  );

  //std::default_random_engine generator;
  std::random_device rd;
  std::default_random_engine generator(rd());
  ParticleContainer particle_container = buildParticlesFromYaml(
    main["Particles"]["Initial Conditions"],
    species_map,
    generator,
    mesh
  );
  const std::string prefix = "particle_moments";
  dumpParticleMoments(particle_operations,particle_container, prefix, 0, 0.0);
  dumpParticles(particle_container, 0.0);

  std::vector<LowFidelityState> low_fidelity_states;
  std::vector<std::unique_ptr<LowFidelityOperations>> low_fidelity_operations;
  std::vector<ElectrostaticFieldState> low_fidelity_field_states;

  YAML::Node euler_fluids = main["Euler Fluids"];
  int dg_euler_order = 0;
  bool push_low_fidelity_with_particle_fields = false;
  if (euler_fluids["Basis Order"]) {
    dg_euler_order = euler_fluids["Basis Order"].as<int>();
  }
  if (euler_fluids["Use Particle Fields"]) {
    push_low_fidelity_with_particle_fields = euler_fluids["Use Particle Fields"].as<bool>();
  }
  Discretization dg_euler_discretization(mesh.get(), dg_euler_order, FETypes::DG, euler::ConservativeVariables::NUM_VARS);

  std::vector<std::unique_ptr<SourceParameters>> list_of_euler_ic_parameters = buildListOfSourceParametersFromYAML(
    euler_fluids["Initial Conditions"], species_map);
  std::vector<std::unique_ptr<SourceParameters>> list_of_euler_source_parameters;
  if (euler_fluids["Sources"]) {
    list_of_euler_source_parameters = buildListOfSourceParametersFromYAML(euler_fluids["Sources"], species_map);
  }

  if (not list_of_euler_ic_parameters.empty()) {
    LowFidelityState dg_euler_state = buildEulerState(dg_euler_discretization, list_of_euler_ic_parameters);
    low_fidelity_states.push_back(dg_euler_state);

    std::vector<Species> species_list = dg_euler_state.getSpeciesList();
    std::unordered_map<int, DGEulerBCType> boundary_attribute_to_bc_type = buildBoundaryAttributeToBCTypeFromYAML(
      euler_fluids["Boundary Conditions"], mesh_dimension);
    std::vector<std::vector<std::unique_ptr<DGBC>>> dg_euler_bcs = buildDGEulerBoundaryConditions(boundary_attribute_to_bc_type, *mesh, species_list);
    std::vector<std::pair<Species, std::unique_ptr<mfem::VectorCoefficient>>> dg_euler_sources = buildListOfSpeciesAndEulerSourceCoefficients(list_of_euler_source_parameters);
    std::unique_ptr<LowFidelityOperations> dg_euler_operations = buildDGEulerOperations(
      dg_euler_discretization,
      electrostatic_discretization,
      species_list,
      dg_euler_bcs,
      dg_euler_sources);
    low_fidelity_operations.push_back(std::move(dg_euler_operations));
    low_fidelity_field_states.emplace_back(electrostatic_discretization);
  }

  for (int i = 0; i < std::ssize(low_fidelity_field_states); ++i) {
    auto* ops = dynamic_cast<DGEulerOperations*>(low_fidelity_operations[i].get());
    if (!ops) throw std::runtime_error("Cannot compute variance reduced moments for low_fidelity_operations[i] that is not DGEulerOperations.");
    const std::string variance_reduced_prefix = "variance_reduced_particle_moments";
    dumpVarianceReducedParticleMoments(particle_operations,particle_container,low_fidelity_states[i],*ops,variance_reduced_prefix, 0, 0.0);
    const std::string low_fidelity_prefix = "low_fidelity_moments";
    dumpLowFidelityMoments(low_fidelity_states[i],*ops,low_fidelity_prefix, 0, 0.0);
  }

  OutputParameters output_parameters;
  if (main["Output"].IsDefined())
    output_parameters = buildOutputParametersFromYAML(main["Output"]);

  MeshDataWriter mesh_data_writer(output_parameters.mesh_output_folder_name, *mesh);
  MeshDataWriter variance_reduced_mesh_data_writer(output_parameters.mesh_output_folder_name + "_VR", *mesh);

  {
    IntegratedCharge integrated_charge = particle_operations.assembleCharge(particle_container);
    electrostatic_field_operations->fieldSolve(particle_electrostatic_field_state, integrated_charge);
  }

  {
    //assume only one low fidelity state
    IntegratedCharge variance_reduced_integrated_charge = particle_operations.assembleVarianceReducedCharge(particle_container,low_fidelity_states[0],*low_fidelity_operations[0]);
    electrostatic_field_operations->fieldSolve(variance_reduced_particle_electrostatic_field_state, variance_reduced_integrated_charge);
  }
  
  for (int i = 0; i < std::ssize(low_fidelity_field_states); ++i) {
    IntegratedCharge integrated_charge = low_fidelity_operations[i]->assembleCharge(low_fidelity_states[i]);
    electrostatic_field_operations->fieldSolve(low_fidelity_field_states[i], integrated_charge);
  }


  mesh_data_writer.output(particle_electrostatic_field_state, low_fidelity_field_states, low_fidelity_states, 0, 0.);
  variance_reduced_mesh_data_writer.output(variance_reduced_particle_electrostatic_field_state, low_fidelity_field_states, low_fidelity_states, 0, 0.);

  const int num_low_fidelity_models = std::ssize(low_fidelity_states);
  TextDataWriter text_data_writer(num_low_fidelity_models,"output");
  text_data_writer.output(
    particle_electrostatic_field_state,
    low_fidelity_field_states,
    *electrostatic_field_operations,
    low_fidelity_states,
    low_fidelity_operations,
    0,
    0.);

  TextDataWriter variance_reduced_text_data_writer(num_low_fidelity_models,"vr_output");
  variance_reduced_text_data_writer.output(
    variance_reduced_particle_electrostatic_field_state,
    low_fidelity_field_states,
    *electrostatic_field_operations,
    low_fidelity_states,
    low_fidelity_operations,
    0,
    0.);

  TimeSteppingParameters time_stepping_parameters = buildTimeSteppingParametersFromYAML(main["Time Stepping"]);
  std::unique_ptr<TimeIntegrator> time_integrator = buildTimeIntegrator(
    time_stepping_parameters,
    electrostatic_discretization,
    push_low_fidelity_with_particle_fields);
  const double smallest_cell_lengthscale = getSmallestCellLengthscale(*mesh);
  for (int i_timestep = 1; i_timestep <= time_stepping_parameters.number_of_timesteps; ++i_timestep) {
    const double timestep_size = time_stepping_parameters.timestep_size;
    const double begin_time = (i_timestep - 1) * timestep_size;
    const double end_time = i_timestep * timestep_size;

    std::cout << "Time Step: " << i_timestep << "    Time: " << begin_time << std::endl;

    time_integrator->advanceTimestep(
      low_fidelity_states,
      low_fidelity_field_states,
      low_fidelity_operations,
      particle_container,
      particle_operations,
      particle_electrostatic_field_state,
      *electrostatic_field_operations,
      timestep_size
    );

    double cfl = 0.;

    for (size_t i_lf_model = 0; i_lf_model < low_fidelity_operations.size(); ++i_lf_model) {
      cfl = fmax(cfl, low_fidelity_operations[i_lf_model]->estimateCFL(timestep_size, smallest_cell_lengthscale));
    }

    std::cout << "    Maximum CFL: " << cfl << std::endl;

    if (main["Particles"]["Sources"].IsDefined()) {
      ParticleContainer source_particles = buildParticlesFromYaml(
      main["Particles"]["Sources"],
      species_map,
      generator,
      mesh
      );

      //F updated based on constant values
      // double source_number_density = 6e19; 
      // double source_temperature = 116045.18121550081; 
      // mfem::Vector source_bulk_velocity(3);
      // source_bulk_velocity = 0.0;
      // particle_operations.updateParticleDistributionFunctionValue(source_particles,source_number_density, source_bulk_velocity, source_temperature);

      //F Updated based on joint particle moments
      //ParticleContainer joint_particles; 
      // joint_particles.addParticles(source_particles);
      // joint_particles.addParticles(particle_container);
      // particle_operations.updateParticleDistributionFunctionValue(source_particles,joint_particles);

      //F Updated based on existing particle moments
      //particle_operations.updateParticleDistributionFunctionValue(source_particles,particle_container);

      particle_container.addParticles(source_particles);
    }

    if (i_timestep % output_parameters.output_stride == 0) {
      const std::string prefix = "particle_moments";
      dumpParticleMoments(particle_operations,particle_container, prefix, i_timestep, end_time);
      dumpParticles(particle_container, end_time, output_parameters.particle_dump_filename);
      for (int i = 0; i < std::ssize(low_fidelity_field_states); ++i) {
        auto* ops = dynamic_cast<DGEulerOperations*>(low_fidelity_operations[i].get());
        if (!ops) throw std::runtime_error("Cannot compute variance reduced moments for low_fidelity_operations[i] that is not DGEulerOperations.");
        const std::string variance_reduced_prefix = "variance_reduced_particle_moments";
        dumpVarianceReducedParticleMoments(particle_operations,particle_container,low_fidelity_states[i],*ops,variance_reduced_prefix, i_timestep, end_time);
        const std::string low_fidelity_prefix = "low_fidelity_moments";
        dumpLowFidelityMoments(low_fidelity_states[i],*ops,low_fidelity_prefix, i_timestep, end_time);
      }

      IntegratedCharge integrated_charge = particle_operations.assembleCharge(particle_container);
      electrostatic_field_operations->fieldSolve(particle_electrostatic_field_state, integrated_charge);
      mesh_data_writer.output(
        particle_electrostatic_field_state,
        low_fidelity_field_states,
        low_fidelity_states,
        i_timestep,
        end_time);

      IntegratedCharge variance_reduced_integrated_charge = particle_operations.assembleVarianceReducedCharge(particle_container,low_fidelity_states[0],*low_fidelity_operations[0]);
      electrostatic_field_operations->fieldSolve(variance_reduced_particle_electrostatic_field_state, variance_reduced_integrated_charge);
      variance_reduced_mesh_data_writer.output(
        variance_reduced_particle_electrostatic_field_state,
        low_fidelity_field_states,
        low_fidelity_states,
        i_timestep,
        end_time);

      text_data_writer.output(
        particle_electrostatic_field_state,
        low_fidelity_field_states,
        *electrostatic_field_operations,
        low_fidelity_states,
        low_fidelity_operations,
        i_timestep,
        end_time);

      variance_reduced_text_data_writer.output(
        variance_reduced_particle_electrostatic_field_state,
        low_fidelity_field_states,
        *electrostatic_field_operations,
        low_fidelity_states,
        low_fidelity_operations,
        i_timestep,
        end_time);
    }
  }

}

}
