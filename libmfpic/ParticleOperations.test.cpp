#include <libmfpic/DGEulerAssembly.hpp>
#include <libmfpic/DGEulerOperations.hpp>
#include <libmfpic/DGEulerInitialConditionsFactory.hpp>
#include <libmfpic/Discretization.hpp>
#include <libmfpic/BuildVarianceReductionParametersFromYaml.hpp>
#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/LoadParticles.hpp>
#include <libmfpic/LowFidelityState.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>
#include <libmfpic/SourcesFactory.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>
#include <mfem/mfem.hpp>

#include <ranges>
#include <unordered_map>

namespace {

using namespace mfpic;

class ZeroElectromagneticFieldsEvaluator : public ElectromagneticFieldsEvaluator {
public:
  virtual mfem::Vector getEFieldAt(const mfem::Vector&, const int) const {
    return mfem::Vector({0.0, 0.0, 0.0});
  }
  virtual mfem::Vector getBFieldAt(const mfem::Vector&, const int) const {
    return mfem::Vector({0.0, 0.0, 0.0});
  }
};

const Species default_species{.charge = 1.0, .mass = 1.0};
const std::vector<std::shared_ptr<ParticleBoundaryFactory>> empty_particle_boundary_factory_list;
const std::shared_ptr<ParticleBoundaryFactory> default_reflecting_particle_boundary_factory
  = std::make_shared<ReflectingParticleBoundaryFactory>();
const std::unordered_map<std::string, Species> one_species {{"one", default_species}};

struct MomentsInCell {
  double number_density;
  mfem::Vector bulk_velocity;
  double temperature;
};

TEST(ParticleOperations, AccelerateDoesNothingWithZeroFields) {
  ParticleContainer static_particles;
  static_particles.addParticle(Particle{
    .element = 0,
    .species = default_species
  });

  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems,1.0);
  constexpr int order = 1;
  Discretization disc(&mesh,order);

  ParticleOperations particle_operations(
    disc,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer accelerated_particles = particle_operations.accelerate(1.0, static_particles, ZeroElectromagneticFieldsEvaluator());

  for (Particle accelerated_particle : accelerated_particles) {
    EXPECT_DOUBLE_EQ(accelerated_particle.velocity[0], 0.0);
  }
}

TEST(ParticleOperations, AccelerateDoesNothingWithDeadParticles) {
  ParticleContainer static_particles;
  static_particles.addParticle(Particle{
    .element = 0,
    .species = default_species,
    .is_alive = false,
  });

  constexpr int num_elems = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems,1.0);
  constexpr int order = 1;
  Discretization disc(&mesh,order);

  ParticleOperations particle_operations(
    disc,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer accelerated_particles = particle_operations.accelerate(1.0, static_particles, ZeroElectromagneticFieldsEvaluator());

  for (Particle accelerated_particle : accelerated_particles) {
    EXPECT_DOUBLE_EQ(accelerated_particle.velocity[0], 0.0);
  }
}

TEST(ParticleOperations, AssembleChargeIgnoresDeadParticles) {
  constexpr int num_elements = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elements);

  constexpr int order = 1;
  Discretization discretization(&mesh, order);

  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector({0.125,0.0,0.0}),
    .velocity = mfem::Vector({0.0,0.0,0.0}),
    .element = 0,
    .species = default_species,
    .weight = 1.0,
    .is_alive = false,
  });


  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  IntegratedCharge charge_state = particle_operations.assembleCharge(particles);

  for (int dof = 0; dof < discretization.getFeSpace().GetNDofs(); dof++) {
    EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(0), 0.0);
  }
}

TEST(ParticleOperations, AssembleChargeWorksin1D) {
  constexpr int num_elements = 4;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elements,2.0);

  constexpr int order = 1;
  Discretization discretization(&mesh,order);

  //Place one particle in the middle of each element
  ParticleContainer static_particles;
  for (int i = 0; i < num_elements; i++) {
    mfem::ElementTransformation *element_transformation = mesh.GetElementTransformation(i);
    mfem::IntegrationPoint integration_point;
    integration_point.Set1w(0.5,1.0);
    mfem::Vector physical_point(1);
    element_transformation->Transform(integration_point, physical_point);

    static_particles.addParticle(Particle{
      .position = mfem::Vector({physical_point(0),0.0,0.0}),
      .velocity = mfem::Vector({0.0,0.0,0.0}),
      .element = i,
      .species = default_species,
      .weight = 1.0,
    });
  }


  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  IntegratedCharge charge_state = particle_operations.assembleCharge(static_particles);

  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(0), 0.5);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(1), 1.0);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(2), 1.0);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(3), 1.0);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(4), 0.5);
}

TEST(ParticleOperations, AssembleChargeWorksin2D) {
  constexpr int num_elems_per_dim = 2;
  constexpr bool generate_edges = true;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );

  // Mesh element ordering
  // *---*---*
  // | 2 | 3 |
  // *---*---*
  // | 0 | 1 |
  // *---*---*

  constexpr int order = 1;
  Discretization discretization(&mesh,order);

  ParticleContainer static_particles;
  //Place one particle in the middle of each element with different charges
  for (int i = 0; i < num_elems_per_dim*num_elems_per_dim; i++) {
    double charge;
    if (i==0) charge = 1.0;
    if (i==1) charge = 2.0;
    if (i==2) charge = 3.0;
    if (i==3) charge = 4.0;

    mfem::ElementTransformation *element_transformation = mesh.GetElementTransformation(i);
    mfem::IntegrationPoint integration_point;
    integration_point.Set2w(0.5,0.5,1.0);
    mfem::Vector physical_point(2);
    element_transformation->Transform(integration_point, physical_point);

    static_particles.addParticle(Particle{
      .position = mfem::Vector({physical_point(0),physical_point(1),0.0}),
      .velocity = mfem::Vector({0.0,0.0,0.0}),
      .element = i,
      .species = {.charge = charge, .mass = 1.0},
      .weight = 1.0,
    });
  }

  //One particle on 0-1 element boundary
  mfem::ElementTransformation *element_transformation = mesh.GetElementTransformation(0);
  mfem::IntegrationPoint integration_point;
  integration_point.Set2w(1.0,0.5,1.0);
  mfem::Vector physical_point(2);
  element_transformation->Transform(integration_point, physical_point);

  static_particles.addParticle(Particle{
    .position = mfem::Vector({physical_point(0),physical_point(1),0.0}),
    .velocity = mfem::Vector({0.0,0.0,0.0}),
    .element = 0,
    .species = {.charge = 1.0, .mass = 1.0},
    .weight = 1.0,
  });


  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  IntegratedCharge charge_state = particle_operations.assembleCharge(static_particles);

  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(0), 0.25);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(1), 1.25);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(2), 0.5);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(3), 1.0);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(4), 3.0);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(5), 1.5);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(6), 0.75);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(7), 1.75);
  EXPECT_DOUBLE_EQ(charge_state.getIntegratedChargeValue(8), 1.0);
}

TEST(ParticleOperations, DeadParticlesDoNotMove) {
  constexpr int num_elems = 1;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.0, 0.0},
    .velocity = mfem::Vector{0.5, 0.0, 0.0},
    .element = 0,
    .species = default_species,
    .is_alive = false,
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.25);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, ParticleCanMoveWithinAnElementIn1D) {
  constexpr int num_elems = 1;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.0, 0.0},
    .velocity = mfem::Vector{0.5, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossOneElementInterfaceIn1D) {
  constexpr int num_elems = 2;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.0, 0.0},
    .velocity = mfem::Vector{0.5, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
    EXPECT_EQ(moved_particle.element, 1);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossMultipleElementInterfacesIn1D) {
  constexpr int num_elems = 6;
  constexpr double dx = 1.0 / num_elems;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{dx / 2.0, 0.0, 0.0},
    .velocity = mfem::Vector{(num_elems - 1) * dx, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 1.0 - dx / 2.0);
    EXPECT_EQ(moved_particle.element, num_elems - 1);
  }
}

TEST(ParticleOperations, InitialElementIsArbitraryWhenParticleStartsOnElementInterfaceIn1D) {
  constexpr int num_elems = 2;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  for (int starting_element = 0; starting_element < num_elems; starting_element++) {
    ParticleContainer particles;
    particles.addParticle(Particle{
      .position = mfem::Vector{0.5, 0.0, 0.0},
      .velocity = mfem::Vector{0.25, 505.0, 87108.0},
      .element = starting_element,
      .species = default_species
    });

    constexpr double dt = 1.0;
    particles = particle_operations.move(dt, particles);

    for (Particle& moved_particle : particles) {
      EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
      EXPECT_EQ(moved_particle.element, 1);
    }
  }
}

TEST(ParticleOperations, ParticleMovesAcrossPeriodicBoundariesIn1D) {
  constexpr int num_elems = 3;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({mfem::Vector({1.0})});
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);
  constexpr int order = 1;
  Discretization discretization(&periodic_mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.1, 0.0, 0.0},
    .velocity = mfem::Vector{1.2, 505.0, 87108.0},
    .element = 0,
  });
  particles.addParticle(Particle{
    .position = mfem::Vector{0.1, 0.0, 0.0},
    .velocity = mfem::Vector{-0.8, 505.0, 87108.0},
    .element = 0,
  });

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.3);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, ParticleCanMoveWithinAnElementIn2D) {
  constexpr int num_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems,
    num_elems,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.25, 7.0},
    .velocity = mfem::Vector{0.5, 0.5, 1213542342.0},
    .element = 0,
    .species = default_species
  });

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.75);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossOneElementInterfaceIn2D) {
  constexpr int num_elems_per_side = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::TRIANGLE;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_side,
    num_elems_per_side,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.1, 0.9, 0.0},
    .velocity = mfem::Vector{0.8, -0.8, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.9);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.1);
    EXPECT_EQ(moved_particle.element, 1);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossMultipleElementInterfacesIn2D) {
  constexpr int num_x_elems = 6;
  constexpr int num_y_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_x_elems,
    num_y_elems,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  ParticleContainer particles;
  constexpr double dx = domain_side_length / num_x_elems;
  particles.addParticle(Particle{
    .position = mfem::Vector{dx / 2.0, domain_side_length / 2.0, 0.0},
    .velocity = mfem::Vector{(num_x_elems - 1) * dx, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], domain_side_length - dx / 2.0);
    EXPECT_EQ(moved_particle.element, num_x_elems - 1);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossCornerInterfaceIn2D) {
  constexpr int num_elems_per_dim = 2;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  ParticleContainer particles;
  constexpr double dx = domain_side_length / num_elems_per_dim;
  particles.addParticle(Particle{
    .position = mfem::Vector{dx / 2.0, dx / 2.0, 0.0},
    .velocity = mfem::Vector{dx, dx, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], domain_side_length - dx / 2.0);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], domain_side_length - dx / 2.0);
    EXPECT_EQ(moved_particle.element, 3);
  }
}

TEST(ParticleOperations, InitialElementIsArbitraryWhenParticleStartsOnElementCornerInterfaceIn2D) {
  constexpr int num_elems_per_dim = 2;
  constexpr mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    generate_edges,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dx = domain_side_length / num_elems_per_dim;
  for (int starting_element = 0; starting_element < num_elems_per_dim*num_elems_per_dim; starting_element++) {
    ParticleContainer particles;
    particles.addParticle(Particle{
      .position = mfem::Vector{dx, dx, 0.0},
      .velocity = mfem::Vector{dx/2.0, -dx/2.0, 505.0},
      .element = starting_element,
      .species = default_species
    });

    constexpr double dt = 1.0;
    particles = particle_operations.move(dt, particles);

    for (Particle& moved_particle : particles) {
      EXPECT_DOUBLE_EQ(moved_particle.position[0], 1.5 * dx);
      EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.5 * dx);
      EXPECT_EQ(moved_particle.element, 1);
    }
  }
}

TEST(ParticleOperations, ParticleMovesAcrossPeriodicBoundariesIn2D) {
  constexpr int num_elems_per_dim = 3;
  constexpr bool generate_edges = true;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  for (const mfem::Element::Type element_type : {mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL}) {
    mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian2D(
      num_elems_per_dim,
      num_elems_per_dim,
      element_type,
      generate_edges,
      domain_side_length,
      domain_side_length,
      space_filling_curve_ordering
    );
    std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({
      mfem::Vector({1.0, 0.0}),
      mfem::Vector({0.0, 1.0})
    });
    mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);
    constexpr int order = 1;
    Discretization discretization(&periodic_mesh, order);
    ParticleOperations particle_operations(
      discretization,
      empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory,
      one_species
    );

    ParticleContainer particles;
    particles.addParticle(Particle{
      .position = mfem::Vector{0.125, 0.25, 0.0},
      .velocity = mfem::Vector{1.0, 1.0, 87108.0},
      .element = 0,
    });
    particles.addParticle(Particle{
      .position = mfem::Vector{0.125, 0.25, 0.0},
      .velocity = mfem::Vector{-1.0, -1.0, 87108.0},
      .element = 0,
    });

    constexpr double dt = 1.0;
    particles = particle_operations.move(dt, particles);

    for (Particle& moved_particle : particles) {
      EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.125);
      EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.25);
      EXPECT_EQ(moved_particle.element, 0);
    }
  }
}

TEST(ParticleOperations, ParticleCanMoveWithinAnElementIn3D) {
  constexpr int num_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems,
    num_elems,
    num_elems,
    element_type
  );
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.25, 0.25},
    .velocity = mfem::Vector{0.5, 0.5, 0.5},
    .element = 0,
    .species = default_species
  });

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.75);
    EXPECT_DOUBLE_EQ(moved_particle.position[2], 0.75);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossOneElementInterfaceIn3D) {
  constexpr int num_x_elems = 2;
  constexpr int num_yz_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_x_elems,
    num_yz_elems,
    num_yz_elems,
    element_type
  );
  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.25, 0.5, 0.5},
    .velocity = mfem::Vector{0.5, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.75);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.5);
    EXPECT_DOUBLE_EQ(moved_particle.position[2], 0.5);
    EXPECT_EQ(moved_particle.element, 1);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossMultipleElementInterfacesIn3D) {
  constexpr int num_x_elems = 6;
  constexpr int num_yz_elems = 1;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_x_elems,
    num_yz_elems,
    num_yz_elems,
    element_type
  );
  ParticleContainer particles;
  constexpr double dx = 1.0 / num_x_elems;
  particles.addParticle(Particle{
    .position = mfem::Vector{dx / 2.0, 0.5, 0.5},
    .velocity = mfem::Vector{(num_x_elems - 1) * dx, 0.0, 0.0},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 1.0 - dx / 2.0);
    EXPECT_EQ(moved_particle.element, num_x_elems - 1);
  }
}

TEST(ParticleOperations, ParticleCanMoveAcrossCornerInterfaceIn3D) {
  constexpr int num_elems_per_dim = 2;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    domain_side_length,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  ParticleContainer particles;
  constexpr double dx = domain_side_length / num_elems_per_dim;
  particles.addParticle(Particle{
    .position = mfem::Vector{dx / 2.0, dx / 2.0, dx / 2.0},
    .velocity = mfem::Vector{dx, dx, dx},
    .element = 0,
    .species = default_species
  });
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], domain_side_length - dx / 2.0);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], domain_side_length - dx / 2.0);
    EXPECT_DOUBLE_EQ(moved_particle.position[2], domain_side_length - dx / 2.0);
    EXPECT_EQ(moved_particle.element, 7);
  }
}

TEST(ParticleOperations, InitialElementIsArbitraryWhenParticleStartsOnElementCornerInterfaceIn3D) {
  constexpr int num_elems_per_dim = 2;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    domain_side_length,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dx = domain_side_length / num_elems_per_dim;
  for (int starting_element = 0; starting_element < num_elems_per_dim*num_elems_per_dim*num_elems_per_dim; starting_element++) {
    ParticleContainer particles;
    particles.addParticle(Particle{
      .position = mfem::Vector{dx, dx, dx},
      .velocity = mfem::Vector{dx/2.0, dx/2.0, dx/2.0},
      .element = starting_element,
      .species = default_species
    });

    constexpr double dt = 1.0;
    particles = particle_operations.move(dt, particles);

    for (Particle& moved_particle : particles) {
      EXPECT_DOUBLE_EQ(moved_particle.position[0], 1.5 * dx);
      EXPECT_DOUBLE_EQ(moved_particle.position[1], 1.5 * dx);
      EXPECT_DOUBLE_EQ(moved_particle.position[2], 1.5 * dx);
      EXPECT_EQ(moved_particle.element, 7);
    }
  }
}

TEST(ParticleOperations, ParticleMovesAcrossPeriodicBoundariesIn3D) {
  constexpr int num_elems_per_dim = 3;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  constexpr double domain_side_length = 1.0;
  constexpr bool space_filling_curve_ordering = false;
  mfem::Mesh non_periodic_mesh = mfem::Mesh::MakeCartesian3D(
    num_elems_per_dim,
    num_elems_per_dim,
    num_elems_per_dim,
    element_type,
    domain_side_length,
    domain_side_length,
    domain_side_length,
    space_filling_curve_ordering
  );
  std::vector<int> periodic_vertex_mapping = non_periodic_mesh.CreatePeriodicVertexMapping({
    mfem::Vector({1.0, 0.0, 0.0}),
    mfem::Vector({0.0, 1.0, 0.0}),
    mfem::Vector({0.0, 0.0, 1.0}),
  });
  mfem::Mesh periodic_mesh = mfem::Mesh::MakePeriodic(non_periodic_mesh, periodic_vertex_mapping);
  constexpr int order = 1;
  Discretization discretization(&periodic_mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  ParticleContainer particles;
  particles.addParticle(Particle{
    .position = mfem::Vector{0.2, 0.2, 0.2},
    .velocity = mfem::Vector{1.1, 1.1, 1.1},
    .element = 0,
  });
  particles.addParticle(Particle{
    .position = mfem::Vector{0.2, 0.2, 0.2},
    .velocity = mfem::Vector{-0.9, -0.9, -0.9},
    .element = 0,
  });

  constexpr double dt = 1.0;
  particles = particle_operations.move(dt, particles);

  for (Particle& moved_particle : particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], 0.3);
    EXPECT_DOUBLE_EQ(moved_particle.position[1], 0.3);
    EXPECT_DOUBLE_EQ(moved_particle.position[2], 0.3);
    EXPECT_EQ(moved_particle.element, 0);
  }
}

TEST(ParticleOperations, PerParticleTimeIntervalsMoveEachParticleByItsOwnInterval) {
  constexpr int num_elems = 6;
  constexpr double dx = 1.0 / num_elems;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  const std::vector<double> per_particle_dt{0.0, 0.5, 2.0};
  for (int i = 0; i < std::ssize(per_particle_dt); i++) {
    particles.addParticle(Particle{
      .position = mfem::Vector{dx / 2.0, 0.0, 0.0},
      .velocity = mfem::Vector{2.0 * dx, 0.0, 0.0},
      .element = 0,
      .species = default_species
    });
  }
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  const ParticleContainer moved_particles = particle_operations.move(per_particle_dt, particles);

  const std::vector<double> expected_positions{dx / 2.0, 3.0 * dx / 2.0, 9.0 * dx / 2.0};
  const std::vector<int> expected_elements{0, 1, 4};
  int i_particle = 0;
  for (const Particle& moved_particle : moved_particles) {
    EXPECT_DOUBLE_EQ(moved_particle.position[0], expected_positions[i_particle]);
    EXPECT_EQ(moved_particle.element, expected_elements[i_particle]);
    i_particle++;
  }
}

TEST(ParticleOperations, PerParticleTimeIntervalsMatchUniformMoveWhenAllEqual) {
  constexpr int num_elems = 6;
  constexpr double dx = 1.0 / num_elems;
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  ParticleContainer particles;
  const std::vector<double> velocities{3.0 * dx, -4.0 * dx, 7.5 * dx};
  for (const double velocity : velocities) {
    particles.addParticle(Particle{
      .position = mfem::Vector{2.5 * dx, 0.0, 0.0},
      .velocity = mfem::Vector{velocity, 0.0, 0.0},
      .element = 2,
      .species = default_species
    });
  }
  constexpr int order = 1;
  Discretization discretization(&mesh, order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  constexpr double dt = 1.0;
  const ParticleContainer uniformly_moved_particles = particle_operations.move(dt, particles);
  const ParticleContainer individually_moved_particles = particle_operations.move(
    std::vector<double>(velocities.size(), dt),
    particles
  );

  for (const auto& [individually_moved_particle, uniformly_moved_particle] : std::views::zip(individually_moved_particles, uniformly_moved_particles)) {
    EXPECT_DOUBLE_EQ(individually_moved_particle.position[0], uniformly_moved_particle.position[0]);
    EXPECT_DOUBLE_EQ(individually_moved_particle.velocity[0], uniformly_moved_particle.velocity[0]);
    EXPECT_EQ(individually_moved_particle.element, uniformly_moved_particle.element);
  }
}

TEST(ParticleOperations, VarianceReducedChargeIsExactForMaxwellianIn3D) {
  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .name="one"};
  const std::unordered_map<std::string, Species> species_map {{"one", species}};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({1.0,2.0,3.0});
  constexpr int num_particles = 20000;
  RandomNumberGenerator generator;

  const int num_elems = 5;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
    num_elems,
    num_elems,
    num_elems,
    element_type
  ));

  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, temperature,bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);
  IntegratedCharge low_fidelity_charge_state = dg_euler_operations.assembleCharge(low_fidelity_state);
  mfem::Vector integrated_charge_vector = low_fidelity_charge_state.getIntegratedCharge();

  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
  };

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(species, source_state_parameters, num_particles),
    generator,
    mesh,
    3
  );

  ParticleOperations particle_operations(
    charge_discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    species_map
  );

  VarianceReductionParameters variance_reduction_parameters;
  variance_reduction_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  particle_operations.setVarianceReductionParameters(variance_reduction_parameters);
  IntegratedCharge variance_reduced_charge_state = particle_operations.assembleVarianceReducedCharge(particles,low_fidelity_state,dg_euler_operations);

  for (int dof = 0; dof < charge_discretization.getFeSpace().GetNDofs(); dof++) {
    EXPECT_DOUBLE_EQ(variance_reduced_charge_state.getIntegratedChargeValue(dof),integrated_charge_vector(dof));
  }
}

ParticleContainer takePrefix(const ParticleContainer& all, int num_particles) {
  ParticleContainer out;
  num_particles = std::max(0, std::min(num_particles, all.numParticles()));
  int count = 0;
  for (auto it = all.begin(); it != all.end() && count < num_particles; ++it, ++count) {
    out.addParticle(*it);
  }
  return out;
}

TEST(ParticleOperations, VarianceReducedChargeAndPICChargeConvergeForKappa) {

  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .name="one"};
  const std::unordered_map<std::string, Species> species_map {{"one", species}};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({1.0, 0.0, 0.0});

  const std::vector<int> num_particles_list = {10,100,1000,10000,100000};

  const int num_elems = 5;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  std::shared_ptr<mfem::Mesh> mesh =
      std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr =
      std::make_shared<DGEulerAssembly>(finite_element_space, species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(
      species, number_density, temperature, bulk_velocity));

  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);
  IntegratedCharge low_fidelity_charge_state =
      dg_euler_operations.assembleCharge(low_fidelity_state);
  mfem::Vector integrated_charge_vector = low_fidelity_charge_state.getIntegratedCharge();

  ParticleOperations particle_operations(
      charge_discretization,
      empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory,
      species_map
    );

  double prev_max_rel_error = std::numeric_limits<double>::infinity();

  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
    .kappa = 2.0
  };
  RandomNumberGenerator gen_for_n(12345);

  ParticleContainer particles_all = loadParticles(
    ConstantSourceParameters(species, source_state_parameters, num_particles_list[4]),
    gen_for_n,
    mesh,
    3
  );

  for (int num_particles : num_particles_list) {
    ParticleContainer particles = takePrefix(particles_all, num_particles);
    IntegratedCharge variance_reduced_charge_state =
        particle_operations.assembleVarianceReducedCharge(
            particles, low_fidelity_state, dg_euler_operations);
    IntegratedCharge charge_state = particle_operations.assembleCharge(particles);

    double max_rel_error = 0.0;
    for (int dof = 0; dof < charge_discretization.getFeSpace().GetNDofs(); dof++) {
      const double vr = variance_reduced_charge_state.getIntegratedChargeValue(dof);
      const double st = charge_state.getIntegratedChargeValue(dof);

      const double denom = std::max(std::abs(st), 1e-300);
      const double rel_error = std::abs(vr - st) / denom;

      max_rel_error = std::max(max_rel_error, rel_error);
    }

    EXPECT_LE(max_rel_error, prev_max_rel_error)
        << "Expected max relative error to decrease with num_particles. "
        << "prev=" << prev_max_rel_error << ", current=" << max_rel_error
        << ", num_particles=" << num_particles;

    prev_max_rel_error = max_rel_error;
  }
}

TEST(ParticleOperations, TemperatureIsZeroWithZeroOrOneParticles) {
  constexpr int num_elems = 1;
  constexpr double domain_length = 1.0;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, domain_length));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  for (int num_particles = 0; num_particles <= 1; num_particles++) {
    ParticleContainer particles;
    for (int iparticle = 0; iparticle < num_particles; iparticle++) {
      particles.addParticle(Particle{.velocity = mfem::Vector({1.0, 2.0, 3.0}), .species = default_species, .weight = 1.0});
    }

    ParticleOperations particle_operations(
      discretization,
      empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory,
      one_species
    );

    ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
    auto temperature_per_species_per_element = particle_moments.temperature;
    const double computed_temperature = temperature_per_species_per_element.at(default_species)(0);

    EXPECT_DOUBLE_EQ(0.0, computed_temperature);
  }
}

TEST(ParticleOperations, TemperatureCalculationIsUnbiasedWithVariableParticleWeights) {
  const int num_elems = 1;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);
  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  RandomNumberGenerator generator;
  constexpr double nominal_temperature = 1.0e25;
  constexpr int num_repetitions = 100;
  double average_computed_temperature = 0.0;
  for (int rep = 0; rep < num_repetitions; rep++) {
    ParticleContainer particles;
    const mfem::Vector nominal_bulk_velocity({0.0, 0.0, 0.0});
    constexpr double weight_lower_bound = 0.0;
    constexpr double weight_upper_bound = 2.0;
    std::uniform_real_distribution<double> weight_distribution(weight_lower_bound, weight_upper_bound);
    constexpr int num_particles_per_repetition = 4;
    for (int iparticle = 0; iparticle < num_particles_per_repetition; iparticle++) {
      particles.addParticle(Particle{
        .velocity = generateMaxwellianVelocity(nominal_bulk_velocity, nominal_temperature, default_species.mass, generator),
        .species = default_species,
        .weight = weight_distribution(generator),
      });
    }

    ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
    auto temperature_per_species_per_element = particle_moments.temperature;
    const double computed_temperature = temperature_per_species_per_element.at(default_species)(0);
    average_computed_temperature += computed_temperature / num_repetitions;
  }

  constexpr double relative_tolerance = 0.01;
  EXPECT_NEAR(average_computed_temperature, nominal_temperature, relative_tolerance * nominal_temperature);
}

TEST(ParticleOperations, ParticleMomentsCorrectForMaxwellian) {
  const int num_elems = 1;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  const mfem::Vector nominal_bulk_velocity({300.0,0.0,0.0});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
  };
  constexpr int num_particles = 200000;
  RandomNumberGenerator generator;

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(default_species, source_state_parameters, num_particles),
    generator,
    mesh,
    3
  );

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    one_species
  );

  ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
  auto pic_number_density = particle_moments.number_density;
  auto pic_bulk_velocity = particle_moments.bulk_velocity;
  auto pic_temperature = particle_moments.temperature;

  mfem::Vector computed_bulk_velocity;
  pic_bulk_velocity.at(default_species).GetColumn(0, computed_bulk_velocity);
  double computed_number_density = pic_number_density.at(default_species)(0);
  double computed_temperature = pic_temperature.at(default_species)(0);

  EXPECT_NEAR(computed_number_density, number_density, 1e-10);

  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(computed_bulk_velocity[i], nominal_bulk_velocity[i], 1e-10);
  }

  EXPECT_NEAR(computed_temperature, temperature, 5e-3*temperature);
}

TEST(ParticleOperations, TemperatureCorrectForMaxwellianWithAnyNumberOfVelocityDimensions) {
  const int num_elems = 1;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  const Species electron{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .name = "electron"};
  const std::unordered_map<std::string, Species> electron_species_map{{"electron", electron}};

  // Bulk velocity is nonzero in every component so components without thermal spread must be excluded.
  const mfem::Vector nominal_bulk_velocity({1.0e5, -2.0e5, 5.0e4});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = nominal_bulk_velocity,
    .temperature = temperature,
  };
  constexpr int num_particles = 200000;

  for (int velocity_dims = 1; velocity_dims <= 3; ++velocity_dims) {
    RandomNumberGenerator generator;
    ParticleContainer particles = loadParticles(
      ConstantSourceParameters(electron, source_state_parameters, num_particles),
      generator,
      mesh,
      velocity_dims
    );

    ParticleOperations particle_operations(
      discretization,
      empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory,
      electron_species_map,
      velocity_dims
    );

    ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
    const double computed_temperature = particle_moments.temperature.at(electron)(0);

    EXPECT_NEAR(computed_temperature, temperature, 1e-2 * temperature) << "velocity_dims = " << velocity_dims;
  }
}

TEST(ParticleOperations, ParticleMomentsCorrectForKnownParticles) {

  const int num_elems_in_one_direction = 2;
  const int num_elems = num_elems_in_one_direction * num_elems_in_one_direction * num_elems_in_one_direction;
  const double lx = 1.0, ly = 2.0, lz = 3.0;
  const double dx = lx / num_elems_in_one_direction, dy = ly / num_elems_in_one_direction, dz = lz / num_elems_in_one_direction;
  const int num_particles_in_cell = 4;
  const double bias_correction = num_particles_in_cell / (num_particles_in_cell - 1.);
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(num_elems_in_one_direction, num_elems_in_one_direction, num_elems_in_one_direction, mfem::Element::HEXAHEDRON, lx, ly, lz));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  Species electron {.mass = constants::electron_mass};
  Species proton {.mass = constants::proton_mass};

  const std::unordered_map<std::string, Species> two_species {{"electron", electron}, {"proton", proton}};

  ParticleContainer particles;
  std::vector<MomentsInCell> exact_moments_electron;
  std::vector<MomentsInCell> exact_moments_proton;
  std::vector<int> elem_ids;
  const double electron_fac = electron.mass / (3. * constants::boltzmann_constant);
  const double proton_fac = proton.mass / (3. * constants::boltzmann_constant);

  for (int i = 0; i < num_elems_in_one_direction; ++i) {
    for (int j = 0; j < num_elems_in_one_direction; ++j) {
      for (int k = 0; k < num_elems_in_one_direction; ++k) {
        double sum_of_weights = 0.;
        double sum_of_inverse = 0.;
        const int i_elem = i + j * num_elems_in_one_direction + k * num_elems_in_one_direction * num_elems_in_one_direction;
        const mfem::Vector center {dx * (i + 1./2.), dy * (j + 1./2.), dz * (k + 1./2.)};
        for (int i_part = 0; i_part < num_particles_in_cell; ++i_part) {
          const double x = dx * i + (i_part + 1) * dx / (num_particles_in_cell + 1);
          const double y = dy * j + (i_part + 1) * dy / (num_particles_in_cell + 1);
          const double z = dz * k + (i_part + 1) * dz / (num_particles_in_cell + 1);
          //const double weight = std::sqrt(x*x + y*y + z*z); // non-constant weight not supported
          constexpr double weight = 12.345;
          const mfem::Vector position {x, y, z};
          mfem::Vector velocity = center;
          mfem::Vector offset(3);
          if (i_part == 0) {
            offset = mfem::Vector {1., -1., 1.};
          } else if (i_part == 1) {
            offset = mfem::Vector {-1., 1., -1.};
          } else if (i_part == 2) {
            offset = mfem::Vector {1., 1., -1.};
          } else if (i_part == 3) {
            offset = mfem::Vector {-1., -1., 1.};
          }
          offset /= weight;

          velocity += offset;

          particles.addParticle({.position = position, .velocity = velocity, .element = i_elem, .species = electron, .weight = weight});
          particles.addParticle({.position = position, .velocity = velocity, .element = i_elem, .species = proton, .weight = weight});

          sum_of_weights += weight;
          sum_of_inverse += 1./weight;

        }

        mfem::Vector bulk = center;
        exact_moments_electron.emplace_back(MomentsInCell{.number_density = sum_of_weights / (dx * dy * dz), .bulk_velocity = bulk, .temperature = bias_correction * 3. * sum_of_inverse / sum_of_weights * electron_fac});
        exact_moments_proton.emplace_back(MomentsInCell{.number_density = sum_of_weights / (dx * dy * dz), .bulk_velocity = bulk, .temperature = bias_correction * 3. * sum_of_inverse / sum_of_weights * proton_fac});
        elem_ids.push_back(i_elem);
      }
    }
  }

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    two_species
  );

  auto check_moments = [&] (MomentsInCell exact, const Species & species, int cell_id) {

    ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
    auto pic_number_density = particle_moments.number_density;
    auto pic_bulk_velocity = particle_moments.bulk_velocity;
    auto pic_temperature = particle_moments.temperature;
    MomentsInCell computed {.number_density = pic_number_density.at(species)(cell_id),
                            .bulk_velocity = mfem::Vector(pic_bulk_velocity.at(species).GetColumn(cell_id),3),
                            .temperature = pic_temperature.at(species)(cell_id) };

    EXPECT_NEAR(exact.number_density, computed.number_density, 1e-12);

    for (int i = 0; i < 3; i++) {
      EXPECT_NEAR(exact.bulk_velocity[i], computed.bulk_velocity[i], 1e-12);
    }

    EXPECT_NEAR(exact.temperature, computed.temperature, 1e-12);
  };

  for (int i_elem = 0; i_elem < num_elems; ++i_elem) {
    check_moments(exact_moments_electron[i_elem], electron, elem_ids[i_elem]);
    check_moments(exact_moments_proton[i_elem], proton, elem_ids[i_elem]);
  }

};


TEST(ParticleOperations, VarianceReducedMomentsAreExactForMaxwellianIn3D) {
  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({293.0,581.0,902.0});
  constexpr int num_particles = 20000;
  RandomNumberGenerator generator;

  const int num_elems = 5;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  constexpr mfem::Element::Type element_type = mfem::Element::HEXAHEDRON;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
    num_elems,
    num_elems,
    num_elems,
    element_type
  ));

  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(finite_element_space, species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, temperature,bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);

  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
  };

  ParticleContainer particles = loadParticles(
    ConstantSourceParameters(species, source_state_parameters, num_particles),
    generator,
    mesh,
    3
  );

  ParticleOperations particle_operations(
    charge_discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    species_map
  );
  VarianceReductionParameters variance_reduction_parameters;
  variance_reduction_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  variance_reduction_parameters.cell_limiting = VarianceReductionParameters::CellLimiting::None;
  variance_reduction_parameters.use_variance_reduced_electric_field = false;
  particle_operations.setVarianceReductionParameters(variance_reduction_parameters);

  ParticleMoments variance_reduced_particle_moments = particle_operations.getVarianceReducedParticleMoments(particles,low_fidelity_state,dg_euler_operations,true);
  auto variance_reduced_number_density = variance_reduced_particle_moments.number_density;
  auto variance_reduced_bulk_velocity = variance_reduced_particle_moments.bulk_velocity;
  auto variance_reduced_temperature = variance_reduced_particle_moments.temperature;

  for (int elem_id = 0; elem_id < num_elems; ++elem_id)
  {
    EXPECT_DOUBLE_EQ(variance_reduced_number_density.at(species)(elem_id),number_density);
    EXPECT_DOUBLE_EQ(variance_reduced_bulk_velocity.at(species)(0,elem_id),bulk_velocity(0));
    EXPECT_DOUBLE_EQ(variance_reduced_bulk_velocity.at(species)(1,elem_id),bulk_velocity(1));
    EXPECT_DOUBLE_EQ(variance_reduced_bulk_velocity.at(species)(2,elem_id),bulk_velocity(2));
    EXPECT_DOUBLE_EQ(variance_reduced_temperature.at(species)(elem_id),temperature);
  }
}

TEST(ParticleOperations, VarianceReducedMomentsAndPICMomentsConvergeForKappa) {

  Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e22;
  constexpr double temperature = 300;
  mfem::Vector bulk_velocity({1.0, 0.0, 0.0});

  const std::vector<int> num_particles_list = {10,100,1000,10000,100000};

  const int num_elems = 5;
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  std::shared_ptr<mfem::Mesh> mesh =
      std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);

  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr =
      std::make_shared<DGEulerAssembly>(finite_element_space, species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(
      species, number_density, temperature, bulk_velocity));

  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);
  IntegratedCharge low_fidelity_charge_state =
      dg_euler_operations.assembleCharge(low_fidelity_state);
  mfem::Vector integrated_charge_vector = low_fidelity_charge_state.getIntegratedCharge();

  ParticleOperations particle_operations(
      charge_discretization,
      empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory,
      species_map
    );


  const SourceStateParameters source_state_parameters{
    .number_density = number_density,
    .bulk_velocity = bulk_velocity,
    .temperature = temperature,
    .kappa = 2.0
  };
  RandomNumberGenerator gen_for_n(12345);

  ParticleContainer particles_all = loadParticles(
    ConstantSourceParameters(species, source_state_parameters, num_particles_list[4]),
    gen_for_n,
    mesh,
    3
  );

  for (int num_particles : num_particles_list) {
    ParticleContainer particles = takePrefix(particles_all, num_particles);

    ParticleMoments variance_reduced_particle_moments = particle_operations.getVarianceReducedParticleMoments(particles,low_fidelity_state,dg_euler_operations,true);
    auto variance_reduced_number_density = variance_reduced_particle_moments.number_density;
    auto variance_reduced_bulk_velocity = variance_reduced_particle_moments.bulk_velocity;
    auto variance_reduced_temperature = variance_reduced_particle_moments.temperature;

    ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
    auto pic_number_density = particle_moments.number_density;
    auto pic_bulk_velocity = particle_moments.bulk_velocity;
    auto pic_temperature = particle_moments.temperature;

    double max_rel_error_number_density = 0.0;
    double max_rel_error_x_bulk_velocity = 0.0;
    double max_rel_error_y_bulk_velocity = 0.0;
    double max_rel_error_z_bulk_velocity = 0.0;
    double max_rel_error_temperature = 0.0;

    double prev_max_rel_error_number_density = std::numeric_limits<double>::infinity();
    double prev_max_rel_error_x_bulk_velocity = std::numeric_limits<double>::infinity();
    double prev_max_rel_error_y_bulk_velocity = std::numeric_limits<double>::infinity();
    double prev_max_rel_error_z_bulk_velocity = std::numeric_limits<double>::infinity();
    double prev_max_rel_error_temperature = std::numeric_limits<double>::infinity();

    for (int elem_id = 0; elem_id < num_elems; ++elem_id)
    {
      double vr = variance_reduced_number_density.at(species)(elem_id);
      double pic = pic_number_density.at(species)(elem_id);
      double denom = std::max(std::abs(pic), 1e-300);
      double rel_error_number_density = std::abs(vr - pic) / denom;
      max_rel_error_number_density = std::max(max_rel_error_number_density, rel_error_number_density);

      vr = variance_reduced_bulk_velocity.at(species)(0,elem_id);
      pic = pic_bulk_velocity.at(species)(0,elem_id);
      denom = std::max(std::abs(pic), 1e-300);
      double rel_error_x_bulk_velocity = std::abs(vr - pic) / denom;
      max_rel_error_x_bulk_velocity = std::max(max_rel_error_x_bulk_velocity, rel_error_x_bulk_velocity);

      vr = variance_reduced_bulk_velocity.at(species)(1,elem_id);
      pic = pic_bulk_velocity.at(species)(1,elem_id);
      denom = std::max(std::abs(pic), 1e-300);
      double rel_error_y_bulk_velocity = std::abs(vr - pic) / denom;
      max_rel_error_y_bulk_velocity = std::max(max_rel_error_y_bulk_velocity, rel_error_y_bulk_velocity);

      vr = variance_reduced_bulk_velocity.at(species)(2,elem_id);
      pic = pic_bulk_velocity.at(species)(2,elem_id);
      denom = std::max(std::abs(pic), 1e-300);
      double rel_error_z_bulk_velocity = std::abs(vr - pic) / denom;
      max_rel_error_z_bulk_velocity = std::max(max_rel_error_z_bulk_velocity, rel_error_z_bulk_velocity);

      vr = variance_reduced_temperature.at(species)(elem_id);
      pic = pic_temperature.at(species)(elem_id);
      denom = std::max(std::abs(pic), 1e-300);
      double rel_error_temperature = std::abs(vr - pic) / denom;
      max_rel_error_temperature = std::max(max_rel_error_temperature, rel_error_temperature);
    }

    EXPECT_LE(max_rel_error_number_density, prev_max_rel_error_number_density);
    prev_max_rel_error_number_density = max_rel_error_number_density;

    EXPECT_LE(max_rel_error_x_bulk_velocity, prev_max_rel_error_x_bulk_velocity);
    prev_max_rel_error_x_bulk_velocity = max_rel_error_x_bulk_velocity;

    EXPECT_LE(max_rel_error_y_bulk_velocity, prev_max_rel_error_y_bulk_velocity);
    prev_max_rel_error_y_bulk_velocity = max_rel_error_y_bulk_velocity;

    EXPECT_LE(max_rel_error_z_bulk_velocity, prev_max_rel_error_z_bulk_velocity);
    prev_max_rel_error_z_bulk_velocity = max_rel_error_z_bulk_velocity;

    EXPECT_LE(max_rel_error_temperature, prev_max_rel_error_temperature);
    prev_max_rel_error_temperature = max_rel_error_temperature;
  }
}

TEST(ParticleOperations, VarianceReducedMomentsMatchPICForSpeciesMissingLowFidelityState) {
  const int num_elems = 1;
  std::shared_ptr<mfem::Mesh> mesh =
    std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int dg_order = 0;
  constexpr int num_equations = 5;
  Discretization dg_discretization(mesh.get(), dg_order, FETypes::DG, num_equations);
  constexpr int charge_order = 1;
  Discretization charge_discretization(mesh.get(), charge_order, FETypes::HGRAD);

  const Species species_a{.charge = 1.0, .mass = 1.0};
  const Species species_b{.charge = 2.0, .mass = 2.0};

  const mfem::Vector nominal_bulk_velocity_a({300.0, 0.0, 0.0});
  const mfem::Vector nominal_bulk_velocity_b({-150.0, 0.0, 0.0});

  constexpr double temperature_a = 11600.0;
  constexpr double temperature_b = 8000.0;

  constexpr double number_density_a = 1.0e18;
  constexpr double number_density_b = 2.5e17;

  constexpr int num_particles_a = 15000;
  constexpr int num_particles_b = 15000;

  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::shared_ptr<DGEulerAssembly> operator_ptr =
      std::make_shared<DGEulerAssembly>(finite_element_space, species_a);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);

  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(
      species_a, number_density_a, temperature_a, nominal_bulk_velocity_a));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);

  RandomNumberGenerator generator(12345);
  const SourceStateParameters source_state_parameters_a{
    .number_density = number_density_a,
    .bulk_velocity  = nominal_bulk_velocity_a,
    .temperature    = temperature_a,
    .kappa = 2.0
  };
  const SourceStateParameters source_state_parameters_b{
    .number_density = number_density_b,
    .bulk_velocity  = nominal_bulk_velocity_b,
    .temperature    = temperature_b,
    .kappa = 3.0
  };

  ParticleContainer particles_a = loadParticles(
    ConstantSourceParameters(species_a, source_state_parameters_a, num_particles_a),
    generator,
    mesh,3);

  ParticleContainer particles_b = loadParticles(
    ConstantSourceParameters(species_b, source_state_parameters_b, num_particles_b),
    generator,
    mesh,3);

  ParticleContainer particles = particles_a;
  particles.addParticles(particles_b);

  const std::unordered_map<std::string, Species> two_species {{"species_a", species_a}, {"species_b", species_b}};
  ParticleOperations particle_operations(
    charge_discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    two_species);

  ParticleMoments variance_reduced_particle_moments = particle_operations.getVarianceReducedParticleMoments(particles,low_fidelity_state,dg_euler_operations,true);
  auto variance_reduced_number_density = variance_reduced_particle_moments.number_density;
  auto variance_reduced_bulk_velocity = variance_reduced_particle_moments.bulk_velocity;
  auto variance_reduced_temperature = variance_reduced_particle_moments.temperature;

  ParticleMoments particle_moments = particle_operations.getParticleMoments(particles);
  auto pic_number_density = particle_moments.number_density;
  auto pic_bulk_velocity = particle_moments.bulk_velocity;
  auto pic_temperature = particle_moments.temperature;

  const double standard_number_density_b = pic_number_density.at(species_b)(0);
  mfem::Vector standard_bulk_velocity_b;
  pic_bulk_velocity.at(species_b).GetColumn(0, standard_bulk_velocity_b);

  const double standard_temperature_b = pic_temperature.at(species_b)(0);

  const double variance_reduced_number_density_b = variance_reduced_number_density.at(species_b)(0);

  mfem::Vector variance_reduced_bulk_velocity_b;
  variance_reduced_bulk_velocity.at(species_b).GetColumn(0, variance_reduced_bulk_velocity_b);

  const double variance_reduced_temperature_b =
  variance_reduced_temperature.at(species_b)(0);

  EXPECT_NEAR(variance_reduced_number_density_b, standard_number_density_b, 1e-6);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(variance_reduced_bulk_velocity_b[i], standard_bulk_velocity_b[i], 1e-6);
  }
  EXPECT_NEAR(variance_reduced_temperature_b, standard_temperature_b,1e-6);
}

struct OneDimensionalMaxwellianVarianceReductionSetup {
  static constexpr int num_elems = 8;
  static constexpr int velocity_dims = 1;
  static constexpr double number_density = 1e18;
  static constexpr double temperature = 11604.5;
  const Species species{.charge = -constants::elementary_charge, .mass = constants::electron_mass, .name = "electron"};
  const mfem::Vector bulk_velocity{2.0e5, 0.0, 0.0};
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, 1e-3));
  Discretization charge_discretization{mesh.get(), 1, FETypes::HGRAD};
  ParticleContainer particles;
  ParticleOperations particle_operations{
    charge_discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory,
    {{"electron", species}},
    velocity_dims};

  OneDimensionalMaxwellianVarianceReductionSetup() {
    RandomNumberGenerator generator(2468);
    const SourceStateParameters source_state_parameters{
      .number_density = number_density,
      .bulk_velocity = bulk_velocity,
      .temperature = temperature,
    };
    particles = loadParticles(
      ConstantSourceParameters(species, source_state_parameters, 4000),
      generator,
      mesh,
      velocity_dims);
  }

  void useUniformMaxwellian(
    const double reference_number_density,
    const VarianceReductionParameters::CellLimiting cell_limiting = VarianceReductionParameters::CellLimiting::None,
    const double reference_temperature_ratio = 1.0
  ) {
    VarianceReductionParameters variance_reduction_parameters;
    variance_reduction_parameters.strategy = VarianceReductionParameters::Strategy::UniformMaxwellian;
    variance_reduction_parameters.reference_number_density = reference_number_density;
    variance_reduction_parameters.reference_bulk_velocity = {bulk_velocity(0), bulk_velocity(1), bulk_velocity(2)};
    variance_reduction_parameters.reference_temperature = reference_temperature_ratio * temperature;
    variance_reduction_parameters.cell_limiting = cell_limiting;
    particle_operations.setVarianceReductionParameters(variance_reduction_parameters);
  }

  mfem::Vector uniformCharge(const double charge_density) {
    mfem::ConstantCoefficient charge_density_coefficient(charge_density);
    mfem::LinearForm charge(&charge_discretization.getFeSpace());
    charge.AddDomainIntegrator(new mfem::DomainLFIntegrator(charge_density_coefficient));
    charge.Assemble();
    return charge;
  }
};

TEST(ParticleOperations, VarianceReducedChargeWithUniformMaxwellianIsExactWithoutLowFidelityModel) {
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  setup.useUniformMaxwellian(setup.number_density);

  const IntegratedCharge variance_reduced_charge = setup.particle_operations.assembleVarianceReducedCharge(setup.particles);
  const mfem::Vector expected_charge = setup.uniformCharge(setup.species.charge * setup.number_density);

  for (int dof = 0; dof < expected_charge.Size(); ++dof) {
    EXPECT_NEAR(variance_reduced_charge.getIntegratedCharge()(dof), expected_charge(dof), 1e-12 * std::abs(expected_charge(dof)));
  }
}

TEST(ParticleOperations, VarianceReducedMomentsWithUniformMaxwellianIgnoreLowFidelityModel) {
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  setup.useUniformMaxwellian(setup.number_density);

  // A fluid that disagrees with the particles; the uniform Maxwellian must be used instead.
  Discretization dg_discretization(setup.mesh.get(), 0, FETypes::DG, 5);
  std::shared_ptr<DGEulerAssembly> operator_ptr = std::make_shared<DGEulerAssembly>(dg_discretization.getFeSpace(), setup.species);
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({operator_ptr});
  DGEulerOperations dg_euler_operations(setup.charge_discretization, dg_operators);
  std::vector<std::unique_ptr<SourceParameters>> list_of_parameters;
  list_of_parameters.push_back(std::make_unique<ConstantSourceParameters>(
    setup.species, 0.5 * setup.number_density, 2.0 * setup.temperature, mfem::Vector({0.0, 0.0, 0.0})));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, list_of_parameters);

  const ParticleMoments with_fluid = setup.particle_operations.getVarianceReducedParticleMoments(
    setup.particles, low_fidelity_state, dg_euler_operations, true);
  const ParticleMoments without_fluid = setup.particle_operations.getVarianceReducedParticleMoments(setup.particles, true);

  for (const ParticleMoments& moments : {with_fluid, without_fluid}) {
    for (int elem_id = 0; elem_id < setup.num_elems; ++elem_id) {
      EXPECT_NEAR(moments.number_density.at(setup.species)(elem_id), setup.number_density, 1e-12 * setup.number_density);
      EXPECT_NEAR(moments.bulk_velocity.at(setup.species)(0, elem_id), setup.bulk_velocity(0), 1e-9 * setup.bulk_velocity(0));
      EXPECT_NEAR(moments.temperature.at(setup.species)(elem_id), setup.temperature, 1e-9 * setup.temperature);
    }
  }
}

TEST(ParticleOperations, VarianceReducedBulkVelocityWithUniformMaxwellianWeightsReferenceVelocityByDensity) {
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  // f_LF = 2 f at every particle, so each particle's noise reducing factor is exactly -1
  const double reference_number_density = 2.0 * setup.number_density;
  setup.useUniformMaxwellian(reference_number_density);

  const ParticleMoments variance_reduced = setup.particle_operations.getVarianceReducedParticleMoments(setup.particles, true);
  const ParticleMoments pic = setup.particle_operations.getParticleMoments(setup.particles);

  for (int elem_id = 0; elem_id < setup.num_elems; ++elem_id) {
    const double pic_number_density = pic.number_density.at(setup.species)(elem_id);
    const double pic_velocity = pic.bulk_velocity.at(setup.species)(0, elem_id);
    const double expected_number_density = reference_number_density - pic_number_density;
    const double expected_velocity =
      (reference_number_density * setup.bulk_velocity(0) - pic_number_density * pic_velocity) / expected_number_density;

    EXPECT_NEAR(variance_reduced.number_density.at(setup.species)(elem_id), expected_number_density, 1e-10 * setup.number_density);
    EXPECT_NEAR(variance_reduced.bulk_velocity.at(setup.species)(0, elem_id), expected_velocity, 1e-8 * std::abs(expected_velocity));
  }
}

TEST(ParticleOperations, VarianceReducedBulkVelocityFallsBackToPICWhereVarianceReducedDensityIsNotPositive) {
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  // f_LF = 3 f at every particle, so each particle's factor is exactly -2 and
  // n_VR = n_ref - 2 n_PIC. With n_ref = n0, n_VR < 0 wherever n_PIC > n0 / 2, i.e. everywhere.
  setup.useUniformMaxwellian(setup.number_density);
  for (Particle& particle : setup.particles) {
    particle.particle_distribution_function_value /= 3.0;
  }

  const ParticleMoments variance_reduced = setup.particle_operations.getVarianceReducedParticleMoments(setup.particles, true);
  const ParticleMoments pic = setup.particle_operations.getParticleMoments(setup.particles);

  for (int elem_id = 0; elem_id < setup.num_elems; ++elem_id) {
    ASSERT_LT(variance_reduced.number_density.at(setup.species)(elem_id), 0.0);
    for (int vel_dim = 0; vel_dim < 3; ++vel_dim) {
      const double velocity = variance_reduced.bulk_velocity.at(setup.species)(vel_dim, elem_id);
      EXPECT_TRUE(std::isfinite(velocity));
      EXPECT_DOUBLE_EQ(velocity, pic.bulk_velocity.at(setup.species)(vel_dim, elem_id));
    }
  }
}

TEST(ParticleOperations, VarianceReducedChargeRespectsCellLimiting) {
  using CellLimiting = VarianceReductionParameters::CellLimiting;
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  // f_LF = 3 f at every particle, so each particle's noise reducing factor is exactly -2.
  // That fails both cell limiting rules: |-2| >= 1, and sum w^2 (-2)^2 = 4 sum w^2.
  const double reference_number_density = 3.0 * setup.number_density;
  const mfem::Vector pic_charge = setup.particle_operations.assembleCharge(setup.particles).getIntegratedCharge();
  const mfem::Vector control_variate_charge = setup.uniformCharge(setup.species.charge * reference_number_density);

  setup.useUniformMaxwellian(reference_number_density, CellLimiting::None);
  const mfem::Vector unlimited_charge = setup.particle_operations.assembleVarianceReducedCharge(setup.particles).getIntegratedCharge();

  for (const CellLimiting cell_limiting : {CellLimiting::EveryParticle, CellLimiting::Variance}) {
    setup.useUniformMaxwellian(reference_number_density, cell_limiting);
    const mfem::Vector limited_charge = setup.particle_operations.assembleVarianceReducedCharge(setup.particles).getIntegratedCharge();

    for (int dof = 0; dof < pic_charge.Size(); ++dof) {
      const double scale = std::abs(control_variate_charge(dof));
      EXPECT_NEAR(limited_charge(dof), pic_charge(dof), 1e-12 * scale);
    }
  }
  for (int dof = 0; dof < pic_charge.Size(); ++dof) {
    const double scale = std::abs(control_variate_charge(dof));
    EXPECT_NEAR(unlimited_charge(dof), -2.0 * pic_charge(dof) + control_variate_charge(dof), 1e-12 * scale);
  }
}

TEST(ParticleOperations, VarianceReducedMomentsRespectCellLimiting) {
  using CellLimiting = VarianceReductionParameters::CellLimiting;
  OneDimensionalMaxwellianVarianceReductionSetup setup;

  for (const CellLimiting cell_limiting : {CellLimiting::EveryParticle, CellLimiting::Variance}) {
    setup.useUniformMaxwellian(3.0 * setup.number_density, cell_limiting);

    const ParticleMoments variance_reduced = setup.particle_operations.getVarianceReducedParticleMoments(setup.particles, true);
    const ParticleMoments pic = setup.particle_operations.getParticleMoments(setup.particles);

    for (int elem_id = 0; elem_id < setup.num_elems; ++elem_id) {
      EXPECT_DOUBLE_EQ(variance_reduced.number_density.at(setup.species)(elem_id), pic.number_density.at(setup.species)(elem_id));
      EXPECT_DOUBLE_EQ(variance_reduced.bulk_velocity.at(setup.species)(0, elem_id), pic.bulk_velocity.at(setup.species)(0, elem_id));
      EXPECT_DOUBLE_EQ(variance_reduced.temperature.at(setup.species)(elem_id), pic.temperature.at(setup.species)(elem_id));
    }
  }
}

TEST(ParticleOperations, VarianceCellLimitingKeepsVarianceReductionWhereEveryParticleRuleDoesNot) {
  using CellLimiting = VarianceReductionParameters::CellLimiting;
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  // A control variate 20% too hot: the mean of (1 - f_LF / f)^2 is about 0.02, so variance reduction
  // works well on average, but f_LF / f >= 2 beyond about 3.1 thermal speeds (~0.2% of particles).
  // With ~500 particles per element most elements contain such a particle.
  constexpr double reference_temperature_ratio = 1.2;
  const mfem::Vector pic_number_density = setup.particle_operations.getParticleMoments(setup.particles).number_density.at(setup.species);

  auto elements_using_variance_reduction = [&](const CellLimiting cell_limiting) {
    setup.useUniformMaxwellian(setup.number_density, cell_limiting, reference_temperature_ratio);
    const mfem::Vector number_density =
      setup.particle_operations.getVarianceReducedParticleMoments(setup.particles, true).number_density.at(setup.species);
    int count = 0;
    for (int elem_id = 0; elem_id < setup.num_elems; ++elem_id) {
      // A VR element differs from PIC; a fallen-back element equals it exactly
      if (number_density(elem_id) != pic_number_density(elem_id)) ++count;
    }
    return count;
  };

  const int every_particle = elements_using_variance_reduction(CellLimiting::EveryParticle);
  const int variance = elements_using_variance_reduction(CellLimiting::Variance);

  EXPECT_EQ(variance, setup.num_elems);
  EXPECT_LT(every_particle, setup.num_elems);
}

TEST(ParticleOperations, UniformMaxwellianOnlyAppliesToListedSpecies) {
  OneDimensionalMaxwellianVarianceReductionSetup setup;
  VarianceReductionParameters variance_reduction_parameters;
  variance_reduction_parameters.strategy = VarianceReductionParameters::Strategy::UniformMaxwellian;
  variance_reduction_parameters.reference_species = {"some_other_species"};
  variance_reduction_parameters.reference_number_density = 3.0 * setup.number_density;
  variance_reduction_parameters.reference_temperature = setup.temperature;
  setup.particle_operations.setVarianceReductionParameters(variance_reduction_parameters);

  const mfem::Vector pic_charge = setup.particle_operations.assembleCharge(setup.particles).getIntegratedCharge();
  const mfem::Vector variance_reduced_charge = setup.particle_operations.assembleVarianceReducedCharge(setup.particles).getIntegratedCharge();

  for (int dof = 0; dof < pic_charge.Size(); ++dof) {
    EXPECT_DOUBLE_EQ(variance_reduced_charge(dof), pic_charge(dof));
  }
}

double measureMomentum(const ParticleContainer& particles) {
  double total = 0.0;
  for (const Particle& particle : particles) total += particle.weight * particle.species.mass * particle.velocity(0);
  return total;
}

double measureKineticEnergy(const ParticleContainer& particles) {
  double total = 0.0;
  for (const Particle& particle : particles) total += 0.5 * particle.weight * particle.species.mass * (particle.velocity * particle.velocity);
  return total;
}

struct BGKVariant {
  const char* name;
};

void measureBGKEnergyChange(double lf_temperature_ratio, double lf_drift_in_thermal_speeds, int particles_per_cell) {
  Species species{.charge = constants::elementary_charge, .mass = constants::proton_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e20;
  constexpr double temperature = 1.16e5;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const mfem::Vector bulk_velocity({-3.0e4, 0.0, 0.0});
  mfem::Vector lf_bulk_velocity = bulk_velocity;
  lf_bulk_velocity(0) += lf_drift_in_thermal_speeds * thermal_speed;

  constexpr int num_elems = 40;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
  Discretization dg_discretization(mesh.get(), 0, FETypes::DG, 5);
  Discretization charge_discretization(mesh.get(), 1, FETypes::HGRAD);
  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({std::make_shared<DGEulerAssembly>(finite_element_space, species)});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);
  std::vector<std::unique_ptr<SourceParameters>> lf_parameters;
  lf_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, lf_temperature_ratio * temperature, lf_bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, lf_parameters);

  ParticleOperations particle_operations(charge_discretization, empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory, species_map);
  VarianceReductionParameters vr_parameters;
  vr_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  particle_operations.setVarianceReductionParameters(vr_parameters);

  const BGKVariant variants[] = {
    {"BGK"},
  };
  constexpr int num_trials = 2000;
  for (const BGKVariant& variant : variants) {
    BGKRelaxationOperations relaxer(1.0, species);
    double sum_before = 0.0, sum_difference = 0.0, sum_difference_squares = 0.0;
    double sum_momentum_before = 0.0, sum_momentum_difference = 0.0, sum_momentum_difference_squares = 0.0;
    for (int trial = 0; trial < num_trials; ++trial) {
      RandomNumberGenerator generator(trial);
      ParticleContainer particles = loadParticles(
        ConstantSourceParameters(species, SourceStateParameters{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature}, particles_per_cell * num_elems),
        generator, mesh, 3);
      const double energy_before = measureKineticEnergy(particles);
      const double momentum_before = measureMomentum(particles);
      relaxer.performCollisions(1.0e15, generator, particles, particle_operations, low_fidelity_state, dg_euler_operations);
      const double difference = measureKineticEnergy(particles) - energy_before;
      sum_before += energy_before; sum_difference += difference; sum_difference_squares += difference * difference;
      const double momentum_difference = measureMomentum(particles) - momentum_before;
      sum_momentum_before += std::abs(momentum_before);
      sum_momentum_difference += momentum_difference;
      sum_momentum_difference_squares += momentum_difference * momentum_difference;
    }
    const double mean_before = sum_before / num_trials;
    const double mean_difference = sum_difference / num_trials;
    const double standard_error = std::sqrt((sum_difference_squares / num_trials - mean_difference * mean_difference) / num_trials);
    const double mean_momentum_before = sum_momentum_before / num_trials;
    const double mean_momentum_difference = sum_momentum_difference / num_trials;
    const double momentum_standard_error = std::sqrt(
      (sum_momentum_difference_squares / num_trials - mean_momentum_difference * mean_momentum_difference) / num_trials);
    std::printf("MEASURE LF T/T=%.2f drift=%+.2f vth N/cell=%3d %s: <dE>/<E> = %+.3f%% +- %.3f%%  <dPx>/<|Px|> = %+.3f%% +- %.3f%%\n",
      lf_temperature_ratio, lf_drift_in_thermal_speeds, particles_per_cell, variant.name,
      100.0 * mean_difference / mean_before, 100.0 * standard_error / mean_before,
      100.0 * mean_momentum_difference / mean_momentum_before, 100.0 * momentum_standard_error / mean_momentum_before);
  }
}

TEST(ParticleOperations, DISABLED_MeasureBGKEnergyChangeWithVRMoments) {
  measureBGKEnergyChange(1.0, 0.0, 40);
  measureBGKEnergyChange(0.9, 0.0, 40);
  measureBGKEnergyChange(1.0, 0.2, 40);
  measureBGKEnergyChange(0.9, 0.2, 40);
  measureBGKEnergyChange(0.9, 0.2, 10);
  measureBGKEnergyChange(0.9, 0.2, 160);
}

struct RunningMean {
  double sum = 0.0;
  double sum_squares = 0.0;
  long count = 0;
  void add(double value) { sum += value; sum_squares += value * value; ++count; }
  double mean() const { return sum / count; }
  double standardError() const { return std::sqrt((sum_squares / count - mean() * mean()) / count); }
};

TEST(ParticleOperations, DISABLED_DiagnoseVRTemperatureBias) {
  Species species{.charge = constants::elementary_charge, .mass = constants::proton_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e20;
  constexpr double temperature = 1.16e5;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const mfem::Vector bulk_velocity({-3.0e4, 0.0, 0.0});
  const double three_k_over_m = 3.0 * constants::boltzmann_constant / species.mass;
  const double true_bulk_speed_squared = bulk_velocity * bulk_velocity;
  const double true_specific_energy = true_bulk_speed_squared + three_k_over_m * temperature;

  for (const auto& [lf_temperature_ratio, lf_drift, particles_per_cell] :
       {std::tuple{1.0, 0.0, 40}, std::tuple{0.9, 0.0, 40}, std::tuple{1.0, 0.2, 40},
        std::tuple{0.9, 0.2, 40}, std::tuple{0.9, 0.2, 10}, std::tuple{0.9, 0.2, 160}}) {
    mfem::Vector lf_bulk_velocity = bulk_velocity;
    lf_bulk_velocity(0) += lf_drift * thermal_speed;
    constexpr int num_elems = 40;
    std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
    Discretization dg_discretization(mesh.get(), 0, FETypes::DG, 5);
    Discretization charge_discretization(mesh.get(), 1, FETypes::HGRAD);
    mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
    std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({std::make_shared<DGEulerAssembly>(finite_element_space, species)});
    DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);
    std::vector<std::unique_ptr<SourceParameters>> lf_parameters;
    lf_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, lf_temperature_ratio * temperature, lf_bulk_velocity));
    LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, lf_parameters);
    ParticleOperations particle_operations(charge_discretization, empty_particle_boundary_factory_list,
      default_reflecting_particle_boundary_factory, species_map);
    VarianceReductionParameters vr_parameters;
    vr_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
    particle_operations.setVarianceReductionParameters(vr_parameters);

    RunningMean energy, weight, sampling, total_plugin, total_variance_corrected, total_energy_conserving, total_raw;
    RunningMean momentum, momentum_sampling, momentum_vr, momentum_corrected, energy_both_corrected;
    RunningMean bulk_weighted, bulk_unweighted, thermal_plugin_weighted, thermal_plugin_unweighted;
    RunningMean thermal_variance_corrected_weighted, thermal_variance_corrected_unweighted, vr_active;
    double max_raw_identity_violation = 0.0;

    constexpr int num_trials = 2000;
    for (int trial = 0; trial < num_trials; ++trial) {
      RandomNumberGenerator generator(trial);
      ParticleContainer particles = loadParticles(
        ConstantSourceParameters(species, SourceStateParameters{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature}, particles_per_cell * num_elems),
        generator, mesh, 3);
      std::vector<double> cell_weight(num_elems, 0.0), cell_energy(num_elems, 0.0), cell_momentum(num_elems, 0.0);
      for (const Particle& particle : particles) {
        cell_weight[particle.element] += particle.weight;
        cell_energy[particle.element] += particle.weight * (particle.velocity * particle.velocity);
        cell_momentum[particle.element] += particle.weight * particle.velocity(0);
      }

      const ParticleMoments vr = particle_operations.getVarianceReducedParticleMoments(particles, low_fidelity_state, dg_euler_operations, true);
      const ParticleMoments raw = particle_operations.getParticleMoments(particles);

      for (int e = 0; e < num_elems; ++e) {
        const double w = cell_weight[e];
        const double s = cell_energy[e];
        double vr_bulk_speed_squared = 0.0, raw_bulk_speed_squared = 0.0;
        for (int d = 0; d < 3; ++d) {
          vr_bulk_speed_squared += std::pow(vr.bulk_velocity.at(species)(d, e), 2);
          raw_bulk_speed_squared += std::pow(raw.bulk_velocity.at(species)(d, e), 2);
        }
        const double plugin_temperature = vr.temperature.at(species)(e) / vr.temperature_bias_correction.at(species)(e);
        const double variance_corrected_temperature = vr.relaxation_temperature.at(species)(e);
        const double energy_conserving_temperature = vr.energy_conserving_temperature.at(species)(e);
        const double raw_plugin_temperature = raw.temperature.at(species)(e) / raw.temperature_bias_correction.at(species)(e);

        const double raw_change = w * (raw_bulk_speed_squared + three_k_over_m * raw_plugin_temperature) - s;
        max_raw_identity_violation = std::max(max_raw_identity_violation, std::abs(raw_change) / s);

        energy.add(s);
        weight.add(w);
        sampling.add(w * true_specific_energy - s);
        total_raw.add(raw_change);
        total_plugin.add(w * (vr_bulk_speed_squared + three_k_over_m * plugin_temperature) - s);
        total_variance_corrected.add(w * (vr_bulk_speed_squared + three_k_over_m * variance_corrected_temperature) - s);
        total_energy_conserving.add(w * (vr_bulk_speed_squared + three_k_over_m * energy_conserving_temperature) - s);
        double corrected_bulk_speed_squared = 0.0;
        for (int d = 0; d < 3; ++d)
          corrected_bulk_speed_squared += std::pow(vr.energy_conserving_bulk_velocity.at(species)(d, e), 2);
        energy_both_corrected.add(w * (corrected_bulk_speed_squared + three_k_over_m * energy_conserving_temperature) - s);
        const double p = cell_momentum[e];
        momentum.add(std::abs(p));
        momentum_sampling.add(w * bulk_velocity(0) - p);
        momentum_vr.add(w * vr.bulk_velocity.at(species)(0, e) - p);
        momentum_corrected.add(w * vr.energy_conserving_bulk_velocity.at(species)(0, e) - p);
        bulk_weighted.add(w * (vr_bulk_speed_squared - true_bulk_speed_squared));
        bulk_unweighted.add(vr_bulk_speed_squared - true_bulk_speed_squared);
        thermal_plugin_weighted.add(w * three_k_over_m * (plugin_temperature - temperature));
        thermal_plugin_unweighted.add(three_k_over_m * (plugin_temperature - temperature));
        thermal_variance_corrected_weighted.add(w * three_k_over_m * (variance_corrected_temperature - temperature));
        thermal_variance_corrected_unweighted.add(three_k_over_m * (variance_corrected_temperature - temperature));
        vr_active.add(vr.number_density.at(species)(e) != raw.number_density.at(species)(e) ? 1.0 : 0.0);
      }
    }

    const double e0 = energy.mean();
    const double w0 = weight.mean();
    auto pct = [&](const RunningMean& m) { return 100.0 * m.mean() / e0; };
    auto pct_se = [&](const RunningMean& m) { return 100.0 * m.standardError() / e0; };
    auto expected_part = [&](const RunningMean& unweighted) { return 100.0 * w0 * unweighted.mean() / e0; };
    std::printf("DIAG ===== LF T/T=%.2f drift=%+.2f vth, N/cell=%3d, VR active in %.0f%% of cells\n",
      lf_temperature_ratio, lf_drift, particles_per_cell, 100.0 * vr_active.mean());
    std::printf("DIAG   raw identity: max |W(u^2+3kT/m)-S|/S = %.1e, mean %+.4f%% +- %.4f%%\n",
      max_raw_identity_violation, pct(total_raw), pct_se(total_raw));
    std::printf("DIAG   sampling check <W eps0 - S>: %+.3f%% +- %.3f%%\n", pct(sampling), pct_se(sampling));
    std::printf("DIAG   TOTAL predicted dE/E: plug-in %+.3f%% +- %.3f%% | var-corrected %+.3f%% +- %.3f%% | energy-conserving %+.3f%% +- %.3f%%\n",
      pct(total_plugin), pct_se(total_plugin), pct(total_variance_corrected), pct_se(total_variance_corrected),
      pct(total_energy_conserving), pct_se(total_energy_conserving));
    std::printf("DIAG   TOTAL predicted dE/E with corrected u and T: %+.3f%% +- %.3f%%\n",
      pct(energy_both_corrected), pct_se(energy_both_corrected));
    const double p0 = momentum.mean();
    std::printf("DIAG   x-momentum dP/|P|: sampling check %+.3f%% +- %.3f%% | VR u %+.3f%% +- %.3f%% | corrected u %+.3f%% +- %.3f%%\n",
      100.0 * momentum_sampling.mean() / p0, 100.0 * momentum_sampling.standardError() / p0,
      100.0 * momentum_vr.mean() / p0, 100.0 * momentum_vr.standardError() / p0,
      100.0 * momentum_corrected.mean() / p0, 100.0 * momentum_corrected.standardError() / p0);
    std::printf("DIAG   bulk    W(|u_VR|^2-|u0|^2)     : %+.3f%% +- %.3f%%  = E[W]E[.] %+.3f%% + cov %+.3f%%\n",
      pct(bulk_weighted), pct_se(bulk_weighted), expected_part(bulk_unweighted), pct(bulk_weighted) - expected_part(bulk_unweighted));
    std::printf("DIAG   thermal W 3k/m(T_plugin-T0)    : %+.3f%% +- %.3f%%  = E[W]E[.] %+.3f%% + cov %+.3f%%\n",
      pct(thermal_plugin_weighted), pct_se(thermal_plugin_weighted), expected_part(thermal_plugin_unweighted),
      pct(thermal_plugin_weighted) - expected_part(thermal_plugin_unweighted));
    std::printf("DIAG   thermal W 3k/m(T_varcorr-T0)   : %+.3f%% +- %.3f%%  = E[W]E[.] %+.3f%% + cov %+.3f%%\n",
      pct(thermal_variance_corrected_weighted), pct_se(thermal_variance_corrected_weighted),
      expected_part(thermal_variance_corrected_unweighted),
      pct(thermal_variance_corrected_weighted) - expected_part(thermal_variance_corrected_unweighted));
  }
}


void measureRepeatedBGKEnergyDrift(double lf_temperature_ratio, double lf_drift_in_thermal_speeds, int particles_per_cell, bool use_vr_moments) {
  Species species{.charge = constants::elementary_charge, .mass = constants::proton_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e20;
  constexpr double temperature = 1.16e5;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const mfem::Vector bulk_velocity({-3.0e4, 0.0, 0.0});
  mfem::Vector lf_bulk_velocity = bulk_velocity;
  lf_bulk_velocity(0) += lf_drift_in_thermal_speeds * thermal_speed;

  constexpr int num_elems = 40;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems));
  Discretization dg_discretization(mesh.get(), 0, FETypes::DG, 5);
  Discretization charge_discretization(mesh.get(), 1, FETypes::HGRAD);
  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({std::make_shared<DGEulerAssembly>(finite_element_space, species)});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);
  std::vector<std::unique_ptr<SourceParameters>> lf_parameters;
  lf_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, lf_temperature_ratio * temperature, lf_bulk_velocity));
  LowFidelityState low_fidelity_state = buildEulerState(dg_discretization, lf_parameters);

  ParticleOperations particle_operations(charge_discretization, empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory, species_map);
  VarianceReductionParameters vr_parameters;
  vr_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  particle_operations.setVarianceReductionParameters(vr_parameters);

  constexpr double collision_frequency = 1.0;
  constexpr double dt = 0.1;
  constexpr int num_calls = 100;
  constexpr int num_trials = 200;
  const BGKRelaxationOperations relaxer(collision_frequency, species, 1.0, 0.0, 0, use_vr_moments);
  const std::vector<std::pair<int, int>> windows{{0, 10}, {10, 30}, {30, 60}, {60, 100}};
  std::vector<RunningMean> energy_change(windows.size()), vr_minus_raw_temperature(windows.size()), active_fraction(windows.size());
  RunningMean final_energy_change;
  for (int trial = 0; trial < num_trials; ++trial) {
    RandomNumberGenerator generator(trial);
    ParticleContainer particles = loadParticles(
      ConstantSourceParameters(species, SourceStateParameters{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature}, particles_per_cell * num_elems),
      generator, mesh, 3);
    const double initial_energy = measureKineticEnergy(particles);
    std::vector<double> window_change(windows.size(), 0.0), window_temperature(windows.size(), 0.0), window_active(windows.size(), 0.0);
    for (int call = 0; call < num_calls; ++call) {
      const ParticleMoments vr = particle_operations.getVarianceReducedParticleMoments(particles, low_fidelity_state, dg_euler_operations, true);
      const ParticleMoments raw = particle_operations.getParticleMoments(particles);
      double temperature_difference = 0.0, active = 0.0;
      for (int e = 0; e < num_elems; ++e) {
        active += vr.number_density.at(species)(e) != raw.number_density.at(species)(e);
        temperature_difference += vr.energy_conserving_temperature.at(species)(e) / raw.relaxation_temperature.at(species)(e) - 1.0;
      }
      const double energy_before = measureKineticEnergy(particles);
      relaxer.performCollisions(dt, generator, particles, particle_operations, low_fidelity_state, dg_euler_operations);
      const double change = (measureKineticEnergy(particles) - energy_before) / initial_energy;
      for (size_t w = 0; w < windows.size(); ++w)
        if (call >= windows[w].first and call < windows[w].second) {
          const double length = windows[w].second - windows[w].first;
          window_change[w] += change / length;
          window_temperature[w] += temperature_difference / num_elems / length;
          window_active[w] += active / num_elems / length;
        }
    }
    for (size_t w = 0; w < windows.size(); ++w) {
      energy_change[w].add(window_change[w]);
      vr_minus_raw_temperature[w].add(window_temperature[w]);
      active_fraction[w].add(window_active[w]);
    }
    final_energy_change.add(measureKineticEnergy(particles) / initial_energy - 1.0);
  }
  std::printf("DRIFT LF T/T=%.2f drift=%+.2f vth N/cell=%3d %s: cumulative dE/E after %d calls = %+.3f%% +- %.3f%%\n",
    lf_temperature_ratio, lf_drift_in_thermal_speeds, particles_per_cell, use_vr_moments ? "VR-moment BGK" : "raw BGK      ",
    num_calls, 100.0 * final_energy_change.mean(), 100.0 * final_energy_change.standardError());
  for (size_t w = 0; w < windows.size(); ++w)
    std::printf("    calls %3d-%3d: dE/E per call = %+.4f%% +- %.4f%%   <T_VR,ec / T_raw,ec - 1> = %+.3f%% +- %.3f%%   VR active %.0f%%\n",
      windows[w].first + 1, windows[w].second, 100.0 * energy_change[w].mean(), 100.0 * energy_change[w].standardError(),
      100.0 * vr_minus_raw_temperature[w].mean(), 100.0 * vr_minus_raw_temperature[w].standardError(), 100.0 * active_fraction[w].mean());
}

TEST(ParticleOperations, DISABLED_RepeatedBGKDriftRaw) { measureRepeatedBGKEnergyDrift(1.0, 0.0, 100, false); }
TEST(ParticleOperations, DISABLED_RepeatedBGKDriftVRMatched) { measureRepeatedBGKEnergyDrift(1.0, 0.0, 100, true); }
TEST(ParticleOperations, DISABLED_RepeatedBGKDriftVRCold5) { measureRepeatedBGKEnergyDrift(0.95, 0.0, 100, true); }
TEST(ParticleOperations, DISABLED_RepeatedBGKDriftVRCold10) { measureRepeatedBGKEnergyDrift(0.9, 0.0, 100, true); }
TEST(ParticleOperations, DISABLED_RepeatedBGKDriftVRCold10Drift) { measureRepeatedBGKEnergyDrift(0.9, 0.2, 100, true); }

void measureVRMomentBiasUnderRepeatedBGK(double lf_temperature_ratio, double lf_drift_in_thermal_speeds, int particles_per_cell, bool stamp_raw_density = false,
  double cells_per_call_at_thermal_speed = 0.0, double final_lf_temperature_ratio = -1.0) {
  Species species{.charge = constants::elementary_charge, .mass = constants::proton_mass};
  const std::unordered_map<std::string, Species> species_map{{"one", species}};
  constexpr double number_density = 1e20;
  constexpr double temperature = 1.16e5;
  const double thermal_speed = std::sqrt(constants::boltzmann_constant * temperature / species.mass);
  const mfem::Vector bulk_velocity({-3.0e4, 0.0, 0.0});
  mfem::Vector lf_bulk_velocity = bulk_velocity;
  lf_bulk_velocity(0) += lf_drift_in_thermal_speeds * thermal_speed;

  constexpr int num_elems = 40;
  mfem::Mesh line_mesh = mfem::Mesh::MakeCartesian1D(num_elems);
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(
    mfem::Mesh::MakePeriodic(line_mesh, line_mesh.CreatePeriodicVertexMapping({mfem::Vector({1.0})})));
  Discretization dg_discretization(mesh.get(), 0, FETypes::DG, 5);
  Discretization charge_discretization(mesh.get(), 1, FETypes::HGRAD);
  mfem::FiniteElementSpace finite_element_space = dg_discretization.getFeSpace();
  std::vector<std::shared_ptr<DGEulerAssembly>> dg_operators({std::make_shared<DGEulerAssembly>(finite_element_space, species)});
  DGEulerOperations dg_euler_operations(charge_discretization, dg_operators);
  constexpr int num_lf_states = 21;
  const double end_ratio = final_lf_temperature_ratio > 0.0 ? final_lf_temperature_ratio : lf_temperature_ratio;
  std::vector<LowFidelityState> lf_states;
  for (int k = 0; k < num_lf_states; ++k) {
    const double ratio = lf_temperature_ratio + (end_ratio - lf_temperature_ratio) * k / (num_lf_states - 1.0);
    std::vector<std::unique_ptr<SourceParameters>> lf_parameters;
    lf_parameters.push_back(std::make_unique<ConstantSourceParameters>(species, number_density, ratio * temperature, lf_bulk_velocity));
    lf_states.push_back(buildEulerState(dg_discretization, lf_parameters));
  }

  ParticleOperations particle_operations(charge_discretization, empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory, species_map);
  VarianceReductionParameters vr_parameters;
  vr_parameters.strategy = VarianceReductionParameters::Strategy::EulerFluid;
  particle_operations.setVarianceReductionParameters(vr_parameters);

  constexpr double dt = 0.1;
  constexpr int num_calls = 100;
  constexpr int num_trials = 100;
  const double move_dt = cells_per_call_at_thermal_speed / num_elems / thermal_speed;
  const double relaxation_probability = -std::expm1(-dt);
  const BGKRelaxationOperations relaxer(1.0, species, 1.0, 0.0, 0, true, false, stamp_raw_density);
  const double three_k_over_m = 3.0 * constants::boltzmann_constant / species.mass;
  const std::vector<std::pair<int, int>> windows{{0, 10}, {10, 30}, {30, 60}, {60, 100}};
  enum Quantity { kEnergyTarget, kEnergyUncorrected, kVelocity, kNumber, kCalibration, kStamped, kMeasured, kEnergyThermal, kEnergyBulk, kNumQuantities };
  std::vector<std::vector<RunningMean>> stats(windows.size(), std::vector<RunningMean>(kNumQuantities));
  for (int trial = 0; trial < num_trials; ++trial) {
    RandomNumberGenerator generator(trial);
    ParticleContainer particles = loadParticles(
      ConstantSourceParameters(species, SourceStateParameters{.number_density = number_density, .bulk_velocity = bulk_velocity, .temperature = temperature}, particles_per_cell * num_elems),
      generator, mesh, 3);
    std::vector<double> initial_pdf;
    for (const Particle& particle : particles) initial_pdf.push_back(particle.particle_distribution_function_value);
    std::vector<std::vector<double>> window_sums(windows.size(), std::vector<double>(kNumQuantities, 0.0));
    for (int call = 0; call < num_calls; ++call) {
      const LowFidelityState& low_fidelity_state = lf_states[(call * (num_lf_states - 1)) / (num_calls - 1)];
      const ParticleMoments vr = particle_operations.getVarianceReducedParticleMoments(particles, low_fidelity_state, dg_euler_operations, true);
      std::vector<double> mass(num_elems, 0.0), momentum(num_elems, 0.0), speed_squared(num_elems, 0.0), calibration(num_elems, 0.0);
      std::vector<double> raw_y_momentum(num_elems, 0.0), raw_z_momentum(num_elems, 0.0);
      double stamped = 0.0, count = 0.0;
      int index = 0;
      for (const Particle& particle : particles) {
        const int e = particle.element;
        mass[e] += particle.weight;
        momentum[e] += particle.weight * particle.velocity(0);
        raw_y_momentum[e] += particle.weight * particle.velocity(1);
        raw_z_momentum[e] += particle.weight * particle.velocity(2);
        speed_squared[e] += particle.weight * (particle.velocity * particle.velocity);
        const mfem::Vector position(particle.position.GetData(), 1);
        calibration[e] += particle.weight * dg_euler_operations.evaluateParticleDistributionFunction(low_fidelity_state, position, particle.velocity, e, 0)
          / particle.particle_distribution_function_value;
        stamped += particle.particle_distribution_function_value != initial_pdf[index++];
        count += 1.0;
      }
      double total_raw_energy = 0.0, target_bias = 0.0, uncorrected_bias = 0.0, velocity_bias = 0.0, number_bias = 0.0, calibration_bias = 0.0;
      double thermal_bias = 0.0, bulk_bias = 0.0;
      for (int e = 0; e < num_elems; ++e) {
        const double volume = mesh->GetElementVolume(e);
        const double raw_specific_energy = speed_squared[e] / mass[e];
        double target_bulk_squared = 0.0, vr_bulk_squared = 0.0;
        for (int d = 0; d < 3; ++d) {
          target_bulk_squared += std::pow(vr.energy_conserving_bulk_velocity.at(species)(d, e), 2);
          vr_bulk_squared += std::pow(vr.bulk_velocity.at(species)(d, e), 2);
        }
        const double target_specific_energy = target_bulk_squared + three_k_over_m * vr.energy_conserving_temperature.at(species)(e);
        const double uncorrected_specific_energy = vr_bulk_squared + three_k_over_m * vr.relaxation_temperature.at(species)(e);
        total_raw_energy += speed_squared[e];
        target_bias += mass[e] * (target_specific_energy - raw_specific_energy);
        const double raw_bulk_squared = std::pow(momentum[e] / mass[e], 2);
        const double raw_thermal = raw_specific_energy - raw_bulk_squared - std::pow(raw_y_momentum[e] / mass[e], 2) - std::pow(raw_z_momentum[e] / mass[e], 2);
        thermal_bias += mass[e] * (three_k_over_m * vr.energy_conserving_temperature.at(species)(e) - raw_thermal);
        bulk_bias += mass[e] * (target_bulk_squared - (raw_specific_energy - raw_thermal));
        uncorrected_bias += mass[e] * (uncorrected_specific_energy - raw_specific_energy);
        velocity_bias += (vr.energy_conserving_bulk_velocity.at(species)(0, e) - momentum[e] / mass[e]) / thermal_speed / num_elems;
        number_bias += (vr.number_density.at(species)(e) * volume / mass[e] - 1.0) / num_elems;
        calibration_bias += (calibration[e] / (number_density * volume) - 1.0) / num_elems;
      }
      const double energy_before = measureKineticEnergy(particles);
      relaxer.performCollisions(dt, generator, particles, particle_operations, low_fidelity_state, dg_euler_operations);
      const double measured = (measureKineticEnergy(particles) - energy_before) / energy_before;
      if (move_dt > 0.0) particles = particle_operations.move(move_dt, particles);
      const double values[kNumQuantities] = {relaxation_probability * target_bias / total_raw_energy, relaxation_probability * uncorrected_bias / total_raw_energy,
        velocity_bias, number_bias, calibration_bias, stamped / count, measured,
        relaxation_probability * thermal_bias / total_raw_energy, relaxation_probability * bulk_bias / total_raw_energy};
      for (size_t w = 0; w < windows.size(); ++w)
        if (call >= windows[w].first and call < windows[w].second)
          for (int q = 0; q < kNumQuantities; ++q) window_sums[w][q] += values[q] / (windows[w].second - windows[w].first);
    }
    for (size_t w = 0; w < windows.size(); ++w)
      for (int q = 0; q < kNumQuantities; ++q) stats[w][q].add(window_sums[w][q]);
  }
  std::printf("BIAS LF T/T=%.2f->%.2f drift=%+.2f vth N/cell=%3d stamp=%s streaming=%.2f cells/call\n", lf_temperature_ratio, end_ratio, lf_drift_in_thermal_speeds, particles_per_cell, stamp_raw_density ? "raw n" : "VR n", cells_per_call_at_thermal_speed);
  for (size_t w = 0; w < windows.size(); ++w) {
    const auto& s = stats[w];
    std::printf("  calls %3d-%3d stamped %3.0f%% | dE/E per call: measured %+.4f%% +- %.4f%%, predicted from target %+.4f%% +- %.4f%%, from uncorrected VR %+.4f%% +- %.4f%%"
      " [thermal %+.4f%% +- %.4f%%, bulk %+.4f%% +- %.4f%%]"
      " | u_target-u_raw %+.4f +- %.4f vth | n_VR/n_raw-1 %+.3f%% +- %.3f%% | sum w f_LF/f_HF / N_LF - 1 %+.3f%% +- %.3f%%\n",
      windows[w].first + 1, windows[w].second, 100.0 * s[kStamped].mean(),
      100.0 * s[kMeasured].mean(), 100.0 * s[kMeasured].standardError(),
      100.0 * s[kEnergyTarget].mean(), 100.0 * s[kEnergyTarget].standardError(),
      100.0 * s[kEnergyUncorrected].mean(), 100.0 * s[kEnergyUncorrected].standardError(),
      100.0 * s[kEnergyThermal].mean(), 100.0 * s[kEnergyThermal].standardError(),
      100.0 * s[kEnergyBulk].mean(), 100.0 * s[kEnergyBulk].standardError(),
      s[kVelocity].mean(), s[kVelocity].standardError(),
      100.0 * s[kNumber].mean(), 100.0 * s[kNumber].standardError(),
      100.0 * s[kCalibration].mean(), 100.0 * s[kCalibration].standardError());
  }
}

TEST(ParticleOperations, DISABLED_VRMomentBiasMatched) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 100); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold10Drift) { measureVRMomentBiasUnderRepeatedBGK(0.9, 0.2, 100); }
TEST(ParticleOperations, DISABLED_VRMomentBiasMatchedRawStamp) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, true); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5RawStamp) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 100, true); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold10DriftRawStamp) { measureVRMomentBiasUnderRepeatedBGK(0.9, 0.2, 100, true); }
TEST(ParticleOperations, DISABLED_VRMomentBiasMatchedRawStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, true, 0.5); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5RawStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 100, true, 0.5); }
TEST(ParticleOperations, DISABLED_VRMomentBiasMatchedVRStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, false, 0.5); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5VRStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 100, false, 0.5); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5RawStampStreaming25) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 25, true, 0.5); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCold5RawStampStatic25) { measureVRMomentBiasUnderRepeatedBGK(0.95, 0.0, 25, true, 0.0); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCoolingRawStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, true, 0.5, 0.95); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCoolingRawStampStatic) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, true, 0.0, 0.95); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCooling10RawStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, true, 0.5, 0.90); }
TEST(ParticleOperations, DISABLED_VRMomentBiasCoolingVRStampStreaming) { measureVRMomentBiasUnderRepeatedBGK(1.0, 0.0, 100, false, 0.5, 0.95); }

} // namespace
