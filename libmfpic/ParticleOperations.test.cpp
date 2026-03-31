#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
#include <libmfpic/LoadUniformMaxwellianParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>
#include <libmfpic/ReflectingParticleBoundary.hpp>
#include <libmfpic/Species.hpp>

#include <gtest/gtest.h>
#include <mfem/mesh/element.hpp>

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

constexpr Species default_species{.charge = 1.0, .mass = 1.0};
const std::vector<std::shared_ptr<ParticleBoundaryFactory>> empty_particle_boundary_factory_list;
const std::shared_ptr<ParticleBoundaryFactory> default_reflecting_particle_boundary_factory
  = std::make_shared<ReflectingParticleBoundaryFactory>();

struct MomentsInCell {
  double number_density;
  std::vector<double> bulk_velocity;
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
      default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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
    default_reflecting_particle_boundary_factory
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

TEST(ParticleOperations, ParticleMomentsCorrectForMaxwellian) {
  const int num_elems = 1;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian1D(num_elems, .234));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  //const mfem::Vector nominal_bulk_velocity({300.0, 600.0, 1000.0});
  const mfem::Vector nominal_bulk_velocity({300.0,0.0,0.0});
  constexpr double temperature = 11600.0;
  constexpr double number_density = 1.0e18;
  constexpr int num_particles = 200000;
  std::mt19937 generator;

  ParticleContainer particles = loadUniformMaxwellianParticles(
    default_species,
    nominal_bulk_velocity,
    temperature,
    number_density,
    num_particles,
    generator,
    mesh
  );

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory
  );

  particle_operations.computeParticleMoments(particles);
  std::vector<double> computed_bulk_velocity = particle_operations.getBulkVelocity(0);
  double computed_number_density = particle_operations.getNumberDensity(0);
  double computed_temperature = particle_operations.getTemperature(0);

  EXPECT_NEAR(computed_number_density, number_density, 1e-10);

  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(computed_bulk_velocity[i], nominal_bulk_velocity[i], 1e-10);
  }

  EXPECT_NEAR(computed_temperature, temperature, 1e-3*temperature);
}

TEST(ParticleOperations, ParticleMomentsCorrectForKnownParticles) {

  const int num_elems = 2;
  const double lx = 1.0, ly = 2.0, lz = 3.0;
  std::shared_ptr<mfem::Mesh> mesh = std::make_shared<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(num_elems, num_elems, num_elems, mfem::Element::HEXAHEDRON, lx, ly, lz));
  constexpr int order = 1;
  Discretization discretization(mesh.get(),order);

  Species electron {.mass = constants::electron_mass};

  ParticleContainer particles;
  particles.addParticle( { .position = mfem::Vector({0.3411759316240717, 0.053821018802222675, 0.33053980915891706}),
               .velocity = mfem::Vector({1.5320330796287964, -0.6599694137918207, -0.31179485646991756}),
               .element = 0,
               .species = electron,
               .weight = 1.8149956155715332 });

  particles.addParticle( { .position = mfem::Vector({0.09218590534933485, 0.17590590108503035, 1.2181417599836606}),
                 .velocity = mfem::Vector({0.337769126558826, -2.2074710981998042, 0.8279214415587369}),
                 .element = 0,
                 .species = electron,
                 .weight = 0.7149484494297222 });
  
  particles.addParticle( { .position = mfem::Vector({0.4616724990135282, 0.27657439779710624, 1.2296318423895032}),
                 .velocity = mfem::Vector({1.541630394690618, 1.126806793265028, 0.7547696443122508}),
                 .element = 0,
                 .species = electron,
                 .weight = 1.648943384621905 });
  
  particles.addParticle( { .position = mfem::Vector({0.44494634655559295, 0.5129704552295319, 0.3674469016031947}),
                 .velocity = mfem::Vector({-0.14597789311522394, 1.2819022270597127, 1.0740306219719435}),
                 .element = 0,
                 .species = electron,
                 .weight = 1.0705371821926435 });
  
  particles.addParticle( { .position = mfem::Vector({0.9593002458867717, 0.625534005735464, 1.3756838588214035}),
                 .velocity = mfem::Vector({-0.6737591499870575, -0.6390601004322052, -0.061361327620372906}),
                 .element = 1,
                 .species = electron,
                 .weight = 1.4391659997195472 });
  
  particles.addParticle( { .position = mfem::Vector({0.9323451255468874, 0.21814287324998793, 1.2991911461073644}),
                 .velocity = mfem::Vector({-0.39278492256994324, 2.2899099473145785, -0.718181147880596}),
                 .element = 1,
                 .species = electron,
                 .weight = 0.5321658136231086 });
  
  particles.addParticle( { .position = mfem::Vector({0.8653759681856269, 0.2778652902989278, 1.1955653297751967}),
                 .velocity = mfem::Vector({0.03260774315697052, 0.028049895585638977, 0.028272122739737816}),
                 .element = 1,
                 .species = electron,
                 .weight = 2.0441238517620315 });
  
  particles.addParticle( { .position = mfem::Vector({0.9326108564218683, 0.29943789563745726, 0.7905631260454454}),
                 .velocity = mfem::Vector({0.05534586195270899, -0.48156285818994926, -0.5834075004641189}),
                 .element = 1,
                 .species = electron,
                 .weight = 0.9247446233849221 });
  
  particles.addParticle( { .position = mfem::Vector({0.41755117365428474, 1.61649125141019, 0.3991258680169122}),
                 .velocity = mfem::Vector({-0.2116374333746801, 0.36396383157911427, 0.9529644919745344}),
                 .element = 2,
                 .species = electron,
                 .weight = 1.2454537779347308 });
  
  particles.addParticle( { .position = mfem::Vector({0.40551105565904155, 1.4994867503030171, 1.1382154815058017}),
                 .velocity = mfem::Vector({1.5195241292412778, 1.7039094489490383, -0.24885870743094268}),
                 .element = 2,
                 .species = electron,
                 .weight = 0.7005794409282443 });
  
  particles.addParticle( { .position = mfem::Vector({0.2830445431633594, 1.437440361557259, 0.5942316666261351}),
                 .velocity = mfem::Vector({-0.4997485911335096, 0.09959750192201784, 0.12834321228831272}),
                 .element = 2,
                 .species = electron,
                 .weight = 1.1601924317153436 });
  
  particles.addParticle( { .position = mfem::Vector({0.011117643959939416, 1.4693507881679626, 0.9353376032371585}),
                 .velocity = mfem::Vector({-0.7342218924448597, -0.620475288234763, 0.8132737204208481}),
                 .element = 2,
                 .species = electron,
                 .weight = 0.8031380706370044 });
  
  particles.addParticle( { .position = mfem::Vector({0.9204110810881805, 1.9504758622500842, 0.4783118681015746}),
                 .velocity = mfem::Vector({0.826273507011688, 0.8503222896230244, -0.5157675930132671}),
                 .element = 3,
                 .species = electron,
                 .weight = 1.9306038956608929 });
  
  particles.addParticle( { .position = mfem::Vector({0.9488841444631287, 1.337529050610342, 1.2181681657703647}),
                 .velocity = mfem::Vector({1.6581133183034882, -0.2972625956392543, -1.3833771199483726}),
                 .element = 3,
                 .species = electron,
                 .weight = 0.5427140001090897 });
  
  particles.addParticle( { .position = mfem::Vector({0.899421799029557, 1.655285177011079, 0.3430551821634081}),
                 .velocity = mfem::Vector({-0.28120450483617765, 0.3600205099410753, -0.2343920149270695}),
                 .element = 3,
                 .species = electron,
                 .weight = 0.8591638156204775 });
  
  particles.addParticle( { .position = mfem::Vector({0.5688372323711341, 1.4243711389270113, 0.2273081295289603}),
                 .velocity = mfem::Vector({2.265520599867209, 0.8553866507799404, 1.7312794387228285}),
                 .element = 3,
                 .species = electron,
                 .weight = 0.5560793173287019 });
  
  particles.addParticle( { .position = mfem::Vector({0.11491071892245852, 0.37065116516756513, 1.9928326312217213}),
                 .velocity = mfem::Vector({0.31450349225488833, 0.07612150815035533, 0.25463292349897276}),
                 .element = 4,
                 .species = electron,
                 .weight = 0.8799360120615429 });
  
  particles.addParticle( { .position = mfem::Vector({0.3063208359090755, 0.544715400191636, 2.624464657336363}),
                 .velocity = mfem::Vector({2.0012311818026616, -0.3052079236753877, -0.5397628579464583}),
                 .element = 4,
                 .species = electron,
                 .weight = 1.4985393263659237 });
  
  particles.addParticle( { .position = mfem::Vector({0.1474771427130201, 0.9605537283750304, 2.3152241379001586}),
                 .velocity = mfem::Vector({1.41363079410976, -0.7057006440785899, 1.719889475677248}),
                 .element = 4,
                 .species = electron,
                 .weight = 1.13545033416701 });
  
  particles.addParticle( { .position = mfem::Vector({0.33971207420475874, 0.6882592484926856, 1.7284466891240025}),
                 .velocity = mfem::Vector({-0.1941978062119582, 0.0738395180442688, 0.7336282811408714}),
                 .element = 4,
                 .species = electron,
                 .weight = 0.8670632210358106 });
  
  particles.addParticle( { .position = mfem::Vector({0.5401938581949015, 0.8547157042203279, 1.717016096057045}),
                 .velocity = mfem::Vector({0.4765285113330109, -0.8532503175609265, -0.3582919723298022}),
                 .element = 5,
                 .species = electron,
                 .weight = 0.9250509143877632 });
  
  particles.addParticle( { .position = mfem::Vector({0.5896205388645837, 0.7846734732992243, 2.8431360849181417}),
                 .velocity = mfem::Vector({-0.0892000537049489, 0.385104174182408, 0.46116198661092983}),
                 .element = 5,
                 .species = electron,
                 .weight = 1.0513840562313006 });
  
  particles.addParticle( { .position = mfem::Vector({0.7303671630361488, 0.038483682202220315, 1.972957890318807}),
                 .velocity = mfem::Vector({0.054223542827762254, -1.9313076540480802, -0.20659269834007993}),
                 .element = 5,
                 .species = electron,
                 .weight = 0.7557422246863237 });
  
  particles.addParticle( { .position = mfem::Vector({0.5076466248462064, 0.3491736851311015, 2.4889650283304525}),
                 .velocity = mfem::Vector({-0.7755394099071095, -0.11134755206603067, 0.9913261626706856}),
                 .element = 5,
                 .species = electron,
                 .weight = 0.5971129520549291 });
  
  particles.addParticle( { .position = mfem::Vector({0.20939189496146027, 1.143017705710068, 2.479004695683661}),
                 .velocity = mfem::Vector({0.06375631404055183, 0.21318941825068694, -0.7507960478659627}),
                 .element = 6,
                 .species = electron,
                 .weight = 2.0781151283277675 });
  
  particles.addParticle( { .position = mfem::Vector({0.22889296342716692, 1.3938301513377254, 2.707009706277322}),
                 .velocity = mfem::Vector({-1.7253795282710642, -0.799603946712547, 1.0796056694391833}),
                 .element = 6,
                 .species = electron,
                 .weight = 0.9093364253485074 });
  
  particles.addParticle( { .position = mfem::Vector({0.22950755520689126, 1.6915683373154118, 2.672983615639859}),
                 .velocity = mfem::Vector({-0.3286666075267842, -0.32670660696567755, 1.5168019034425484}),
                 .element = 6,
                 .species = electron,
                 .weight = 2.8037905768120863 });
  
  particles.addParticle( { .position = mfem::Vector({0.30744062912078085, 1.9613908947704324, 2.3263170711324914}),
                 .velocity = mfem::Vector({0.4935377921925973, 0.7023083586133827, 0.8494486593868219}),
                 .element = 6,
                 .species = electron,
                 .weight = 0.7820965290029297 });
  
  particles.addParticle( { .position = mfem::Vector({0.5399070095380727, 1.1369895530226886, 1.9040233862539042}),
                 .velocity = mfem::Vector({1.130119646001526, 2.62894656687513, 0.4528764326586353}),
                 .element = 7,
                 .species = electron,
                 .weight = 2.3658743609885664 });
  
  particles.addParticle( { .position = mfem::Vector({0.8886869751678957, 1.229526008796011, 2.1326903722053556}),
                 .velocity = mfem::Vector({0.23403931476303633, 0.6986958268525457, 0.7286285854727402}),
                 .element = 7,
                 .species = electron,
                 .weight = 1.4965009591409915 });
  
  particles.addParticle( { .position = mfem::Vector({0.7290682674871196, 1.3096859921393007, 1.7655163540098673}),
                 .velocity = mfem::Vector({0.13491633484718182, -0.3738158660469447, 0.36175920966930913}),
                 .element = 7,
                 .species = electron,
                 .weight = 1.2458860130475344 });
  
  particles.addParticle( { .position = mfem::Vector({0.6340153936882715, 1.8523543812200487, 2.9911084120705307}),
                 .velocity = mfem::Vector({-1.371789209545305, -1.8824265432952796, 1.3609386407215076}),
                 .element = 7,
                 .species = electron,
                 .weight = 0.31090296242633314 });

  ParticleOperations particle_operations(
    discretization,
    empty_particle_boundary_factory_list,
    default_reflecting_particle_boundary_factory
  );

  particle_operations.computeParticleMoments(particles);

  auto check_moments = [&] (MomentsInCell exact, int cell_id) {

    MomentsInCell computed {.number_density = particle_operations.getNumberDensity(cell_id),
                            .bulk_velocity = particle_operations.getBulkVelocity(cell_id),
                            .temperature = particle_operations.getTemperature(cell_id) };

    EXPECT_NEAR(exact.number_density, computed.number_density, 1e-12);

    for (int i = 0; i < 3; i++) {
      EXPECT_NEAR(exact.bulk_velocity[i], computed.bulk_velocity[i], 1e-12);
    }

    EXPECT_NEAR(exact.temperature, computed.temperature, 1e-12);
  };

  MomentsInCell m0 {.number_density = 6.999232842421073,
                  .bulk_velocity = {1.0301904862978897,0.08654186050955474,0.4610747178478363},
                  .temperature = 5.2474391444070144e-08};

  MomentsInCell m1 {.number_density = 6.586933717986145,
                    .bulk_velocity = {-0.2147368229633683,-0.018033070767927425,-0.19274761437963095},
                    .temperature = 2.011379915422906e-08};

  MomentsInCell m2 {.number_density = 5.21248496162043,
                    .bulk_velocity = {-0.09426683972107544,0.3233899658465975,0.4641680775506456},
                    .temperature = 3.012395980852108e-08};

  MomentsInCell m3 {.number_density = 5.1847480382922155,
                    .bulk_velocity = {0.9034952905212801,0.5825513231401855,-0.2533516889623397},
                    .temperature = 3.5162259049005926e-08};

  MomentsInCell m4 {.number_density = 5.841318524840383,
                    .bulk_velocity = {1.075645812495728,-0.257395886909774,0.4574665763224785},
                    .temperature = 3.5884831747027144e-08};

  MomentsInCell m5 {.number_density = 4.439053529813756,
                    .bulk_velocity = {-0.022550147753564707,-0.5738356591704337,0.17698154361898316},
                    .temperature = 2.6388007584081555e-08};

  MomentsInCell m6 {.number_density = 8.764451545988388,
                    .bulk_velocity = {-0.2999960502440965,-0.09900917166706724,0.660034434776269},
                    .temperature = 3.444938074170518e-08};

  MomentsInCell m7 {.number_density = 7.2255523941379005,
                    .bulk_velocity = {0.5103292775296365,1.1467399708160815,0.5601735888875488},
                    .temperature = 5.6206806378162676e-08};

  check_moments(m0, 0);
  check_moments(m1, 1);
  check_moments(m2, 2);
  check_moments(m3, 3);
  check_moments(m4, 4);
  check_moments(m5, 5);
  check_moments(m6, 6);
  check_moments(m7, 7);

};

} // namespace
