#include <libmfpic/Discretization.hpp>
#include <libmfpic/ElectromagneticFieldsEvaluator.hpp>
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

} // namespace
