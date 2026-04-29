#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <hdf5.h>

#include <libmfpic/DumpParticles.hpp>
#include <libmfpic/ParticleContainer.hpp>
#include <libmfpic/ParticleOperations.hpp>

#include <unordered_map>

namespace mfpic {

static bool file_is_created = false;

void dumpParticles(const ParticleContainer& particles, double simulation_time, const std::string filename) {
  hid_t file;
  if (file_is_created) {
    file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  } else {
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_is_created = true;
  }

  std::vector<double> x, y, z, vx, vy, vz, weight, particle_distribution_function_value;
  std::vector<int> element;
  std::vector<std::string> species_name; 
  for (const Particle& particle : particles) {
    if (particle.is_alive) {
      x.push_back(particle.position[0]);
      y.push_back(particle.position[1]);
      z.push_back(particle.position[2]);
      vx.push_back(particle.velocity[0]);
      vy.push_back(particle.velocity[1]);
      vz.push_back(particle.velocity[2]);
      weight.push_back(particle.weight);
      element.push_back(particle.element);
      particle_distribution_function_value.push_back(particle.particle_distribution_function_value);
      species_name.push_back(particle.species.name);  
    }
  }

  H5G_info_t top_level_group_info;
  H5Gget_info_by_name(file, "/", &top_level_group_info, H5P_DEFAULT);
  const int step_to_write = top_level_group_info.nlinks;

  std::string step_name = "Step#" + std::to_string(step_to_write);
  hid_t step_group = H5Gcreate(file, step_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  constexpr hsize_t one = 1;
  hid_t timevalue_dataspace = H5Screate_simple(1, &one, NULL);
  hid_t timevalue_attribute = H5Acreate(step_group, "TimeValue", H5T_NATIVE_DOUBLE, timevalue_dataspace, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(timevalue_attribute, H5T_NATIVE_DOUBLE, &simulation_time);
  H5Aclose(timevalue_attribute);
  H5Sclose(timevalue_dataspace);

  const hsize_t num_particles = x.size();
  hid_t dataspace = H5Screate_simple(1, &num_particles, NULL);
  hid_t dataset = H5Dcreate(step_group, "x", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "y", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "z", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vx", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vx.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vy.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "vz", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vz.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "weight", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weight.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "particle_distribution_function_value", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_distribution_function_value.data());
  H5Dclose(dataset);
  dataset = H5Dcreate(step_group, "element", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, element.data());
  H5Dclose(dataset);

  hid_t strtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype, H5T_VARIABLE);
  std::vector<const char*> cstrs;
  cstrs.reserve(species_name.size());
  for (auto &s : species_name) { cstrs.push_back(s.c_str()); }
  dataset = H5Dcreate(step_group, "species_name", strtype, dataspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cstrs.data());
  H5Dclose(dataset);
  H5Tclose(strtype);

  H5Sclose(dataspace);
  H5Gclose(step_group);
  H5Fclose(file);
}

static bool fileIsEmpty(const std::string& filename)
{
  std::ifstream in(filename, std::ios::binary | std::ios::ate);
  return (!in) || (in.tellg() == 0);
}

void dumpParticleMoments(
  ParticleOperations& particle_operations,
  const ParticleContainer& particles,
  const std::string& file_prefix,
  const int step,
  const double time) 
{
  std::unordered_map<Species, mfem::Vector> number_density     = particle_operations.getNumberDensity(particles);
  std::unordered_map<Species, mfem::DenseMatrix> bulk_velocity = particle_operations.getBulkVelocity(particles,true);
  std::unordered_map<Species, mfem::Vector> temperature        = particle_operations.getTemperature(particles,false,false);

  mfem::Mesh& mesh = particle_operations.getMesh();
  const int nelem = mesh.GetNE();

  for (const auto & [species, _] : number_density) {
    std::string filename = file_prefix + "_" + species.name + ".csv";

    std::ofstream out;
    if (step > 0)
      out.open(filename, std::ios::out | std::ios::app);
    else
      out.open(filename, std::ios::out | std::ios::trunc);  
    if (!out) throw std::runtime_error("Failed to open CSV file: " + filename);

    out.setf(std::ios::scientific);
    out << std::setprecision(17);

    const bool need_header = fileIsEmpty(filename);
    if (need_header) {
      out << "step,time,elem,x,y,z,number_density,temperature,bulk_velocity_0,bulk_velocity_1,bulk_velocity_2\n";
    }

    for (int e = 0; e < nelem; ++e) {
      mfem::Vector element_point(3);
      element_point = 0.0;
      const int dim = mesh.SpaceDimension();
      mfem::Vector element_point_view(element_point.GetData(), dim);   
      mesh.GetElementCenter(e, element_point_view);

      const mfem::Vector bulk_velocity_in_element(bulk_velocity.at(species).GetColumn(e), 3);

      out << step << ","
          << time << ","
          << e << ","
          << element_point(0) << "," << element_point(1) << "," << element_point(2) << ","
          << number_density.at(species)(e) << ","
          << temperature.at(species)(e) << ","
          << bulk_velocity_in_element(0) << "," << bulk_velocity_in_element(1) << "," << bulk_velocity_in_element(2) << "\n";
    }
  }
}

void dumpVarianceReducedParticleMoments(
  ParticleOperations& particle_operations,
  const ParticleContainer& particles,
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations,
  const std::string& file_prefix,
  const int step,
  const double time) 
{

  mfem::DenseMatrix variance_reduced_number_density = particle_operations.getVarianceReducedNumberDensity(particles,low_fidelity_state,low_fidelity_operations);
  mfem::DenseTensor variance_reduced_bulk_velocity = particle_operations.getVarianceReducedBulkVelocity(particles,low_fidelity_state,low_fidelity_operations);
  mfem::DenseMatrix variance_reduced_temperature = particle_operations.getVarianceReducedTemperature(particles,low_fidelity_state,low_fidelity_operations);

  mfem::Mesh& mesh = particle_operations.getMesh();
  const int nelem = mesh.GetNE();

  for (int s = 0; s < particle_operations.getNumSpecies(); ++s) {
    std::string filename = file_prefix + "_species_" + std::to_string(s) + ".csv";

    std::ofstream out;
    if (step > 0)
      out.open(filename, std::ios::out | std::ios::app);
    else
      out.open(filename, std::ios::out | std::ios::trunc);  
    if (!out) throw std::runtime_error("Failed to open CSV file: " + filename);

    out.setf(std::ios::scientific);
    out << std::setprecision(17);

    const bool need_header = fileIsEmpty(filename);
    if (need_header) {
      out << "step,time,elem,x,y,z,number_density,temperature,bulk_velocity_0,bulk_velocity_1,bulk_velocity_2\n";
    }

    for (int e = 0; e < nelem; ++e) {
      mfem::Vector element_point(3);
      element_point = 0.0;
      const int dim = mesh.SpaceDimension();
      mfem::Vector element_point_view(element_point.GetData(), dim);   
      mesh.GetElementCenter(e, element_point_view);

      const mfem::Vector bulk_velocity_in_element(variance_reduced_bulk_velocity(s).GetColumn(e), 3);

      out << step << ","
          << time << ","
          << e << ","
          << element_point(0) << "," << element_point(1) << "," << element_point(2) << ","
          << variance_reduced_number_density(e, s) << ","
          << variance_reduced_temperature(e, s) << ","
          << bulk_velocity_in_element(0) << "," << bulk_velocity_in_element(1) << "," << bulk_velocity_in_element(2) << "\n";
    }
  }
}

void dumpLowFidelityMoments(
  const LowFidelityState& low_fidelity_state,
  const DGEulerOperations& low_fidelity_operations,
  const std::string& file_prefix,
  const int step,
  const double time) 
{

  mfem::DenseMatrix number_density = low_fidelity_operations.getCellAveragedNumberDensity(low_fidelity_state);
  mfem::DenseTensor bulk_velocity = low_fidelity_operations.getCellAveragedBulkVelocity(low_fidelity_state);
  mfem::DenseMatrix temperature = low_fidelity_operations.getCellAveragedTemperature(low_fidelity_state);

  mfem::Mesh& mesh = low_fidelity_operations.getMesh();
  const int nelem = mesh.GetNE();

  for (int s = 0; s < low_fidelity_state.numSpecies(); ++s) {
    std::string filename = file_prefix + "_species_" + std::to_string(s) + ".csv";

    std::ofstream out;
    if (step > 0)
      out.open(filename, std::ios::out | std::ios::app);
    else
      out.open(filename, std::ios::out | std::ios::trunc);  
    if (!out) throw std::runtime_error("Failed to open CSV file: " + filename);

    out.setf(std::ios::scientific);
    out << std::setprecision(17);

    const bool need_header = fileIsEmpty(filename);
    if (need_header) {
      out << "step,time,elem,x,y,z,number_density,temperature,bulk_velocity_0,bulk_velocity_1,bulk_velocity_2\n";
    }

    for (int e = 0; e < nelem; ++e) {
      mfem::Vector element_point(3);
      element_point = 0.0;
      const int dim = mesh.SpaceDimension();
      mfem::Vector element_point_view(element_point.GetData(), dim);   
      mesh.GetElementCenter(e, element_point_view);

      const mfem::Vector bulk_velocity_in_element(bulk_velocity(s).GetColumn(e), 3);

      out << step << ","
          << time << ","
          << e << ","
          << element_point(0) << "," << element_point(1) << "," << element_point(2) << ","
          << number_density(e, s) << ","
          << temperature(e, s) << ","
          << bulk_velocity_in_element(0) << "," << bulk_velocity_in_element(1) << "," << bulk_velocity_in_element(2) << "\n";
    }
  }
}

} // namespace mfpic
