#include <libmfpic/BGKRelaxationOperations.hpp>
#include <libmfpic/BuildCollisionOperationsFromYaml.hpp>
#include <libmfpic/Errors.hpp>
#include <libmfpic/ParseHelpers.hpp>

namespace mfpic {

std::vector<std::unique_ptr<CollisionOperations>> buildCollisionOperationsFromYaml(
  const YAML::Node& collisions_nodes,
  const std::unordered_map<std::string, Species>& species_map
) {
  std::vector<std::unique_ptr<CollisionOperations>> collision_operations;
  if (collisions_nodes.IsSequence()) {
    for (const YAML::Node& collisions_node : collisions_nodes) {
      if (collisions_node["BGK Relaxation"]) {
        const YAML::Node bgk_relaxation_node = collisions_node["BGK Relaxation"];

        const double collision_frequency = parseRequiredPositiveDouble(bgk_relaxation_node, "Collision Frequency");
        const std::string species_name = parseSpeciesName(bgk_relaxation_node["Species"], species_map);

        const double averaging_timesteps = bgk_relaxation_node["Moment Averaging Timesteps"].as<double>(1.0);
        if (averaging_timesteps < 1.0) {
          errorWithUserMessage(formatParseMessage(bgk_relaxation_node["Moment Averaging Timesteps"], "Moment Averaging Timesteps must be at least 1!"));
        }
        const double smoothing_length = bgk_relaxation_node["Moment Smoothing Length"].as<double>(0.0);
        if (smoothing_length < 0.0) {
          errorWithUserMessage(formatParseMessage(bgk_relaxation_node["Moment Smoothing Length"], "Moment Smoothing Length must be nonnegative!"));
        }
        const int smoothing_order = bgk_relaxation_node["Moment Smoothing Order"].as<int>(0);
        if (smoothing_order < 0 or smoothing_order > 2) {
          errorWithUserMessage(formatParseMessage(bgk_relaxation_node["Moment Smoothing Order"], "Moment Smoothing Order must be 0, 1, or 2!"));
        }
        const bool use_vr_moments = bgk_relaxation_node["Use VR Moments"].as<bool>(false);
        const bool log_cell_energy = bgk_relaxation_node["Log Cell Energy Change"].as<bool>(false);
        const bool stamp_raw_density = bgk_relaxation_node["Stamp Raw Density"].as<bool>(false);
        const bool ratio_stamp = bgk_relaxation_node["Ratio Stamp"].as<bool>(false);
        const int density_scale_stride = bgk_relaxation_node["Density Scale Stride"].as<int>(0);
        const bool invariant_constraint = bgk_relaxation_node["Invariant Constraint"].as<bool>(false);
        const int kde_neighbors = bgk_relaxation_node["Weight KDE Neighbors"].as<int>(0);
        const bool relaxed_mass_renormalization = bgk_relaxation_node["Relaxed Mass Renormalization"].as<bool>(false);
        const bool kde_unit_scale = bgk_relaxation_node["KDE Unit Scale"].as<bool>(false);
        collision_operations.push_back(std::make_unique<BGKRelaxationOperations>(
          collision_frequency,
          species_map.at(species_name),
          averaging_timesteps,
          smoothing_length,
          smoothing_order,
          use_vr_moments,
          log_cell_energy,
          stamp_raw_density,
          ratio_stamp,
          density_scale_stride,
          invariant_constraint,
          kde_neighbors,
          relaxed_mass_renormalization,
          kde_unit_scale
        ));
      }
      else {
        errorWithUserMessage(formatParseMessage(collisions_node, "Unrecognized collision type!"));
      }
    }
  }

  return collision_operations;
}

} // namespace mfpic
