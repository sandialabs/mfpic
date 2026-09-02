import numpy as np


def arrange_particle_moment_data_by_timestep(particle_moment_data):
    rearranged_particle_moment_data = dict()

    num_steps = int(np.max(particle_moment_data['step']) + 1)
    for i in range(num_steps):
        mask = particle_moment_data['step'] == i
        rearranged_particle_moment_data[i] = particle_moment_data[mask]

    return rearranged_particle_moment_data


def generate_empty_particle_distribution_function_data(cell_edges, velocity_bins, times):
    num_cells = cell_edges.shape[0] - 1
    cell_key_width = int(np.log10(num_cells - 1)) + 1

    particle_distribution_function_data = dict()
    for i_timestep, time in enumerate(times):
        timestep_key = f'Step#{i_timestep}'

        timestep_dict = dict()
        timestep_dict['time'] = time
        timestep_dict['Number Of Runs'] = 0
        timestep_dict['cell_edges'] = cell_edges
        timestep_dict['number_density'] = np.zeros(num_cells)
        timestep_dict['velocity_bins'] = velocity_bins
        timestep_dict['cell_data'] = dict()
        for i_cell in range(num_cells):
            cell_dict = dict()
            cell_dict['f_v_weighted_sum'] = np.zeros(velocity_bins.shape[0] - 1)
            cell_dict['number_density_sum'] = 0.
            # put this in so don't have to do regular expressions to pull out of cell key
            cell_dict['i_cell'] = i_cell

            cell_key = f"Cell#{i_cell:0{cell_key_width}d}"
            timestep_dict['cell_data'][cell_key] = cell_dict

        particle_distribution_function_data[timestep_key] = timestep_dict
    return particle_distribution_function_data


def update_particle_distribution_function_data(particle_distribution_function_data, particle_data, particle_moment_data, species_name):
    for i, (timestep_key, particle_data_i) in enumerate(particle_data.items()):
        timestep_dict = particle_distribution_function_data[timestep_key]
        number_of_runs = timestep_dict["Number Of Runs"][0]

        present_number_density_average = np.array(timestep_dict['number_density'])
        number_density_to_add = particle_moment_data[i]['number_density']
        new_number_density_average = (number_of_runs * present_number_density_average + number_density_to_add) / (number_of_runs + 1)
        timestep_dict['number_density'][...] = new_number_density_average

        velocity_bins = np.array(timestep_dict['velocity_bins'])

        for cell_dict in timestep_dict['cell_data'].values():
            i_cell = cell_dict['i_cell'][0]

            mask = np.logical_and(np.array(particle_data_i['element']) == i_cell, np.array(particle_data_i['species_name']) == species_name)
            if np.sum(mask) > 0:
                velocities = np.array(particle_data_i['vx'][mask])
                weights = np.array(particle_data_i['weight'][mask])
                if (velocities.shape[0] > 0):
                    if np.min(velocities) < velocity_bins[0]:
                        print(f"np.min(velocities) = {np.min(velocities)}")
                        print(f"velocity_bins[0] = {velocity_bins[0]}")
                        print("A particle is outside of the velocity bins")
                    if np.max(velocities) > velocity_bins[-1]:
                        print(f"np.max(velocities) = {np.max(velocities)}")
                        print(f"velocity_bins[-1] = {velocity_bins[-1]}")
                        print("A particle is outside of the velocity bins")

                f_v, _ = np.histogram(velocities, bins=velocity_bins, weights=weights, density=True)

                present_f_v_weighted_sum = np.array(cell_dict['f_v_weighted_sum'])
                present_number_density_sum = cell_dict['number_density_sum'][0]

                number_density_i_cell = number_density_to_add[i_cell]
                new_f_v_weighted_sum = present_f_v_weighted_sum + number_density_i_cell * f_v
                new_number_density_sum = present_number_density_sum + number_density_i_cell

                cell_dict['f_v_weighted_sum'][...] = new_f_v_weighted_sum
                cell_dict['number_density_sum'][0] = new_number_density_sum
        timestep_dict["Number Of Runs"][0] = number_of_runs + 1

    return particle_distribution_function_data


def write_dict_to_hdf5_file(h5_group, data):
    """Recursively saves a dictionary of numpy arrays to an HDF5 group."""
    for key, value in data.items():
        if isinstance(value, dict):
            # Create a group for a nested sub-dictionary
            sub_group = h5_group.create_group(key)
            write_dict_to_hdf5_file(sub_group, value)
        else:
            # Create a dataset for the array
            h5_group.create_dataset(key, data=np.atleast_1d(value))


def compute_error_in_stored_f_vs_exact_f(particle_data, species_name, exact_particle_distribution_function, times):
    error_over_time = []
    for i, particle_data_i in enumerate(particle_data.values()):
        time = times[i]

        mask = np.array(particle_data_i["species_name"]) == species_name
        num_macroparticles = np.sum(mask)

        particle_positions = np.array(particle_data_i['x'][mask])
        particle_velocities = np.array(particle_data_i['vx'][mask])
        particle_f = np.array(particle_data_i['particle_distribution_function_value'][mask])

        exact_f = exact_particle_distribution_function(particle_positions, particle_velocities, time)

        error = np.sum(np.abs(particle_f - exact_f)) / num_macroparticles
        error_over_time.append(error)

    return np.array(error_over_time)


def compute_error_in_stored_f_vs_histogram_f(particle_data, species_name, particle_distribution_function_data):
    error_over_time = []
    for timestep_key, particle_data_i in particle_data.items():
        timestep_dict = particle_distribution_function_data[timestep_key]
        number_density = np.array(timestep_dict['number_density'])
        velocity_bins = np.array(timestep_dict['velocity_bins'])

        error = 0
        species_mask = np.array(particle_data_i["species_name"]) == species_name
        num_macroparticles = np.sum(species_mask)

        for cell_dict in timestep_dict['cell_data'].values():
            i_cell = cell_dict['i_cell'][0]

            mask = np.logical_and(species_mask, np.array(particle_data_i['element']) == i_cell)

            particle_velocities = np.array(particle_data_i['vx'][mask])
            particle_f = np.array(particle_data_i['particle_distribution_function_value'][mask])

            number_density_sum = cell_dict['number_density_sum'][0]
            f_v_weighted_sum = np.array(cell_dict['f_v_weighted_sum'])
            f_v = f_v_weighted_sum / number_density_sum
            f = number_density[i_cell] * f_v

            particle_bin_indices = np.digitize(particle_velocities, velocity_bins) - 1
            histogram_f = f[particle_bin_indices]

            error += np.sum(np.abs(particle_f - histogram_f)) / num_macroparticles

        error_over_time.append(error)

    return np.array(error_over_time)