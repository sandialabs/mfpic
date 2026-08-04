import h5py
import matplotlib.pyplot as plt
import numpy as np
import os


def arrange_particle_moment_data_by_timestep(particle_moment_data):
    rearranged_particle_moment_data = dict()

    num_steps = int(np.max(particle_moment_data['step']) + 1)
    for i in range(num_steps):
        mask = particle_moment_data['step'] == i
        rearranged_particle_moment_data[i] = particle_moment_data[mask]

    return rearranged_particle_moment_data


def generate_particle_distribution_function_data(particle_data, particle_moment_data):
    num_cells = int(np.max(particle_moment_data[0]['elem']) + 1)
    cell_size = particle_moment_data[0]['x'][1] - particle_moment_data[0]['x'][0]
    domain_length = particle_moment_data[0]['x'][-1] + 0.5 * cell_size
    cell_edges = np.linspace(0, domain_length, num_cells + 1)

    num_particles = particle_data['Step#0']['element'].shape[0]
    approximate_num_particles_per_cell = num_particles / num_cells
    # num_velocity_bins = int(approximate_num_particles_per_cell / 10)
    num_velocity_bins = 30

    particle_distribution_function_data = dict()
    for i, (timestep_key, particle_data_i) in enumerate(particle_data.items()):
        timestep_dict = dict()
        timestep_dict['cell_edges'] = cell_edges
        timestep_dict['number_density'] = particle_moment_data[i]['number_density']
        timestep_dict['cell_data'] = dict()
        for i_cell in range(num_cells):
            cell_dict = dict()
            mask = np.array(particle_data_i['element']) == i_cell
            velocities = particle_data_i['vx'][mask]
            weights = particle_data_i['weight'][mask]
            f_v, bin_edges = np.histogram(velocities, bins=num_velocity_bins, weights=weights, density=True)
            cell_dict['f_v'] = f_v
            cell_dict['bin_edges'] = bin_edges

            cell_key = f"Cell#{i_cell}"
            timestep_dict['cell_data'][cell_key] = cell_dict

        particle_distribution_function_data[timestep_key] = timestep_dict
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


def plot_particle_distribution_function(particle_data, particle_distribution_function_data):
    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    for timestep_key, timestep_dict in particle_distribution_function_data.items():
        figures_subdirectory = f"Figures/{timestep_key}"
        os.makedirs(figures_subdirectory, exist_ok=True)

        number_density = timestep_dict['number_density']
        fig, axes = plt.subplots()
        axes.stairs(number_density, timestep_dict['cell_edges'], baseline=None)
        axes.set_xlabel("x")
        axes.set_ylabel("n")
        fig.savefig(f"{figures_subdirectory}/NumberDensity.png")
        plt.close(fig)

        particle_data_i = particle_data[timestep_key]

        for i_cell, (cell_key, cell_dict) in enumerate(timestep_dict['cell_data'].items()):
            mask = np.array(particle_data_i['element']) == i_cell

            fig, axes = plt.subplots()
            axes.stairs(number_density[i_cell] * cell_dict['f_v'], cell_dict['bin_edges'])
            axes.plot(particle_data_i['vx'][mask], particle_data_i['particle_distribution_function_value'][mask], 'bo')
            axes.set_xlabel("vx")
            axes.set_ylabel("f")
            axes.set_yscale('log')
            fig.savefig(f"{figures_subdirectory}/{cell_key}.png")
            plt.close(fig)


particle_data = h5py.File("particles.h5part")
particle_moment_data = np.genfromtxt('particle_moments_N2.csv', names=True, delimiter=',')
particle_moment_data = arrange_particle_moment_data_by_timestep(particle_moment_data)
particle_distribution_function_data = generate_particle_distribution_function_data(particle_data, particle_moment_data)

with h5py.File("particle_distribution_function_data.h5", 'w') as hf:
    write_dict_to_hdf5_file(hf, particle_distribution_function_data)

plot_particle_distribution_function(particle_data, particle_distribution_function_data)

# particle_distribution_function_data = h5py.File("particle_distribution_function_data.h5")
# particle_data = h5py.File("particles.h5part")
# plot_particle_distribution_function(particle_data, particle_distribution_function_data)