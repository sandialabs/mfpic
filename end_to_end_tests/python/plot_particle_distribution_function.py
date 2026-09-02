import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.integrate


def plot_particle_distribution_function(particle_distribution_function_data, exact_particle_distribution_function=None, v_max=None):
    for timestep_key, timestep_dict in particle_distribution_function_data.items():
        figures_directory = f"Figures/{timestep_key}/HistogramF"
        os.makedirs(figures_directory, exist_ok=True)

        number_density = np.array(timestep_dict['number_density'])
        cell_edges = np.array(timestep_dict['cell_edges'])
        velocity_bins = np.array(timestep_dict['velocity_bins'])
        time = timestep_dict['time'][0]

        for cell_key, cell_dict in timestep_dict['cell_data'].items():
            i_cell = cell_dict['i_cell'][0]

            number_density_sum = cell_dict['number_density_sum'][0]
            f_v_weighted_sum = np.array(cell_dict['f_v_weighted_sum'])
            f_v = f_v_weighted_sum / number_density_sum

            fig, axes = plt.subplots()
            axes.stairs(number_density[i_cell] * f_v, velocity_bins, linewidth=2, label='Histogram')
            if exact_particle_distribution_function:
                cell_center = 0.5 * (cell_edges[i_cell] + cell_edges[i_cell + 1])
                v_plot = np.linspace(velocity_bins[0], velocity_bins[-1], 400)
                f_v_exact = exact_particle_distribution_function(cell_center, v_plot, time)
                axes.plot(v_plot, f_v_exact, label=f'Exact at x = {cell_center:.2e}')

            axes.legend()
            axes.set_xlabel("vx")
            axes.set_ylabel("f")
            file_postfix = ""
            if v_max:
                axes.set_xlim((-v_max, v_max))
                file_postfix += "_zoom"

            fig.savefig(f"{figures_directory}/{cell_key}{file_postfix}.png")

            axes.set_yscale('log')
            file_postfix += "_log"
            fig.savefig(f"{figures_directory}/{cell_key}{file_postfix}.png")

            plt.close(fig)


def plot_particle_distribution_function_vs_stored_f(particle_data, particle_distribution_function_data, species_name, exact_particle_distribution_function=None):
    min_v = np.min([np.min(particle_data_i['vx']) for particle_data_i in particle_data.values()])
    max_v = np.max([np.max(particle_data_i['vx']) for particle_data_i in particle_data.values()])
    v_plot = np.linspace(min_v, max_v, 200)

    for timestep_key, particle_data_i in particle_data.items():
        figures_directory = f"Figures/{timestep_key}/ParticleStoredfVsExactf"
        os.makedirs(figures_directory, exist_ok=True)

        timestep_dict = particle_distribution_function_data[timestep_key]
        time = timestep_dict['time'][0]
        number_density = np.array(timestep_dict['number_density'])
        cell_edges = np.array(timestep_dict['cell_edges'])
        velocity_bins = np.array(timestep_dict['velocity_bins'])

        for cell_key, cell_dict in timestep_dict['cell_data'].items():
            i_cell = cell_dict['i_cell'][0]

            mask = np.logical_and(np.array(particle_data_i['element']) == i_cell, np.array(particle_data_i['species_name']) == species_name)
            particle_velocities = particle_data_i['vx'][mask]
            particle_f = particle_data_i['particle_distribution_function_value'][mask]

            number_density_sum = cell_dict['number_density_sum'][0]
            f_v_weighted_sum = np.array(cell_dict['f_v_weighted_sum'])
            f_v = f_v_weighted_sum / number_density_sum

            fig, axes = plt.subplots()
            axes.stairs(number_density[i_cell] * f_v, velocity_bins, linewidth=2, label='Histogram')
            if exact_particle_distribution_function:
                cell_center = 0.5 * (cell_edges[i_cell] + cell_edges[i_cell + 1])
                f_v_exact = exact_particle_distribution_function(cell_center, v_plot, time)
                axes.plot(v_plot, f_v_exact, label=f'Exact at x = {cell_center:.2e}')
            axes.plot(particle_velocities, particle_f, 'o', label="Particles")

            axes.legend()
            axes.set_xlabel("v")
            axes.set_ylabel("f")
            axes.set_xlim((min_v, max_v))
            fig.savefig(f"{figures_directory}/{cell_key}.png")
            plt.close(fig)


def plot_number_density(particle_distribution_function_data, exact_number_density=None):
    for timestep_key, timestep_dict in particle_distribution_function_data.items():
        figures_directory = f"Figures/{timestep_key}"
        os.makedirs(figures_directory, exist_ok=True)

        time = timestep_dict['time'][0]
        number_density = timestep_dict['number_density']
        cell_edges = timestep_dict['cell_edges']

        fig, axes = plt.subplots()
        axes.stairs(number_density, cell_edges, baseline=None, linewidth=2, label="Numerical")
        if exact_number_density:
            x_plot = np.linspace(cell_edges[0], cell_edges[-1], 100)
            n_exact = np.array([exact_number_density(x, time) for x in x_plot])
            axes.plot(x_plot, n_exact, label="Exact")
        axes.legend()
        axes.set_xlabel("x")
        axes.set_ylabel("number density")
        fig.savefig(f"{figures_directory}/NumberDensity.png")
        plt.close(fig)