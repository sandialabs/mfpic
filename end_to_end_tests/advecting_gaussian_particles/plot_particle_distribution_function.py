import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.integrate
from scipy.constants import Boltzmann as boltzmann_constant

def exact_number_density(x):
    gaussian_standard_deviation = 0.1
    gaussian_center = 4 * gaussian_standard_deviation
    domain_length = 9 * gaussian_standard_deviation

    number_density_offset = 1e26
    perturbation = 1.0
    number_density_height = perturbation * number_density_offset

    exact_number_density = number_density_height * np.exp(-0.5 * np.power(x - gaussian_center, 2) / np.power(gaussian_standard_deviation, 2)) + number_density_offset
    return exact_number_density


def exact_number_density_cell_average(x_left, x_right):
    cell_size = x_right - x_left
    integral, _ = scipy.integrate.quad(exact_number_density, x_left, x_right)
    cell_average = integral / cell_size
    return cell_average


def maxwellian_1d(velocity, number_density):
    pressure = 27613
    temperature = pressure / (number_density * boltzmann_constant)
    mass = 4.65e-26
    sigma = np.sqrt(boltzmann_constant * temperature / mass)
    inv_sq_sigma = 1.0 / (sigma * sigma)
    bulk_velocity = 161.1830422788935
    difference = velocity - bulk_velocity
    exponent = inv_sq_sigma * np.power(difference, 2)
    norm = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
    probability_density_function = norm * np.exp(-0.5 * exponent)
    return probability_density_function


def plot_particle_distribution_function(particle_data, particle_distribution_function_data):
    figures_directory = "Figures"
    os.makedirs(figures_directory, exist_ok=True)

    for timestep_key, timestep_dict in particle_distribution_function_data.items():
        figures_subdirectory = f"Figures/{timestep_key}"
        os.makedirs(figures_subdirectory, exist_ok=True)

        number_density = timestep_dict['number_density']
        cell_edges = timestep_dict['cell_edges']

        fig, axes = plt.subplots()
        axes.stairs(number_density, cell_edges, baseline=None, label='Numerical')
        if timestep_key == "Step#0":
            x = np.linspace(cell_edges[0], cell_edges[-1], 200)
            axes.plot(x, exact_number_density(x), label='Exact')
        axes.legend()
        axes.set_xlabel("x")
        axes.set_ylabel("n")
        fig.savefig(f"{figures_subdirectory}/NumberDensity.png")
        plt.close(fig)

        particle_data_i = particle_data[timestep_key]

        for i_cell, (cell_key, cell_dict) in enumerate(timestep_dict['cell_data'].items()):
            mask = np.array(particle_data_i['element']) == i_cell
            bin_edges = cell_dict['bin_edges']

            n = exact_number_density_cell_average(cell_edges[i_cell], cell_edges[i_cell + 1])
            v_linspace = np.linspace(bin_edges[0], bin_edges[-1], 200)
            f_v = maxwellian_1d(v_linspace, n)

            fig, axes = plt.subplots()
            axes.stairs(number_density[i_cell] * cell_dict['f_v'], cell_dict['bin_edges'], linewidth=2, label='Refine')
            axes.plot(particle_data_i['vx'][mask], particle_data_i['particle_distribution_function_value'][mask], 'o', label='Stored')
            if timestep_key == "Step#0":
                axes.plot(v_linspace, n * f_v, label='Exact')
            axes.legend()
            axes.set_xlabel("vx")
            axes.set_ylabel("f")
            fig.savefig(f"{figures_subdirectory}/{cell_key}.png")
            plt.close(fig)

            fig, axes = plt.subplots()
            axes.stairs(number_density[i_cell] * cell_dict['f_v'], cell_dict['bin_edges'], linewidth=2, label='Refined')
            axes.plot(particle_data_i['vx'][mask], particle_data_i['particle_distribution_function_value'][mask], 'o', label='Stored')
            if timestep_key == "Step#0":
                axes.plot(v_linspace, n * f_v, label='Exact')
            axes.legend()
            axes.set_xlabel("vx")
            axes.set_ylabel("f")
            axes.set_yscale('log')
            fig.savefig(f"{figures_subdirectory}/{cell_key}Log.png")
            plt.close(fig)

particle_distribution_function_data = h5py.File("particle_distribution_function_data.h5")
particle_data = h5py.File("particles.h5part")
plot_particle_distribution_function(particle_data, particle_distribution_function_data)