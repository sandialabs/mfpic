import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate


def create_1D_interpolater(field_data, points):
    # points is (2 * n_cells,), and is (x_0, x_1, x_1, x_2, x_2, x_3, ...)
    # field_data is (2 * n_cells,) and is (f(x_0), f_l(x_1), f_r(x_1), ...)
    # i.e. field data could be discontinuous at points

    n_cells = int(0.5 * points.shape[0])
    points_by_cell = points.reshape(n_cells, 2)
    field_by_cell = field_data.reshape(n_cells, 2)

    def interpolater(x):
        result = 0 if np.isscalar(x) else np.zeros(x.shape)
        for i_cell in range(n_cells):
            x_left = points_by_cell[i_cell, 0]
            x_right = points_by_cell[i_cell, 1]
            mask = (x >= x_left) * (x < x_right) if i_cell < n_cells - 1 else (x >= x_left) * (x <= x_right)
            f_left = field_by_cell[i_cell, 0]
            f_right = field_by_cell[i_cell, 1]
            line = f_left * (x - x_right) / (x_left - x_right) + f_right * (x - x_left) / (x_right - x_left)
            result += line * mask
        return result

    return interpolater


def integrate_over_mesh_1D(function, points):
    # function should be function of x
    # points is (2 * n_cells) and looks like (x_0, x_1, x_1, x_2, ...)
    n_cells = int(0.5 * points.shape[0])
    points_by_cell = points.reshape(n_cells, 2)

    integral = 0
    for i_cell in range(n_cells):
        x_left = points_by_cell[i_cell, 0]
        x_right = points_by_cell[i_cell, 1]
        cell_integral, _ = scipy.integrate.fixed_quad(function, x_left, x_right)
        integral += cell_integral

    return integral


def compute_L1_norm_1D(function, points):
    function_to_integrate = lambda x: np.abs(function(x))
    L1_norm = integrate_over_mesh_1D(function_to_integrate, points)
    return L1_norm


def compute_L1_error_1D(numerical_solution, exact_solution, points):
    error_function = lambda x: numerical_solution(x) - exact_solution(x)
    return compute_L1_norm_1D(error_function, points)


def compute_L1_relative_error_1D(numerical_solution, exact_solution, points):
    L1_error = compute_L1_error_1D(numerical_solution, exact_solution, points)
    exact_solution_norm = compute_L1_norm_1D(exact_solution, points)
    L1_relative_error = L1_error / exact_solution_norm
    return L1_relative_error


def compute_L2_norm_1D(function, points):
    # function should be function of x
    # points is (2 * n_cells) and looks like (x_0, x_1, x_1, x_2, ...)

    function_to_integrate = lambda x: np.power(function(x), 2)
    integrand = integrate_over_mesh_1D(function_to_integrate, points)
    norm = np.sqrt(integrand)
    return norm


def compute_L2_error_1D(numerical_solution, exact_solution, points):
    # numerical_solution and exact_solution should be functions of x
    # points is (2 * n_cells) and looks like (x_0, x_1, x_1, x_2, ...)

    error_function = lambda x: numerical_solution(x) - exact_solution(x)
    return compute_L2_norm_1D(error_function, points)


def compute_L2_relative_error_1D(numerical_solution, exact_solution, points):
    L2_error = compute_L2_error_1D(numerical_solution, exact_solution, points)
    exact_solution_norm = compute_L2_norm_1D(exact_solution, points)
    L2_relative_error = L2_error / exact_solution_norm
    return L2_relative_error


def compute_convergence_rates(errors, h_array):
    errors_np = np.array(errors)
    h_np = np.array(h_array)
    return np.log(errors_np[:-1] / errors_np[1:]) / np.log(h_np[:-1] / h_np[1:])


def plot_errors_and_expected_convergence_rate(refinement_quantities, errors, expected_rate, figname="./error_convergence.png"):
    h_log_mean = np.exp(np.mean(np.log(refinement_quantities)))
    error_log_mean = np.exp(np.mean(np.log(errors)))
    c = error_log_mean / np.power(h_log_mean, expected_rate) # error = c*h^(expected_rate)
    h_begin = refinement_quantities[0]
    h_end = refinement_quantities[-1]
    error_begin = c * np.power(h_begin, expected_rate)
    error_end = c * np.power(h_end, expected_rate)
    fig, ax = plt.subplots()
    ax.loglog([h_begin, h_end], [error_begin, error_end], label=f"Expected Rate of {expected_rate}", color='tab:orange', linewidth=2)
    ax.loglog(refinement_quantities, errors, label='Observed Rate', color='tab:blue')
    ax.legend()
    ax.set_xlabel("Refinement Quantity (h)")
    ax.set_ylabel("Error")
    fig.savefig(figname, bbox_inches='tight')
    plt.close(fig)


def check_total_fluid_energy_positive_and_constant():
    txt_data = np.genfromtxt('output_lf_0.csv', names=True)
    total_fluid_energy = txt_data['Total_Fluid_Energy']

    assert(np.min(total_fluid_energy) >= 0)

    change_in_total_fluid_energy = np.max(total_fluid_energy) - np.min(total_fluid_energy)

    if (np.max(total_fluid_energy) > 0):
        relative_change_in_total_fluid_energy = change_in_total_fluid_energy / np.max(total_fluid_energy)
    else:
        relative_change_in_total_fluid_energy = 0

    relative_tolerance = 1e-14
    assert(relative_change_in_total_fluid_energy < relative_tolerance)