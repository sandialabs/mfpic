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
            mask = (
                (x >= x_left) * (x < x_right)
                if i_cell < n_cells - 1
                else (x >= x_left) * (x <= x_right)
            )
            f_left = field_by_cell[i_cell, 0]
            f_right = field_by_cell[i_cell, 1]
            line = f_left * (x - x_right) / (x_left - x_right) + f_right * (
                x - x_left
            ) / (x_right - x_left)
            result += line * mask
        return result

    return interpolater

def compute_L2_norm_1D(function, points):
    # function should be function of x
    # points is (2 * n_cells) and looks like (x_0, x_1, x_1, x_2, ...)
    n_cells = int(0.5 * points.shape[0])
    points_by_cell = points.reshape(n_cells, 2)

    function_to_integrate = lambda x: np.power(function(x), 2)
    integrand = 0
    for i_cell in range(n_cells):
        x_left = points_by_cell[i_cell, 0]
        x_right = points_by_cell[i_cell, 1]
        cell_integrand, _ = scipy.integrate.quad(function_to_integrate, x_left, x_right)
        integrand += cell_integrand

    norm = np.sqrt(integrand)
    return norm

def compute_L2_error_1D(numerical_solution, exact_solution, points):
    # numerical_solution and exact_solution should be functions of x
    # points is (2 * n_cells) and looks like (x_0, x_1, x_1, x_2, ...)

    error_function = lambda x: numerical_solution(x) - exact_solution(x)
    return compute_L2_norm_1D(error_function)

def compute_L2_relative_error_1D(numerical_solution, exact_solution, points):
    L2_error = compute_L2_error_1D(numerical_solution, exact_solution, points)
    exact_solution_norm = compute_L2_norm_1D(exact_solution, points)
    L2_relative_error = L2_error / exact_solution_norm
    return L2_relative_error

def compute_convergence_rates(errors, h_array):
    return np.log(errors[:-1] / errors[1:]) / np.log(h_array[:-1] / h_array[1:])
