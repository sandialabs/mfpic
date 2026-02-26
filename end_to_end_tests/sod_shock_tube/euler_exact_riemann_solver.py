import euler

import numpy as np
from scipy.optimize import fsolve

# The exact riemann solution to the Euler equations consists of three waves.
# There is either a shock or rarefaction on the left, a contact discontinuity in the
# middle, and a shock or rarefaction on the right.
# The pressure and velocity is constant through the contact discontinuity, i.e. only the
# density jumps over the contact discontinuity.
# Exactly solving the riemann problem for the Euler equations, thus consists of solving
# for the pressure and velocity in the middle state, solving for the density to the left
# and right of the contact discontinuity, determining if the left and right waves are
# shocks or rarefactions, and what speeds they are propogating at.
# With all of this information the exact solution to the Riemann problem can be
# reconstructed.
# See Finite Volume Methods for Hyperbolic Problems by Randall Leveque section 14.11 for
# more details


def unpack_state(primitive_state, species):
    mass_density = euler.get_mass_density_from_primitive_state(primitive_state, species)
    velocity = euler.get_bulk_velocity_from_primitive_state(primitive_state)[0]
    pressure = euler.get_pressure_from_primitive_state(primitive_state)
    return mass_density, velocity, pressure


# The middle state has a constant pressure and velocity
# which can be solved for using the left and right states.
def solve_for_middle_pressure_and_velocity(
    left_state_primitive, right_state_primitive, species
):
    mass_density_left, velocity_left, pressure_left = unpack_state(
        left_state_primitive, species
    )
    sound_speed_left = euler.speed_of_sound(species, mass_density_left, pressure_left)

    mass_density_right, velocity_right, pressure_right = unpack_state(
        right_state_primitive, species
    )
    sound_speed_right = euler.speed_of_sound(
        species, mass_density_right, pressure_right
    )

    gamma = species.specific_heat_ratio

    # if there is a rarefaction middle velocity is given by integral curves
    exponent = (gamma - 1.0) / (2.0 * gamma)
    integral_curve_left = lambda pressure: velocity_left + 2 * sound_speed_left / (
        gamma - 1.0
    ) * (1.0 - np.power(pressure / pressure_left, exponent))
    integral_curve_right = lambda pressure: velocity_right - 2 * sound_speed_right / (
        gamma - 1.0
    ) * (1.0 - np.power(pressure / pressure_right, exponent))

    # if there is a shock then middle velocity is given by hugoniot locus
    beta = (gamma + 1.0) / (gamma - 1)
    denominator = np.sqrt(2 * gamma * (gamma - 1.0))
    hugoniot_locus_left = (
        lambda pressure: velocity_left
        + 2
        * sound_speed_left
        / denominator
        * (
            (1 - pressure / pressure_left)
            / (np.sqrt(1 + beta * pressure / pressure_left))
        )
    )
    hugoniot_locus_right = (
        lambda pressure: velocity_right
        - 2
        * sound_speed_right
        / denominator
        * (
            (1 - pressure / pressure_right)
            / (np.sqrt(1 + beta * pressure / pressure_right))
        )
    )

    # middle velocity must be the same when computed from left state and right state
    middle_velocity_left = lambda pressure: (
        hugoniot_locus_left(pressure)
        if pressure >= pressure_left
        else integral_curve_left(pressure)
    )
    middle_velocity_right = lambda pressure: (
        hugoniot_locus_right(pressure)
        if pressure >= pressure_right
        else integral_curve_right(pressure)
    )

    func = lambda pressure: middle_velocity_left(pressure) - middle_velocity_right(
        pressure
    )
    initial_guess = 0.5 * (pressure_left + pressure_right)
    sol = fsolve(func, initial_guess, xtol=1e-14)
    pressure_middle = sol[0]

    velocity_middle = middle_velocity_left(pressure_middle)
    return pressure_middle, velocity_middle


def check_is_rarefaction(pressure_middle, pressure):
    return pressure_middle < pressure


def mass_density_if_rarefaction(pressure_middle, pressure, mass_density, species):
    exponent = 1.0 / species.specific_heat_ratio
    return np.power(pressure_middle / pressure, exponent) * mass_density


def mass_density_if_shock(pressure_middle, pressure, mass_density, species):
    beta = (species.specific_heat_ratio + 1.0) / (species.specific_heat_ratio - 1)
    return (
        (1.0 + beta * pressure_middle / pressure) / (pressure_middle / pressure + beta)
    ) * mass_density


def compute_middle_mass_density(primitive_state, species, pressure_middle):
    mass_density, _, pressure = unpack_state(primitive_state, species)
    is_rarefaction = check_is_rarefaction(pressure_middle, pressure)

    mass_density_rarefaction = mass_density_if_rarefaction(
        pressure_middle, pressure, mass_density, species
    )
    mass_density_shock = mass_density_if_shock(
        pressure_middle, pressure, mass_density, species
    )
    mass_density_middle = (
        mass_density_rarefaction if is_rarefaction else mass_density_shock
    )
    return mass_density_middle


def shock_speed(mass_density_left, velocity_left, mass_density_right, velocity_right):
    return (mass_density_left * velocity_left - mass_density_right * velocity_right) / (
        mass_density_left - mass_density_right
    )


def compute_wavespeeds(
    left_state_primitive,
    right_state_primitive,
    species,
    pressure_middle,
    velocity_middle,
    mass_density_middle_left,
    mass_density_middle_right,
):
    # return 5 wavespeeds associated with exact riemann solution in ascending order
    # (smallest_left, biggest_left, contact_speed, smallest_right, biggest_right)
    # If the left/right wave is a shock,
    # then smallest_left/right = biggest_left/right = shock_speed
    # If the left/right wave is a rarefaction,
    # then smallest_left/right and biggest_left/right
    # give the speeds of the extents of the rarefaction

    contact_speed = velocity_middle

    mass_density_left = euler.get_mass_density_from_primitive_state(
        left_state_primitive, species
    )
    velocity_left = euler.get_bulk_velocity_from_primitive_state(left_state_primitive)[
        0
    ]
    pressure_left = euler.get_pressure_from_primitive_state(left_state_primitive)

    mass_density_right = euler.get_mass_density_from_primitive_state(
        right_state_primitive, species
    )
    velocity_right = euler.get_bulk_velocity_from_primitive_state(
        right_state_primitive
    )[0]
    pressure_right = euler.get_pressure_from_primitive_state(right_state_primitive)

    shock_speed_left = shock_speed(
        mass_density_left, velocity_left, mass_density_middle_left, velocity_middle
    )

    sound_speed_left = euler.speed_of_sound(species, mass_density_left, pressure_left)
    sound_speed_middle_left = euler.speed_of_sound(
        species, mass_density_middle_left, pressure_middle
    )
    rarefaction_smallest_left = velocity_left - sound_speed_left
    rarefaction_biggest_left = velocity_middle - sound_speed_middle_left

    shock_speed_right = shock_speed(
        mass_density_middle_right, velocity_middle, mass_density_right, velocity_right
    )

    sound_speed_right = euler.speed_of_sound(
        species, mass_density_right, pressure_right
    )
    sound_speed_middle_right = euler.speed_of_sound(
        species, mass_density_middle_right, pressure_middle
    )
    rarefaction_smallest_right = velocity_middle + sound_speed_middle_right
    rarefaction_biggest_right = velocity_right + sound_speed_right

    is_rarefaction_left = check_is_rarefaction(pressure_middle, pressure_left)
    is_rarefaction_right = check_is_rarefaction(pressure_middle, pressure_right)

    smallest_left = (
        rarefaction_smallest_left if is_rarefaction_left else shock_speed_left
    )
    biggest_left = rarefaction_biggest_left if is_rarefaction_left else shock_speed_left

    smallest_right = (
        rarefaction_smallest_right if is_rarefaction_right else shock_speed_right
    )
    biggest_right = (
        rarefaction_biggest_right if is_rarefaction_right else shock_speed_right
    )

    return smallest_left, biggest_left, contact_speed, smallest_right, biggest_right

def compute_max_wavespeed(left_state_primitive, right_state_primitive, species):
    pressure_middle, velocity_middle = solve_for_middle_pressure_and_velocity(
        left_state_primitive, right_state_primitive, species
    )
    mass_density_middle_left = compute_middle_mass_density(
        left_state_primitive, species, pressure_middle
    )
    mass_density_middle_right = compute_middle_mass_density(
        right_state_primitive, species, pressure_middle
    )

    wavespeeds = compute_wavespeeds(
        left_state_primitive,
        right_state_primitive,
        species,
        pressure_middle,
        velocity_middle,
        mass_density_middle_left,
        mass_density_middle_right,
    )

    return np.max(np.abs(wavespeeds))

def form_rarefaction_solutions(
    primitive_state, species, discontinuity_location, is_left_going
):
    gamma = species.specific_heat_ratio
    mass_density, velocity, pressure = unpack_state(primitive_state, species)
    sound_speed = euler.speed_of_sound(species, mass_density, pressure)

    shift = lambda x, t: (x - discontinuity_location) / t if t > 0 else 0

    sign = 1 if is_left_going else -1
    velocity_rarefaction = lambda x, t: (
        (gamma - 1.0) * velocity + 2.0 * (shift(x, t) + sign * sound_speed)
    ) / (gamma + 1.0)

    exponent = 1.0 / (gamma - 1.0)
    scale = np.power(mass_density, gamma) / (gamma * pressure)
    mass_density_rarefaction = lambda x, t: np.power(
        scale * np.power(shift(x, t) - velocity_rarefaction(x, t), 2), exponent
    )
    pressure_rarefaction = lambda x, t: pressure * np.power(
        mass_density_rarefaction(x, t) / mass_density, gamma
    )

    return mass_density_rarefaction, velocity_rarefaction, pressure_rarefaction


def form_exact_solutions(
    left_state_primitive, right_state_primitive, species, discontinuity_location
):
    mass_density_left, velocity_left, pressure_left = unpack_state(
        left_state_primitive, species
    )
    mass_density_right, velocity_right, pressure_right = unpack_state(
        right_state_primitive, species
    )

    pressure_middle, velocity_middle = solve_for_middle_pressure_and_velocity(
        left_state_primitive, right_state_primitive, species
    )
    mass_density_middle_left = compute_middle_mass_density(
        left_state_primitive, species, pressure_middle
    )
    mass_density_middle_right = compute_middle_mass_density(
        right_state_primitive, species, pressure_middle
    )

    wavespeeds = compute_wavespeeds(
        left_state_primitive,
        right_state_primitive,
        species,
        pressure_middle,
        velocity_middle,
        mass_density_middle_left,
        mass_density_middle_right,
    )

    (
        mass_density_rarefaction_left,
        velocity_rarefaction_left,
        pressure_rarefaction_left,
    ) = form_rarefaction_solutions(
        left_state_primitive, species, discontinuity_location, True
    )
    (
        mass_density_rarefaction_right,
        velocity_rarefaction_right,
        pressure_rarefaction_right,
    ) = form_rarefaction_solutions(
        right_state_primitive, species, discontinuity_location, False
    )

    d = discontinuity_location
    w = wavespeeds

    def mass_density_solution(x, t):
        left = mass_density_left * (x - d <= w[0] * t)
        left_rarefaction = (
            mass_density_rarefaction_left(x, t)
            * (x - d > w[0] * t)
            * (x - d <= w[1] * t)
        )
        middle_left = (
            mass_density_middle_left * (x - d > w[1] * t) * (x - d <= w[2] * t)
        )
        middle_right = (
            mass_density_middle_right * (x - d > w[2] * t) * (x - d <= w[3] * t)
        )
        right_rarefaction = (
            mass_density_rarefaction_right(x, t)
            * (x - d > w[3] * t)
            * (x - d <= w[4] * t)
        )
        right = mass_density_right * (x - d > w[4] * t)

        return (
            left
            + left_rarefaction
            + middle_left
            + middle_right
            + right_rarefaction
            + right
        )

    def velocity_solution(x, t):
        left = velocity_left * (x - d <= w[0] * t)
        left_rarefaction = (
            velocity_rarefaction_left(x, t) * (x - d > w[0] * t) * (x - d <= w[1] * t)
        )
        middle = velocity_middle * (x - d > w[1] * t) * (x - d <= w[3] * t)
        right_rarefaction = (
            velocity_rarefaction_right(x, t) * (x - d > w[3] * t) * (x - d <= w[4] * t)
        )
        right = velocity_right * (x - d > w[4] * t)
        return left + left_rarefaction + middle + right_rarefaction + right

    def pressure_solution(x, t):
        left = pressure_left * (x - d <= w[0] * t)
        left_rarefaction = (
            pressure_rarefaction_left(x, t) * (x - d > w[0] * t) * (x - d <= w[1] * t)
        )
        middle = pressure_middle * (x - d > w[1] * t) * (x - d <= w[3] * t)
        right_rarefaction = (
            pressure_rarefaction_right(x, t) * (x - d > w[3] * t) * (x - d <= w[4] * t)
        )
        right = pressure_right * (x - d > w[4] * t)
        return left + left_rarefaction + middle + right_rarefaction + right

    return mass_density_solution, velocity_solution, pressure_solution
