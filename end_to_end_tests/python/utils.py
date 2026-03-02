import numpy as np


def compute_timestepping_that_satisfies_cfl(max_cfl, dx, max_wavespeed, final_time):
    """Compute the timestep size and number of timesteps that will get to final_time
    while satisfying the cfl condition given by max_cfl, dx, and max_wavespeed"""
    max_dt = max_cfl * dx / max_wavespeed

    num_time_steps = int(np.ceil(final_time / max_dt))
    dt = final_time / num_time_steps

    return dt, num_time_steps
