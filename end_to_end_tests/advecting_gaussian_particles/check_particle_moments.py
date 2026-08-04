import h5py
import numpy as np

particle_data = h5py.File("particles.h5part")
particle_moments_from_file = np.genfromtxt('particle_moments_N2.csv', names=True, delimiter=',')

num_steps = int(np.max(particle_moments_from_file['step']) + 1)

temp = dict()
for i in range(num_steps):
    mask = particle_moments_from_file['step'] == i
    temp[i] = particle_moments_from_file[mask]
particle_moments_from_file = temp

num_cells = int(np.max(particle_moments_from_file[0]['elem']) + 1)
cell_size = particle_moments_from_file[0]['x'][1] - particle_moments_from_file[0]['x'][0]
print(cell_size)

particle_moments_from_data = dict()
for i in range(num_steps):
    sum_of_weights_in_cell = np.zeros(num_cells)
    particle_data_i = particle_data[f"Step#{i}"]
    num_particles = int(particle_data_i['element'].shape[0])
    for i_particle in range(num_particles):
        particle_weight = particle_data_i['weight'][i_particle]
        i_cell = particle_data_i['element'][i_particle]
        sum_of_weights_in_cell[i_cell] += particle_weight

    number_density = sum_of_weights_in_cell / cell_size
    particle_moments_from_data[i] = dict()
    particle_moments_from_data[i]['number_density'] = number_density

for i in range(num_steps):
    number_density_difference = particle_moments_from_file[0]['number_density'] - particle_moments_from_data[0]['number_density']
    relative_error = np.linalg.norm(number_density_difference) / np.linalg.norm(particle_moments_from_data[0]['number_density'])
    print(relative_error)