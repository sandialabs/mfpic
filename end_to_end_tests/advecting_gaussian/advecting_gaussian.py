import utils

import subprocess

domain_length = 2.0
base_num_elements = 50

temperature = 300
velocity = 3.0
max_cfl = 0.9

final_time = 0.5 * domain_length / velocity


def get_input_deck(refinement_level):
    num_elements = base_num_elements * refinement_level
    dx = domain_length / num_elements
    dt, num_time_steps = utils.compute_timestepping_that_satisfies_cfl(max_cfl, dx, velocity, final_time)

    input_deck_contents = f"""
Mesh:
  Type: line
  Lengths: [{domain_length}]
  Number of Elements: [{num_elements}]
  Periodic Dimensions: [x]

Time Stepping:
  Number of Time Steps: {num_time_steps}
  Time Step Size: {dt}
"""
    return input_deck_contents


def run(mfpic_executable):
    input_deck_contents = get_input_deck(1)
    yaml = "advecting_gaussian.yaml"
    with open(yaml, "w") as input_deck:
        input_deck.write(input_deck_contents)

    result = subprocess.run([mfpic_executable, "-i", yaml])
    result.check_returncode()


def analyze():
    pass


def plot():
    pass


if __name__ == "__main__":
    import sys

    if "run" in sys.argv[1:]:
        run(sys.argv[2])
    elif "plot" in sys.argv[1:]:
        plot()
    else:
        analyze()
