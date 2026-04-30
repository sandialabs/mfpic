import sys

sys.path.append("../../python")
sys.path.append("..")

import sod_shock_tube

if __name__ == "__main__":
    import sys

    time_integrator = "Verlet"
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        sod_shock_tube.run(mfpic_executable, time_integrator)
    elif "plot" in sys.argv[1:]:
        sod_shock_tube.plot()
    else:
        sod_shock_tube.analyze()