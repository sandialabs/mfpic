import sys

sys.path.append("../../python")
sys.path.append("..")

import advecting_gaussian

if __name__ == "__main__":
    import sys

    time_integrator = "Verlet"
    if "run" in sys.argv[1:]:
        mfpic_executable = sys.argv[2]
        advecting_gaussian.run(mfpic_executable, time_integrator)
    elif "plot" in sys.argv[1:]:
        advecting_gaussian.plot()
    else:
        advecting_gaussian.analyze()