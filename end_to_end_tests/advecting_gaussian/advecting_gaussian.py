
def run(mfpic_executable):
    pass

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