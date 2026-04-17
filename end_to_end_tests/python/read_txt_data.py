import numpy as np

def read_txt_data(filename):
    with open(filename) as f:
        header = f.readline()
        stripped_header = header.strip('#\n ')
        labels = stripped_header.split(' ')
        data = np.loadtxt(filename)
        result = {}
        for i, label in enumerate(labels):
            result[label] = data[:, i]

    return result