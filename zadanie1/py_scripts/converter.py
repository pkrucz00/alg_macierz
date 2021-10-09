import numpy as np


def converter():
    return np.loadtxt(open("riga.csv", "rb"), delimiter=",", skiprows=1)