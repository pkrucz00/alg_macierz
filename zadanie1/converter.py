import numpy as np

print(np.loadtxt(open("riga.csv", "rb"), delimiter=",", skiprows=1))
