import numpy as np
import matplotlib.pyplot as plt

def spy(matrix):
    mask = matrix == 0
    plt.matshow(mask)
    plt.show()