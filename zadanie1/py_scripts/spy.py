import numpy as np
import matplotlib.pyplot as plt

def spy(matrix):
    mask = matrix == 0
    plt.matshow(mask)
    plt.show()


spy(np.array([[0,1,1], [0,1,1], [1,0,1]]))