import numpy as np

def norm(x):
    norm = np.linalg.norm(x)
    if (norm != 0):
        return x.copy()/norm
    return np.zeros(np.shape(x))