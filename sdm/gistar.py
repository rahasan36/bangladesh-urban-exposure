# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import numpy as np
from scipy.ndimage import uniform_filter
def gistar(x, m, k=3):
    x = np.where(m, x, 0.0).astype(float); n = m.sum(); xb = x[m].mean(); s = np.sqrt((x[m] ** 2).mean() - xb ** 2)
    sw = uniform_filter(x, k, mode="constant") * k * k; wi = uniform_filter(m.astype(float), k, mode="constant") * k * k
    with np.errstate(all="ignore"): zz = (sw - xb * wi) / (s * np.sqrt((n * wi - wi ** 2) / (n - 1)))
    return np.where(m, zz, np.nan)
