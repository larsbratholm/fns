from fns.fn import iffn
import numpy as np


x = np.random.random((5,21))

print(iffn(x, npartitions=2))
