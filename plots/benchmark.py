from fns.fn import iffn
import numpy as np
import matplotlib.pyplot as plt
import qml.distance
import time

def unpickle(file):
    import pickle
    with open(file, 'rb') as fo:
        dict = pickle.load(fo, encoding='bytes')
    return dict

def largest(x, order):
    d = np.zeros(order.shape)
    dist_mat = qml.distance.l2_distance(x,x)
    for i in range(1, order.size):
        d[i] = np.max(dist_mat[np.ix_(order[:i]-1, order[i:i+1]-1)])
    return d[1:]


x = unpickle('cifar-10-batches-py/test_batch')[b'data'][:200]

t = time.time()
exact = iffn(x, npartitions=1, exact = True, seed=1)
print(time.time()-t)
t = time.time()
exact2 = iffn(x, npartitions=1, exact = True, brute = True, seed=1)
print(time.time()-t)
t = time.time()
exact3 = iffn(x, npartitions=16, exact = True, seed=1)
print(time.time()-t)
t = time.time()
exact4 = iffn(x, npartitions=16, exact = True, brute = True, seed=1)
print(time.time()-t)
t = time.time()
exact5 = iffn(x, npartitions=128, exact = True, seed=1)
print(time.time()-t)
t = time.time()
exact4 = iffn(x, npartitions=128, exact = True, brute = True, seed=1)
print(time.time()-t)
t = time.time()
