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


np.random.seed(1)
#x = np.exp(np.random.random((5,201)))

x = unpickle('cifar-10-batches-py/test_batch')[b'data'][:150]

t = time.time()
exact = iffn(x, npartitions=1, exact = True, seed=1)
print(time.time()-t)
t = time.time()
exact2 = iffn(x, npartitions=1, exact = True, brute = True, seed=1)
print(time.time()-t, np.array_equal(exact,exact2))
t = time.time()
exact3 = iffn(x, npartitions=16, exact = True, seed=1)
print(time.time()-t, np.array_equal(exact,exact3))
t = time.time()
exact4 = iffn(x, npartitions=16, exact = True, brute = True, seed=1)
print(time.time()-t, np.array_equal(exact,exact4))
t = time.time()
approx1 = iffn(x, npartitions=1, exact = False)
print(time.time()-t)
t = time.time()
approx16 = iffn(x, npartitions=16, exact = False)
print(time.time()-t)
t = time.time()
approx128 = iffn(x, npartitions=128, exact = False)
print(time.time()-t)

for y in approx1, approx16, approx128:
    plt.scatter(largest(x,exact), largest(x,y))
plt.show()
