from fns.fn import iffn, brutefn, fastfn

import numpy as np
import matplotlib.pyplot as plt
import time
from collections import Counter

def unpickle(file):
    import pickle
    with open(file, 'rb') as fo:
        dict = pickle.load(fo, encoding='bytes')
    return dict

def largest(x, order, metric = "l1"):
    d = np.zeros(order.shape)
    if metric == "l1":
        dist_mat = np.sum(abs(x[:,None,:] - x[None,:,:]), axis=2)
    else:
        dist_mat = np.sum((x[None,:,:] - x[:,None,:])**2, axis=2)
    for i in range(1, order.size):
        d[i] = np.max(dist_mat[np.ix_(order[:i]-1, order[i:i+1]-1)])
    return d[1:]

n = 200
x = unpickle('../tests/test_batch')[b'data'][:n]
x = np.asarray(x, dtype=float)

dist_mat_l1 = np.sum(abs(x[:,None,:] - x[None,:,:]), axis=2)
dist_mat_l2 = np.sum((x[None,:,:] - x[:,None,:])**2, axis=2)

def get_error(x, dist_mat, fn):
    error = []
    for seed in range(1, x.shape[0]+1):
        order = eval(fn)[1]
        order_dist = dist_mat[seed-1, order]
        dist = dist_mat[seed-1]
        maxdist =  np.max(dist)
        error.append(sum((x > order_dist) for x in dist))

    c = Counter(error)
    return c, np.mean(error)

#print("Mean error in order for %d points for l1 approximation:" %n)
#for i in range(0,9):
#    npart = 2**i
#    fn = "fastfn(x, seed=seed, nmax = 2, npartitions=%d, metric='l1')" % npart
#    #fn = "brutefn(x, seed=seed, nmax = 2, metric='l1')"
#    c,e = get_error(x, dist_mat_l1, fn)
#    print("{:s} {:3d} {:s} {:.4f}".format("Mean error with", npart, "partitions:", e))
#print("\nMean error in order for %d points for l2 approximation:" %n)
#for i in range(0,9):
#    npart = 2**i
#    fn = "fastfn(x, seed=seed, nmax = 2, npartitions=%d, metric='l2')" % npart
#    c,e = get_error(x, dist_mat_l2, fn)
#    print("{:s} {:3d} {:s} {:.4f}".format("Mean error with", npart, "partitions:", e))

t = time.time()
brutefn(x, metric = "l1", seed=1)
print("{:20s}: {:6.3} seconds".format("\nTime spent on l1 brute force furthest neighbours", time.time()-t))
t = time.time()
fastfn(x, metric = "l1", seed=1, npartitions=1)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 approximation by   1 partition ", time.time()-t))
t = time.time()
fastfn(x, metric = "l1", seed=1, npartitions=16)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 approximation by  16 partitions", time.time()-t))
t = time.time()
fastfn(x, metric = "l1", seed=1, npartitions=64)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 approximation by  64 partitions", time.time()-t))
t = time.time()
fastfn(x, metric = "l1", seed=1, npartitions=256)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 approximation by 256 partitions", time.time()-t))
t = time.time()
iffn(x, metric = "l1", seed=1, npartitions=1)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 exact sped up by 1 partitions", time.time()-t))
t = time.time()
iffn(x, metric = "l1", seed=1, npartitions=16)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 exact sped up by 16 partitions", time.time()-t))
t = time.time()
iffn(x, metric = "l1", seed=1, npartitions=64)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 exact sped up by 64 partitions", time.time()-t))
t = time.time()
iffn(x, metric = "l1", seed=1, npartitions=256)
print("{:20s}: {:1.4} seconds".format("Time spent on l1 exact sped up by 256 partitions", time.time()-t))


t = time.time()
brutefn(x, metric = "l2", seed=1)
print("{:20s}: {:6.3} seconds".format("\nTime spent on l2 brute force furthest neighbours", time.time()-t))
t = time.time()
iffn(x, seed=1, npartitions=1)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 exact sped up by 1 partitions", time.time()-t))
t = time.time()
iffn(x, seed=1, npartitions=16)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 exact sped up by 16 partitions", time.time()-t))
t = time.time()
iffn(x, seed=1, npartitions=64)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 exact sped up by 64 partitions", time.time()-t))
t = time.time()
iffn(x, seed=1, npartitions=256)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 exact sped up by 256 partitions", time.time()-t))
t = time.time()
fastfn(x, metric = "l2", seed=1, npartitions=1)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 approximation by   1 partition ", time.time()-t))
t = time.time()
fastfn(x, metric = "l2", seed=1, npartitions=16)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 approximation by  16 partitions", time.time()-t))
t = time.time()
fastfn(x, metric = "l2", seed=1, npartitions=64)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 approximation by  64 partitions", time.time()-t))
t = time.time()
fastfn(x, metric = "l2", seed=1, npartitions=256)
print("{:20s}: {:1.4} seconds".format("Time spent on l2 approximation by 256 partitions", time.time()-t))

