import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.special import erf
import scipy.stats

def get_distance_matrix(x):
    d = np.sum(((x[:,None,:] - x[None,:,:]))**2,axis=2)
    return d

def get_distance_matrix_1(x):
    means = x.mean(1)
    std = x.std(1)
    d = x.shape[1] * ((means[:,None] - means[None,:])**2 + std[None,:]**2 + std[:,None]**2)
    return d

def get_distance_matrix_2(X,n):
    x = X - X.mean(0)
    means = np.zeros((x.shape[0],n))
    std = np.zeros((x.shape[0],n))
    bin_size = x.shape[1] // n
    if bin_size * n != x.shape[1]:
        bin_size += 1
    d = np.zeros((x.shape[0],)*2)
    for i in range(n):
        start = i * bin_size
        end = (i+1) * bin_size
        means[:,i] = x[:,start:end].mean(1)
        std[:,i] = x[:,start:end].std(1)
        d += x[0,start:end].size * ((means[:,None,i] - means[None,:,i])**2 + std[None,:,i]**2 + std[:,None,i]**2)
    return d

def get_distance_matrix_3(y,n):
    m = y.mean(0)
    s = y.std(0)
    order = np.argsort(m**2+s**2)
    x = y[:,order]
    means = np.zeros((x.shape[0],n))
    std = np.zeros((x.shape[0],n))
    bin_size = x.shape[1] // n
    if bin_size * n != x.shape[1]:
        bin_size += 1
    d = np.zeros((x.shape[0],)*2)
    for i in range(n):
        start = i * bin_size
        end = (i+1) * bin_size
        means[:,i] = x[:,start:end].mean(1)
        std[:,i] = x[:,start:end].std(1)
        d += x[0,start:end].size * ((means[:,None,i] - means[None,:,i])**2 + std[None,:,i]**2 + std[:,None,i]**2)
    return d

def get_mdistance_matrix(x):
    d = np.sum(abs((x[:,None,:] - x[None,:,:])),axis=2)
    return d

def get_mdistance_matrix_1(x):
    means = x.mean(1)
    std = x.std(1)
    d = np.zeros((x.shape[0],)*2)
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            mu = means[i] - means[j]
            s2 = std[i]**2 + std[j]**2
            d[i,j] = x.shape[1] * (np.sqrt(2/np.pi) * s2**0.5 * np.exp(-0.5*mu**2 / s2) + abs(mu)*erf(abs(mu)/2**0.5 / s2**0.5))
    return d

def get_mdistance_matrix_2(x,n):
    means = np.zeros((x.shape[0],n))
    std = np.zeros((x.shape[0],n))
    bin_size = x.shape[1] // n
    if bin_size * n != x.shape[1]:
        bin_size += 1
    d = np.zeros((x.shape[0],)*2)
    for k in range(n):
        start = k * bin_size
        end = (k+1) * bin_size
        means[:,k] = x[:,start:end].mean(1)
        std[:,k] = x[:,start:end].std(1)
        for i in range(x.shape[0]):
            for j in range(x.shape[0]):
                mu = means[i,k] - means[j,k]
                s2 = std[i,k]**2 + std[j,k]**2
                d[i,j] += x[0,start:end].size  * (np.sqrt(2/np.pi) * s2**0.5 * np.exp(-0.5*mu**2 / s2) + abs(mu)*erf(abs(mu)/2**0.5 / s2**0.5))
    return d

def get_mdistance_matrix_3(x,n,m):
    means = np.zeros((x.shape[0],n))
    std = np.zeros((x.shape[0],n))
    bin_size = x.shape[1] // n
    if bin_size * n != x.shape[1]:
        bin_size += 1
    d = np.zeros((x.shape[0],)*2)
    for k in range(n):
        start = k * bin_size
        end = (k+1) * bin_size
        means[:,k] = x[:,start:end].mean(1)
        std[:,k] = x[:,start:end].std(1)
        for i in range(x.shape[0]):
            for j in range(x.shape[0]):
                mu = means[i,k] - means[j,k]
                s2 = std[i,k]**2 + std[j,k]**2
                d[i,j] += x[0,start:end].size  * (np.sqrt(2/np.pi) * s2**0.5 * np.exp(-0.5*mu**2 / s2) + abs(mu)*fast_erf(abs(mu)/2**0.5 / s2**0.5))
    return d

def get_mdistance_matrix_4(y,n):
    m = y.mean(0)
    s = y.std(0)
    order = np.argsort(m*s)
    x = y[:,order]
    means = np.zeros((x.shape[0],n))
    std = np.zeros((x.shape[0],n))
    bin_size = x.shape[1] // n
    if bin_size * n != x.shape[1]:
        bin_size += 1
    d = np.zeros((x.shape[0],)*2)
    for k in range(n):
        start = k * bin_size
        end = (k+1) * bin_size
        means[:,k] = x[:,start:end].mean(1)
        std[:,k] = x[:,start:end].std(1)
        for i in range(x.shape[0]):
            for j in range(x.shape[0]):
                mu = means[i,k] - means[j,k]
                s2 = std[i,k]**2 + std[j,k]**2
                d[i,j] += x[0,start:end].size  * (np.sqrt(2/np.pi) * s2**0.5 * np.exp(-0.5*mu**2 / s2) + abs(mu)*erf(abs(mu)/2**0.5 / s2**0.5))
    return d

def fast_exp(x, n):
    if x < -2**n:
        return 0
    for i in range(n):
        x /= 2
    y = 1 + x
    for i in range(n-1):
        y *= y
    return y

def faster_erf(x):
    return 1 - np.exp(-16/23 * x**2 - 2/np.sqrt(np.pi) * x)

def fast_erf(x):
    w1 = 928/1175
    a = 653/720
    b1 = 4532/6043
    w2 = (1-w1)
    b2 = (2/np.sqrt(np.pi) - b1*w1) / w2
    return 1 - w1 * np.exp(-a*x**2-b1*x) - w2 * np.exp(-a*x**2-b2*x)

def fastest_erf(x,n):
    return 1 - fast_exp(-16/23 * x**2 - 2/np.sqrt(np.pi) * x,n)

def unpickle(file):
    with open(file, 'rb') as fo:
        dict = pickle.load(fo, encoding='bytes')
    return dict

x = unpickle('../tests/test_batch')[b'data'][:100].astype(np.float)

idx = np.triu_indices(x.shape[0],1)

#print(get_distance_matrix(x)[0,1])
#print(get_distance_matrix_1(x)[0,1])


z = get_distance_matrix(x)[idx]
y = get_distance_matrix_2(x,2)[idx]
print(scipy.stats.pearsonr(z,y))
y = get_distance_matrix_2(x,8)[idx]
print(scipy.stats.pearsonr(z,y))
y = get_distance_matrix_2(x,32)[idx]
print(scipy.stats.pearsonr(z,y))
plt.scatter(z, y)
plt.show()
