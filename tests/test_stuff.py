# MIT License
#
# Copyright (c) 2017 Lars Andersen Bratholm
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function

from fns.fn import iffn, brutefn
import numpy as np

def unpickle(file):
    with open(file, 'rb') as fo:
        try:
            import cPickle as pickle
            dict = pickle.load(fo)
        except:
            import pickle
            dict = pickle.load(fo, encoding='bytes')
    return dict

def tests():
    x = unpickle('test_batch')[b'data'][:50]

    brute = brutefn(x, seed=1)
    exact1 = iffn(x, npartitions=1, seed=1)
    exact16 = iffn(x, npartitions=16, seed=1)

    assert(np.array_equal(brute, exact1))
    assert(np.array_equal(brute, exact16))


if __name__ == "__main__":
    tests()

