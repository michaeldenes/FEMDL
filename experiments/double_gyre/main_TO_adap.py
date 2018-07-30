""" Double Gyre Matlab to Python conversion
- MATLAB code from FEMDL Matlab
"""
# Add path to FEMDL code
import sys
sys.path.append('../../src')

# Import Statements
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from experiments.double_gyre.double_gyre import rotating_double_gyre_vf
from src.normed import normed
from src.flow_map import flow_map
from src.integrator import integrator

# Flow map
t0 = 0; tf = 1; nt = 6; tspan = np.linspace(t0, tf, nt)

# Data points
n = 25; x = np.linspace(0,1,n)
[X,Y] = np.meshgrid(x,x, indexing = 'ij')
p = [] #think about np arrays, this is a list! hence I think we're getting confused...
p_0 = np.concatenate((np.concatenate((X), axis=None), np.concatenate((Y), axis=None)), axis=None)
p.append(np.reshape(p_0, (-1,np.square(n))))

pb = np.array([np.arange(1,np.square(n)+1), np.arange(1,np.square(n)+1)])

# Time Integration
def T(x):
    return flow_map(rotating_double_gyre_vf, x, tspan)

P = T(p[0])
for k in range(1,nt):
    p_val = np.concatenate((P[0][k][:],  P[1][k][:]), axis=None)
    p.append(np.reshape(p_val, (-1,np.square(n))))

# Assembly
pm = 1
D = sp.sparse.csr_matrix((np.square(n), np.square(n)))
M = sp.sparse.csr_matrix((np.square(n), np.square(n)))

for k in range(0, nt):
    r = np.random.permutation(np.arange(0, np.square(n)))[0:int(np.floor(pm*np.square(n)))]
    tr = sp.spatial.Delaunay(np.transpose(p[k][:,r]))
    t = [r[tr.simplices[:,0]], r[tr.simplices[:,1]], r[tr.simplices[:,2]]]
    A = np.kron([1, 0, 1], np.ones((len(t[0]), 1)))





#Look at Shane/Erics code, there are a number of numpy functions that you could use that you are not currently using...
