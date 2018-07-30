""" norm """
from numpy.linalg import norm
from numpy import inf
def normed(x):
    y = x/norm(x,inf)
    return y
