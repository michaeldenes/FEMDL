""" FLOW MAP """
import sys
sys.path.append('../../src')
from src.integrator import integrator
from numpy import conj, transpose

def flow_map(v, x, tspan):

    x0 = x[0,:]
    y0 = x[1,:]

    [Fx,Fy] = integrator(v,x0,y0,tspan)

    y =  np.array((conj(transpose(Fx)), conj(transpose(Fy))))

    return y
