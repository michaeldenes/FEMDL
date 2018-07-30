""" DOUBLE GYRE - MATLAB TO PYTHON CONVERSION"""

from numpy import square, pi, cos, sin, shape, zeros_like
def rotating_double_gyre_vf(t,z):
    n = int(len(z)/2)
    x = z[0:n]
    y = z[n:]

    st = (t > 0 and t < 1)*square(t)*(3-2*t) + (t > 1)*1

    dxPsi_P = 2*pi*cos(2*pi*x)*sin(pi*y)
    dyPsi_P = pi*sin(2*pi*x)*cos(pi*y)
    dxPsi_F = pi*cos(pi*x)*sin(2*pi*y)
    dyPsi_F = 2*pi*sin(pi*x)*cos(2*pi*y)
    dxPsi = (1-st)*dxPsi_P + st*dxPsi_F
    dyPsi = (1-st)*dyPsi_P + st*dyPsi_F
    dz = zeros_like(z)
    dz[0:n]= -dyPsi
    dz[n:] = dxPsi

    return dz
