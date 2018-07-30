"""Integreator Function """
from scipy import integrate
from numpy import concatenate

def integrator(v, x0, y0, tspan):
    #Ps = integrate.odeint(rotating_double_gyre_vf, np.concatenate((x0,y0)), tspan, tfirst=True)
    tminmax = np.array((tspan[0],tspan[-1]))
    Ps = integrate.solve_ivp(v, tminmax, concatenate((x0,y0)), method='RK45', vectorized=True, atol=1e-3, rtol=1e-3, t_eval=tspan)
    print('Integrator: ' + Ps.message)
    F = Ps.y

    xt = F[0:len(x0),:]
    yt = F[len(x0):,:]

    return xt, yt
