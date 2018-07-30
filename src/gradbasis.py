""" GRAD BASIS """

# Based on below:
#
# [dphi,area] = GRADBASIS(p,t) computes the gradients of the shape
# functions on the standard simplex and the areas of the triangles
#   p: (n x 2), one node per row
#   t: (m x 3), integers, each row defines a triangle by indexing into p
#
# based on code from ifem by Long Chen
#
# (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT

from numpy import asarray, array
def gradbasis(p,t):

    v = []
    #inserting with positions rather than append to make sure order is correct
    v.insert(0, p[:,t[:,2]] - p[:,t[:,1]])
    v.insert(1, p[:,t[:,0]] - p[:,t[:,2]])
    v.insert(2, p[:,t[:,1]] - p[:,t[:,0]])

    v = np.asarray(v) #make as np array so indexing is normal

    area = 0.5*(-v[2,0,:]*v[1,1,:] + v[2,1,:]*v[1,0,:])

    dphi = []
    dphi.insert(0, -v[:,1,:]/(2*area))
    dphi.insert(1, v[:,0,:]/(2*area))
    dphi = np.asarray(dphi)

    area = abs(area)

    return dphi, area
