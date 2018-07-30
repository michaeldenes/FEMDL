""" ASSEMBLE """
# Add path to FEMDL code
import sys
sys.path.append('../../src')
from src.gradbasis import gradbasis
from numpy import max,concatenate, ones
from scipy.sparse import csc_matrix
def assemble(p,t,pb=None,G=None):
    n = max(pb[1,:]) + 1
    m = len(t)

    dphi, area = gradbasis(p,t)

    D = csc_matrix((n, n))
    M = csc_matrix((n, n))


    for i in range(0,3):
        for j in range(0,3):
            #Is the formula for Dij wrong?
            Dij = -area*(dphi[0,i,:]*G[:,0]*dphi[0,j,:] + dphi[0,i,:]*G[:,1]*dphi[1,j,:] + dphi[1,i,:]*G[:,1]*dphi[0,j,:] + dphi[1,i,:]*G[:,2]*dphi[1,j,:])
            Mij = area/12*np.ones(np.shape(dphi[0,i]))
            I = pb[1, t[:,i]]
            J = pb[1, t[:,j]]

            if (i == j):
                D = D + csc_matrix((Dij, (I,J)), shape=(n,n))
                M = M + csc_matrix((Mij+ area/12, (I,J)), shape=(n,n))

            else:
                D = D + csc_matrix((concatenate((Dij,Dij)), (concatenate((I,J)), concatenate((J,I)))), shape=(n,n))
                M = M + csc_matrix((concatenate((Mij,Mij)), (concatenate((I,J)), concatenate((J,I)))), shape=(n,n))

    return [D, M]
