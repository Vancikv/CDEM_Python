'''
Created on 7. 11. 2015

@author: Kapsak
'''

import numpy as np

GRAVITY = 0.

def get_plane_stress_stiffness_matrix(E, nu):
    shear = E / (2 + 2 * nu)
    lame = E * nu / ((1 + nu) * (1 - 2 * nu))
    C = np.array([[2 * shear + lame, lame, 0.],
                  [lame, 2 * shear + lame, 0.],
                  [0., 0., shear]])
    return C

def get_triangle_stiffness_matrix(C, x, y, t):
        rot = lambda l, n: l[n:] + l[:n]
        ids = [0, 1, 2]
        a = [x[j] * y[k] - x[k] * y[j] for i, j, k in [rot(ids, n) for n in ids]]
        b = [y[j] - y[k] for i, j, k in [rot(ids, n) for n in ids]]
        c = [x[k] - x[j] for i, j, k in [rot(ids, n) for n in ids]]
        A = 0.5 * np.linalg.det(np.transpose(np.vstack((np.ones(3), x, y))))

        # du/dx, du/dy, B matrix
        b1 = 1 / 2. / A * np.array(b)
        b2 = 1 / 2. / A * np.array(c)
        B = np.vstack((np.hstack((b1, np.zeros(3))), np.hstack((np.zeros(3), b2)), np.hstack((b2, b1))))

        K = A * t * np.dot(np.dot(np.transpose(B), C) , B)
        K = K[[0, 3, 1, 4, 2, 5], :]
        K = K[:, [0, 3, 1, 4, 2, 5]]
        return (K, A)

A, B = .8, .5
c1 = (B - A * A) / (A * A * A - A * A)
c2 = 1. - c1

def load_function(tau):
    '''
    Gets the relative time <0,1>, returns the load coefficient <0,1>
    '''
    return c1 * (tau ** 3.) + c2 * (tau ** 2.) if tau < 1. else 1.
    # return tau ** 3 - 3 * tau ** 2 + 3 * tau
    # return 1.

if __name__ == '__main__':
    pass
