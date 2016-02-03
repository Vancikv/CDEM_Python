'''
Created on Jan 4, 2016

@author: vancik
'''

import numpy as np
from basic import load_function, GRAVITY

class Node(object):
    def __init__(self, x, y, F_ext=[0., 0.], neighbors=[0, 0], supports=[0, 0], nnodedofs=2):
        self.x = x
        self.y = y
        self.F_ext = np.array(F_ext)
        self.neighbors = neighbors
        self.n_vects = [0, 0]
        self.l_edges = [0, 0]
        self.free_dofs = [i for i, j in enumerate(supports) if j == 0]
        self.supports = supports
        self.nnodedofs = nnodedofs

        self.v_disp = np.array([0., 0.])
        self.v_velo = np.array([0., 0.])
        self.v_acce = np.array([0., 0.])

    def init_vals(self, maxiter):
        # Init for dynamic relaxation
        tau = 1 / float(maxiter)
        self.F_g = load_function(tau) * np.array([0., -self.mass * GRAVITY])
        self.v_acce[self.free_dofs] = load_function(tau) / self.mass * (self.F_ext + self.F_g)

    def set_codes(self, maxcode, verbose=False):
        self.v_code = []
        for s, i in zip(self.supports, range(self.nnodedofs)):
            if s == 0:
                maxcode += 1
                self.v_code.append(maxcode)
            else:
                self.v_code.append(0)
        if verbose: print self.v_code
        return maxcode

    def set_disp(self, val):
        self.v_disp[self.free_dofs] = val[self.free_dofs]

    def set_velo(self, val):
        self.v_velo[self.free_dofs] = val[self.free_dofs]

    def set_acce(self, val):
        self.v_acce[self.free_dofs] = val[self.free_dofs]
