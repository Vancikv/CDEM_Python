'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *

dom = get_2D_quadrangle_domain(ni=3, nj=10, lx=100., ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., alfaC=100., stiffdef='isoparametric',
                             c_n=50000., c_s=50000.,
                             load=[((75., 105.), (25., 35.), [0., -50.])],
                             supports=[((-5., 5.), (-5., 35.), [1, 1])])

dom.solve(dt=0.01, tol=0.00000001, maxiter=100, node_to_plot=119)
dom.plot(magnitude=500.)
