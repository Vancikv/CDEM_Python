'''
Created on 10. 11. 2015

@author: Kapsak

'''

from CDEM import *
import numpy as np

dom = get_2D_quadrangle_domain(ni=3, nj=20, lx=2., ly=0.2, E=25000000000., nu=0.2, density=2500.,
                             thickness=0.2, stiffdef='isoparametric',
                             spring_stiffness_factor=10., domtype='CDEMStatic',
                             load=[((0.99, 1.01), (0.19, 0.21), [0., -200000.])],
                             supports=[((-0.01, 0.01), (-0.01, 0.01), [1, 1]), ((1.99, 2.01), (-0.01, 0.01), [0, 1])])

dom.solve()
dom.plot(magnitude=10.)
