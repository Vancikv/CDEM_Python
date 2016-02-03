'''
Created on Nov 25, 2015

@author: Werner
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np

dom2 = get_2D_quadrangle_domain(ni=3, nj=40, lx=400.,ly=30., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             c_n=100000.,c_s=100000., domtype='CDEMStatic',
                             load=[((199., 201.), (29., 31.), [0., -50.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((399., 401.), (-1., 1.), [0, 1])])
dom2.solve()
dom2.plot(magnitude=10.)
print 'CDEM static deflection: %.6f, CDEM/FEM = \n' % (dom2.nodes[238].v_disp[1])