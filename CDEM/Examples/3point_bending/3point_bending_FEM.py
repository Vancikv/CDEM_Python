'''
Created on 10. 11. 2015

@author: Kapsak
'''

from CDEM import *
import numpy as np

dom = get_2D_quadrangle_FEM_domain(ni=3, nj=30, lx=100., ly=10., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             load=[((49., 51.), (9., 11.), [0., -20.])],
                             supports=[((-1., 1.), (-1., 1.), [1, 1]), ((99., 101.), (-1., 1.), [0, 1])])

dom.solve()
dom.plot(magnitude=5.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[62].x, dom.nodes[47].y)
print 'FEM deflection: %.6f' % dom.nodes[62].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'FEM L2 norm: %.6f' % np.linalg.norm(v_disp)
print dom.get_stress_norm(ids=[0])
