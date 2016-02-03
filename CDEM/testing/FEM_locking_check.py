'''
Created on 22. 11. 2015

@author: Kapsak
'''

from CDEM.basic import *
from CDEM.utils import *
import numpy as np

ni, nj = 3, 10
dom = get_2D_quadrangle_FEM_domain(ni=ni, nj=nj, lx=50., ly=10., E=10000., nu=0.2, density=1.,
                             thickness=1., stiffdef='isoparametric',
                             load=[((49., 51.), (-1., 11.), [10., 0.])],
                             supports=[((-1., 1.), (-1., 11.), [1, 1])])

dom.solve()
dom.plot(magnitude=5.)
print 'Reference node coordinates: x = %f, y = %f' % (dom.nodes[(nj + 1) * (ni + 1) - 1].x, dom.nodes[(nj + 1) * (ni + 1) - 1].y)
print 'FEM deflection: %.6f' % dom.nodes[(nj + 1) * (ni + 1) - 1].v_disp[1]
v_disp = np.hstack([n.v_disp for n in dom.nodes])
print 'FEM L2 norm: %.6f' % np.linalg.norm(v_disp)
print dom.elements[10].get_stress()